#!/usr/bin/python -u
# -*- coding: utf-8 -*-
##### nohup ./infiniCharges.py > output.dat &

################################################################################
#
#    InfiniCharges: A program to generate partial charges for periodic systems
#    Copyright (C) 2015  Andrea Gabrieli and Marco Sant                        
# 
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
# 
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor,
#    Boston, MA  02110-1301, USA.
# 
#    See also http://www.fsf.org/licensing/licenses/gpl.html
# 
#    InfiniCharges, including its sources and pointers to the authors
#    can be found at http://www.physchem.uniss.it/cpc/
#
#    Contact Address:
#    agabrieli@uniss.it
#    msant@uniss.it
#
################################################################################

################################################################################
#
#    InfiniCharges is a computer program for generating reliable
#    partial charges for molecular simulations in periodic systems.
#    It relies on the DM-REPEAT method where the stability of the
#    resulting charges, over a large set of fitting regions, is obtained
#    through the simultaneous fit of multiple electrostatic potential (ESP)
#    configurations together with the total dipole fluctuations (TDF).
#
#    This program performs the following kinds of fits:
#    M-REPEAT (also standard REPEAT)
#    DM-REPEAT (also D-REPEAT)
#    PARABOLIC RESTRAINED M-REPEAT
#    "RAPPE-GODDARD LIKE" RESTRAINED M-REPEAT
#
################################################################################

import numpy as np
from numpy.linalg import lstsq
import ctypes
from os import getcwd, getpid, path, makedirs
from time import time
from datetime import datetime, timedelta
from collections import Counter
from dir_atautils.atautils import atautils
from settings import settings, Tiny, ata_prefix

print "________     ___________      __________________                                     "
print "____  _/________  __/_(_)________(_)_⚡ ____/_/ /________ ______________ _____________"
print " __  /__  __ \_  /_ _  /_  __ \_  /_/ /   __/ __ \  ՛__ `/_/ ___/_  __ `/  _ \_  ___/"
print "__/ / _  / / /  __/_  /_  / / /  / / /___ _/ / / // /_/ /_/ / ___  /_/ //  __/(__  ) "
print "/___/ /_/ /_//_/   /_/ /_/ /_//_/  \____/ /_/ /_/ \__,_/ /_/   __\__, / \___//____/  "
print "                                                                /____/               "
print ""
print "Version 1.0"
print ""

def read_xyz_tr_data(file_unit, starting_line, end_line, names_col, x_col, y_col, z_col):
    """
    Return data taken from a snapshot of an xyz like file.

    Parameters
    ----------
    file_unit : pointer to an already opened file.
    starting_line : integer
        it is the first line in the file to be read (1 based)
    end_line : integer
        it is the last line in the file to be read (1 based)
    names_col : integer
        it is the column containing the atom types (1 based)
    x_col : integer 
        it is the column containing the x coordinates (1 based)
    y_col : integer 
        it is the column containing the y coordinates (1 based)
    z_col : integer 
        it is the column containing the z coordinates (1 based)

    Returns
    -------
    names : list of size N  = (end_line - starting_line + 1) 
        it contains the atom types
    x : list of size N  = (end_line - starting_line + 1) 
        it contains the x coordinates
    y : list of size N  = (end_line - starting_line + 1) 
        it contains the y coordinates
    z : list of size N  = (end_line - starting_line + 1) 
        it contains the z coordinates

    Examples
    --------
    Read from inf (already opened) from line 3 throug 100 (included), taking
    from column 1 the atom types and from column 2 throug 4 x, y, and z 
    coordinates, respectively.
    
    names,x,y,z = read_xyz_tr_data(inf,3,100,1,2,3,4)
    """
    names = []
    x = []
    y = []
    z = []
    
    for i in range(end_line):
        line = file_unit.readline()
        if (i >= (starting_line - 1) and i < (end_line)):
            s_line = line.split()
            names.append(s_line[names_col-1])
            x.append(float(s_line[x_col-1]))
            y.append(float(s_line[y_col-1]))
            z.append(float(s_line[z_col-1]))
    
    return names,x,y,z

def read_full_traj_tdf_auto(data):
    """
    Return the ATAs for the TDF needed to solve the linear system.

    This function requires only the names of the files to be processed
    which are contained in "data".

    Parameters
    ----------
    data : instance of class Settings
         object containing the filename of a trajectory (data.primitive_pos),
         the corresponding total dipole moments (data.primitive_tdf),
         and a logical variable telling if the latter
         comes from CP2K (data.tdf_from_cp2k).

    Returns
    -------
    atb : (N,) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and the reference data column matrix (of size M)
    ata : (N,N) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and itself (of size M,N)
    btb : float 
        dot product between the transposed reference data column matrix and 
        itself, this is the magnitude of the reference data
    n_atoms : int
        number of system atoms 
    nstep :
        number of system configurations
    
    Notes
    -----
    This is a wrapper to the read_full_traj_tdf function which 
    actually processes the files.
    """

    num_lines_p = sum(1 for line in open(data.primitive_pos,"r"))
    num_lines_b = sum(1 for line in open(data.primitive_tdf,"r"))
    n_atoms = int(open(data.primitive_pos,"r").readline())
    nstep   = num_lines_p / (n_atoms + 2) 
    
    # check consistency between trajectory and dipoles
    if data.tdf_from_cp2k: # in CP2K there are 10 lines for each dipole
        assert (num_lines_b / 10) == nstep
    else: # in .tdf format there is an header of 4 lines and then the dipoles
        assert (num_lines_b - 4) == nstep
    
    atb, ata, btb = read_full_traj_tdf(n_atoms,nstep,data)
    
    return atb, ata, btb, n_atoms, nstep


def read_full_traj_tdf(n_atoms, nstep, data):
    """
    Return the ATAs for the TDF needed to solve the linear system.

    Parameters
    ----------
    n_atoms : int
        number of system atoms 
    nstep : int
        number of system configurations
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    atb : (N,) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and the reference data column matrix (of size M)
    ata : (N,N) float ndarray
        dot product between the transposed model data matrix (of size N,M)
        and itself (of size M,N)
    btb : float 
        dot product between the transposed reference data column matrix and 
        itself, this is the magnitude of the reference data
    """
    
    in_f_pos = open(data.primitive_pos,"r")
    
    traj_p = np.zeros(nstep*3*n_atoms)
    traj_f = np.zeros(nstep*3)
    atoms_per_par = data.atoms_per_par
   
    #A_ij elements for TDF (model data)
    current_tdf = np.zeros([len(atoms_per_par),nstep*3])
    #loop over all frames
    for i in xrange(nstep):
        names_tmp,x_t,y_t,z_t = read_xyz_tr_data(in_f_pos,3,n_atoms+2,1,2,3,4)
        if i == 0:
            #check consistency between dipoles trajectory and atoms_per_par 
            if data.check_types(names_tmp):
                names_f = names_tmp
            else:
                message = "WARNING: atom types in " + data.primitive_pos + \
                          " do not correspond to those in " + data.coor + \
                          " This is not a problem if the atoms are sorted " + \
                          "in the same way. Continue at your own risk."
                print message
                # Replacing the types read from tdf with the types of data.coor
                names_f = data.names
        
        #loop over all atoms positions in given frame
        #and accumulate them together according to their charge type
        for k, at in enumerate(names_f):
            n_idx = atoms_per_par.keys().index(at)
            #n_idx = name_idxs[at]
            offset=i*3
            current_tdf[n_idx][offset+0] += x_t[k]
            current_tdf[n_idx][offset+1] += y_t[k]
            current_tdf[n_idx][offset+2] += z_t[k]

    #average atoms positions, three independent values (x,y,z) for each parameter
    if len(current_tdf[0]) > 3:
        avg_traj = np.zeros([len(current_tdf),3])
        for i,col in enumerate(current_tdf):
            avg_traj[i][0] = np.average(col[0::3])
            avg_traj[i][1] = np.average(col[1::3])
            avg_traj[i][2] = np.average(col[2::3])
            current_tdf[i][0::3] = np.subtract(current_tdf[i][0::3], avg_traj[i][0])
            current_tdf[i][1::3] = np.subtract(current_tdf[i][1::3], avg_traj[i][1])
            current_tdf[i][2::3] = np.subtract(current_tdf[i][2::3], avg_traj[i][2])

    # The entire trajectory is on one single frame xyz with n_atom_tdf entries
    n_atom_tdf = nstep  

    # B_j elements for TDF (reference data)
    x_f,y_f,z_f = atautils.tdf_smooth(data.primitive_tdf,
                  data.tdf_from_cp2k, nstep, print_tdf=False)

    # The average is computed only if there is more than one dipole 
    if len(current_tdf[0]) > 3:
        x_avg, y_avg, z_avg = np.average(x_f), np.average(y_f), np.average(z_f)
        x_f = np.subtract(x_f, x_avg)
        y_f = np.subtract(y_f, y_avg)
        z_f = np.subtract(z_f, z_avg)

    traj_f = np.array(zip(x_f,y_f,z_f)).flatten()
    
    in_f_pos.close()
    
    # current tdf is already the transpose (at) of the A matrix
    at = current_tdf
    b = traj_f
    ata = np.dot(at, np.transpose(at))
    atb = np.dot(at,b)
    btb = np.dot(np.transpose(b),b)

    return atb, ata, btb
 

def iterative_round_charges(charges, sorted_n_atoms_per_charge, rnd=4):
    """
    Round all numbers in the working array to the wanted decimal place,
    trying to keep the total sum equal to zero.
    
    The function converges when the total sum of the product of
    the working array times the corresponding multiplicity of each element becomes
    smaller than 1.0e-10. After 1000 iterations the function moves the
    rounding to the next decimal place on the right. If there is no convergence
    before reaching the 15th decimal place the function gives up.

    Parameters
    ----------
    charges : (N,) float ndarray
        working array of numbers to be rounded
    sorted_n_atoms_per_charge : (N,) int ndarray
        multiplicity of each element in the working array
        (e.g., number of atoms sharing same charge value), the order follows
        the one of "charges"
    rnd : int
        decimal place position to be rounded

    Returns
    -------
    new_uc_charge : float
        sum value of the product between the rounded working array times
        the corresponding multiplicity of each element
    new_ch : (N,) float ndarray
        rounded working array
    rnd : int
        decimal place position that has been actually rounded
    """

    rounded_charges = np.around(charges,rnd)
    rounded_charges_2 = np.zeros(len(charges))
    #take the values of working array and round them to the rnd decimal
    #place, generate as well a vector with the opposite rounding
    #(e.g., 0.04337 becomes 0.0434 in one array and 0.0433 in the other)
    for i,ch in enumerate(charges):
        offset = np.power(10.0,(-rnd))
        if ch - rounded_charges[i] < 0: # ch = 0.12337 rounded = 0.1234 
            offset = -offset
        rounded_charges_2[i] = rounded_charges[i] + offset
    
    fused = zip(rounded_charges, rounded_charges_2)
        
    select = np.zeros(2*len(charges),dtype=int)
    select[:len(charges)] = 1
    assert np.sum(select) == len(charges)
    
    counter = 0
    while True:
        #randomly choose the up or down rounding for each value
        np.random.shuffle(select)

        new_ch = np.zeros(len(charges))
        for i in xrange(len(charges)):
            new_ch[i] = fused[i][select[i]]
    
        new_uc_charge = np.sum(np.multiply(sorted_n_atoms_per_charge,new_ch))
        counter += 1
        #break if no neutrality found
        if rnd > 14:
            print "Iterative rounding uc charge not converged."
            break
        
        #neutrality found
        if abs(new_uc_charge) < 1.0e-10:
            break
        
        #move rounding one decimal place to the right
        if counter >= 1000:
            rnd +=1
            new_uc_charge, new_ch, rnd = iterative_round_charges(charges, sorted_n_atoms_per_charge, rnd)
            break

    return new_uc_charge, new_ch, rnd


def add_lambda_elements(atb, ata, data):
    """
    Return the ATA and ATB matrices to solve the linear system
    ATAq = ATB after adding the lagrange multiplier to constrain
    the total charge


    Parameters
    ----------
    atb : (N,) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and the reference data column matrix (of size M)
    ata : (N,N) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and itself (of size M,N)
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    B : (N+1,) float ndarray 
        equal to atb with one extra element containing the desired total charge value
    A : (N+1,N+1) float ndarray 
        equal to ata with one extra row and column having values corresponding to the
        number of atoms belonging to the same type (e.g., without kind symmetry the
        values are set to 1). The special element A[N+1][N+1] is set to 0
        (it refers to the Lagrange multiplier).
    """
    
    #charge neutrality equation line
    neu_line = np.array(data.natoms_per_par, ndmin=2, dtype=float)
    neu_prm = np.concatenate((ata, neu_line), axis=0)
    
    # Create extra column and add an extra element (row) to this column
    column_extra = np.concatenate((np.transpose(neu_line), np.zeros((1,1), dtype=float)), axis=0)

    #add  the extra column (axis 1) for lambda parameter
    A = np.concatenate((neu_prm, column_extra), axis=1)
    charge_wanted=np.array([data.q_tot])
    B = np.concatenate((atb, charge_wanted))
    
    return B,A


def compute_chi(btb, atb, ata, q, nvals):
    """
    Return the root of the mean squared (RMS) error between reference and
    model data for a given set of parameters q.

    Parameters
    ----------
    atb : (N,) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and the reference data column matrix (of size M)
    ata : (N,N) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and itself (of size M,N)
    btb : float 
        dot product between the transposed reference data column matrix and 
        itself, this is the magnitude of the reference data
    q : (N,) float ndarray
        array of parameters
    nvlas : number of equations contributing to the RMS

    Returns
    -------
    RMS : float
        the root mean squared error
    """

    #X = || B - Aq ||**2
    #RMS = sqrt( BTB - 2 qT ATB + qT ATA q )
    qt = np.transpose(q)
    return np.sqrt((btb - 2.0*np.dot(qt,atb) + np.dot(np.dot(qt,ata),q))/nvals)

def precompute_working_matrices_vdW(data, cube_fn, miss_g):
    """
    Return the ATA and ATB matrices together with BTB and the average model
    data per parameter for a set of different vdW scaling radius gamma.
    Data useful to restart the calculations (number of valid grid points and
    the sum of the squared ESP values of the full grid) are also returned.

    Parameters
    ----------
    data : instance of class Settings
        contains input and check variables
    cube_fn : string
        the name of the cube file containing the ESP data
    miss_g : float list
        the values of gamma for which the matrices are to be computed

    Returns
    -------
    results : list of dict (one for each value of miss_g)
        each dict contains:
        atb : (N,) float ndarray 
            dot product between the transposed model data matrix (of size N,M)
            and the reference data column matrix (of size M)
        ata : (N,N) float ndarray
            dot product between the transposed model data matrix (of size N,M)
            and itself (of size M,N)
        btb : float 
            dot product between the transposed reference data column matrix and 
            itself, this is the magnitude of the reference data
        avg_prm_traj_f : (N,) float ndarray
            contains the average value per column of matrix A (M,N)
        n_pts: int
            number of valid grid points
        e2_tot: float
            sum of the squared ESP values of the full grid
    """
    
    atautils.init(cube_fn)
    results = []

    # compute the ESP model data only for the values of gamma contained in miss_g
    for gamma in miss_g:
        # remove points lying inside a sphere of radius  vdW * gamma
        ngood_grid_pts, e2_tot = atautils.filter(gamma, False)
        # compute model data
        e_atb, e_ata, e_btb, avg_prm_traj_f = atautils.get_ata_matrix(15, data.locut, -1.0, ngood_grid_pts, data.natoms_per_par, data.kind_idx, False, data.esp_sign)
        # store 
        matrices = dict(zip(['tail','atb', 'ata', 'btb', 'avg_prm_traj_f','n_pts', 'e2_tot'],[gamma,e_atb, e_ata, e_btb, avg_prm_traj_f,ngood_grid_pts,e2_tot]))
        results.append(matrices)
    
    # get atomic numbers for this cube
    tmp_z = atautils.get_zatom(data.natoms)
    assert save_z_and_check(data,tmp_z), "ERROR: atomic numbers are not consistent across cubes"
    print cube_fn
    # cleanup of the stored data in the fortran module
    atautils.release()

    return results

def precompute_working_matrices_frame(data):
    """
    Return the ATA and ATB matrices together with BTB accumulated over a set of
    different system configurations for a set of different vdW scaling radius
    gamma. 
    Data useful to restart the calculations (sum of the number of valid
    grid points and sum of the sum of the squared ESP values of the full grid)
    are also returned.

    Parameters
    ----------
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    results : list of dict (one for each value of gamma)
        each dict contains:
        atb : (N,) float ndarray 
            sum over system configurations of 
            dot product between the transposed model data matrix (of size N,M)
            and the reference data column matrix (of size M) 
        ata : (N,N) float ndarray
            sum over system configurations of 
            dot product between the transposed model data matrix (of size N,M)
            and itself (of size M,N)
        btb : float 
            sum over system configurations of 
            dot product between the transposed reference data column matrix and 
            itself, this is the magnitude of the reference data
        avg_prm_traj_f : None
            this quantity is meaningless for multiple configurations computation
        n_pts: int
            sum over system configurations of the
            number of valid grid points
        e2_tot: float
            sum over system configurations of the
            sum of the squared ESP values of the full grid
    """

    # computing the set of gamma over which the computations will be performed
    # and store them inside data
    g_stop = data.g_start + data.g_step * (data.n_gammas - 1)
    data.gammas = np.linspace(data.g_start,g_stop,data.n_gammas)

    # check which gammas are already computed and read them,
    # also store which ones are missing (miss_g) 
    miss_g = []
    matrices = []
    for gamma in data.gammas:
        directory = ata_prefix + str(gamma) + "/"
        if path.exists(directory):
            mat = read_ata(directory)
            matrices.append(mat)
            if "atomic_data" in mat:
                # get and check atomic number
                tmp_z = mat["atomic_data"]
                assert save_z_and_check(data,tmp_z), "ERROR: atomic numbers are not consistent across ata directories"
            print "Reading ESP data from: ", directory
        else:
            miss_g.append(gamma)
            matrices.append("missing")

    # compute the (possibly) missing values
    if miss_g:
        print "Computing ESP data for gammas: "
        for ga in miss_g:
            print ga
        assert data.cube_files, "ERROR: cube files not found while generating missing data for gamma " + str(ga) 
        cube_files = data.cube_files
        nframes = len(cube_files)
        tmp_results = []
        # actually compute missing atas
        for cube_fn in cube_files:
	    tmp_results.append(precompute_working_matrices_vdW(data,cube_fn,miss_g))
        print "Number of cubes is: ", nframes 

        # store results for first frame 
        miss_res = tmp_results[0]
        # loop over remaining frames and accumulate over the first
        for tres in tmp_results[1:]:
            # loop over gammas
            for i,gres in enumerate(miss_res):
                gres['atb'] += tres[i]['atb']
                gres['ata'] += tres[i]['ata']
                gres['btb'] += tres[i]['btb']
                gres['n_pts'] += tres[i]['n_pts']
                # avg_prm_traj_f now are meaningless so they are nullified
                gres['avg_prm_traj_f'] = None
                gres['e2_tot'] += tres[i]['e2_tot']

        # set e2_tot reference for check
        try:
            cmp_e2 = matrices[0]['e2_tot']
        except:
            cmp_e2 = miss_res[0]['e2_tot']
        # fuse the missing matrices with the read ones
        cnt = 0
        for i,mat in enumerate(matrices):
            if mat == "missing":
                matrices[i] = miss_res[cnt]
                cnt += 1
            # check that all gammas have data taken from same cubes
            assert abs(cmp_e2 - matrices[i]['e2_tot']) < 10e-4,\
            "ERROR: it seems that cubes used to generate ESP data for gamma " + str(matrices[0]['tail']) + " and " +\
            str(matrices[i]['tail']) + " do not coincide. "+\
            "Total sum of squared ESP data is: %s" % str(cmp_e2) + " vs " + str(matrices[i]['e2_tot'])
        # if everything is ok save the missing matrices
        for mat in miss_res:
            save_ata(mat, data.zatoms)
    else:
        # check that all gammas have data taken from same cubes (only for the read ata case)
        cmp_e2 = matrices[0]['e2_tot']
        for i,mat in enumerate(matrices):
            assert abs(cmp_e2 - matrices[i]['e2_tot']) < 10e-4,\
            "ERROR: it seems that cubes used to generate ESP data for gamma " + str(matrices[0]['tail']) + " and " +\
            str(matrices[i]['tail']) + " do not coincide. "+\
            "Total sum of squared ESP data is: %s" % str(cmp_e2) + " vs " + str(matrices[i]['e2_tot'])

    return matrices

def solve(atbs, atas, btbs, neqs, weights, data):
    """
    Return the set of parameters minimizing the sum of the squares
    of the errors (least squares) between reference and model data.

    The overall set of equations to be solved is made up by C (e.g., one or
    more) independent contributions, whose influence over the total
    merit function is controlled through opportune weights.

    Parameters
    ----------
    atbs : (C,(N,)) float ndarray of (N,) float ndarray
        where (N,) comes from the dot product between the transposed model
        data matrix (of size N,M)
        and the reference data column matrix (of size M)
    atas : (C,(N,N)) float ndarray of (N,N) float ndarray 
        where (N,N) comes from the dot product between the transposed model
        data matrix (of size N,M) and itself (of size M,N)
    btbs : (C,) float ndarray
        where each element comes from the dot product between the transposed
        reference data column matrix and itself,
        this is the magnitude of the reference data
    neqs : (C,) int ndarray
        array containing the number of equations for each independent
        contribution
    weights : (C,) float ndarray
        array containing the weight of each independent contribution over the
        merit function
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    res : (N,) float ndarray
        array of optimized parameters
    res_full : (N+1,) float ndarray
        array of optimized parameters plus the value of extra Lagrange multiplier
        for unit cell charge neutralization
    hogs : (N,) float ndarray
        array containing the relative root mean squared error (RRMS) for the
        resulting parameters, evaluated independently for each contribution
    hchis : (N,) float ndarray
        array containing the root mean squared error (RMS) for each contribution
    rchis: (N,) float ndarray
        array containing the square of the mean magnitude of the reference data, for each
        contribution
    """

    # fuse together the various contributions, taking into account their
    # weight and magnitude (the latter is needed to normalize each
    # contribution
    mgnts = weights**2/btbs
    atb = np.zeros_like(atbs[0])
    ata = np.zeros_like(atas[0])
    for i in range(len(atbs)):
        atb += atbs[i]*mgnts[i]
        ata += atas[i]*mgnts[i]

    # add the u.c. total charge neutralization
    B,A =  add_lambda_elements(atb, ata, data)
    res_full = np.linalg.solve(A,B)
    res = res_full[:-1]
    
    # evaluate match independently for each contribution
    hchis = [compute_chi(btbs[i], atbs[i], atas[i], res, neqs[i]) for i in range(len(atas))]
    rchis = np.sqrt(btbs/neqs)
    hogs = hchis/rchis
    return res, res_full, hogs, hchis, rchis
 
def apply_restraints(data, e_atb, e_ata, e_btb):
    """
    Return the set of parameters minimizing the sum of the squares
    of the errors (least squares) between ESP reference and model data,
    applying simultaneously parabolic restraints to the parameters.

    One or more parameters can be restrained, balancing
    their contribution to the figure of merit through global or
    individual weights.
    X = ... + w_1**2(Q_1-q_1)**2 + w_2**2(Q_2-q_2)**2 + ...
    where w_i is the weight of the ith restraint, Q_i the restrained
    value itself, and q_i the parameter to be restrained.

    Parameters
    ----------
    atb : (N,) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and the reference data column matrix (of size M)
    ata : (N,N) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and itself (of size M,N)
    btb : float 
        dot product between the transposed reference data column matrix and 
        itself, this is the magnitude of the reference data
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    res : (N,) float ndarray
        array of optimized parameters

    Notes
    -----
    The ESP part contribution to the merit function is normalized by the
    ESP magnitude to allow the use of restraint weights which are independent
    from the system size.
    Moreover, the restraint contribution is not normalized by the restraint
    magnitude to allow the possible restraining of all parameters to 0.0.
    """

    # vector containing all weights
    r_weight = np.zeros(data.npar, dtype = "float")
    # vector containing all restraints
    r_b = np.zeros(data.npar, dtype = "float")
    for i,key in enumerate(data.atoms_per_par.keys()):
        if key in data.restraint_weight:
            r_weight[i] = data.restraint_weight[key]
            r_b[i] = data.restraints[key]

    # restraints must be multiplied by their atom-type symmetry population
    r_a = np.identity(data.npar,dtype=float)
    r_atb = r_a.T.dot(r_b)*data.natoms_per_par
    r_ata = r_a.T.dot(r_a)*data.natoms_per_par
    # the magnitude is not used to allow restraining to 0.0
    #r_btb = np.sum(r_b*r_b*data.natoms_per_par)

    # fuse restraints and ESP
    atb = r_atb*r_weight**2 + e_atb/e_btb
    ata = r_ata*r_weight**2 + e_ata/e_btb

    # solve
    B,A =  add_lambda_elements(atb, ata, data)
    res_full = np.linalg.solve(A,B)
    res = res_full[:-1]

    return res


def save_z_and_check(data, tmp_z):
    """
    Check the consistency of atomic numbers across .cube files and atoms
    already processed.
    Return True if the atomic number contained in tmp_z coincide with
    the ones saved inside data. If no atomic numbers are contained in 
    data then tmp_z is copied inside it.
    
    Parameters
    ----------
    data : instance of class Settings
        contains input and check variables
    tmp_z : (N,) int ndarray
        contains the atomic numbers of the system atoms

    Returns
    -------
    bool :
        True if tmp_z coincide with the ndarray zatoms contained in data or
        if zatoms is not present in data and tmp_z is saved.
    """
    
    if not hasattr(data, 'zatoms'):
        data.zatoms = tmp_z
        assert type_to_z(data), "ERROR: atomic numbers are not consistent across each atom type group"
        return True
    
    return np.array_equal(data.zatoms, tmp_z)

def type_to_z(data):
    """
    This function associates an atomic nuber z to each atom type
    and returns true if all elements of the same type have the same z.
    
    Parameters
    ----------
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    bool  :
        True if all elements of the same type have the same z.

    Notes
    -----
    Cube files should have atoms in the same order of .xyz (or .pdb) file 
    provided in the input (i.e., the atom types should coincide).
    This cannot be automatically checked but it is possible to compare .xyz atom type
    with corresponding atomic number in .cube file (i.e., the atom kind can
    be checked).
    """

    print "Atom types and atomic number:"
    # associate one atomic number Z to each key
    atom_grp_index = 0 # index of atoms when these are grouped into types
    z_per_par=np.zeros(data.npar,dtype=int)
    for i,key in enumerate(data.atoms_per_par.keys()):
        # we want the index of the first atom of each type
        # found using the "atom_grp_index"
        # fed into the "keep_idx" array
        # kind_idx[atom_grp_index] index of atom in pdb file
        z_per_par[i]=(data.zatoms[data.kind_idx[atom_grp_index]-1])
        print "{:<4s}-->{:>4d}".format(key,z_per_par[i])

        # check if all atoms belonging to same atom-type symmetry group
        # have same Z
        for j in range(atom_grp_index, atom_grp_index+data.atoms_per_par[key]):
            if (z_per_par[i] != data.zatoms[data.kind_idx[j]-1]):
                return False

        # update index
        atom_grp_index += data.atoms_per_par[key]

    data.z_per_par = z_per_par

    return True

def apply_rg(data, e_atb, e_ata, e_btb):
    """
    Return the set of parameters minimizing the sum of the squares
    of the errors (least squares) between ESP reference and model data,
    applying simultaneously to the parameters "Rappe-Goddard like" restraints
    exactly as implemented in REPEAT.

    All parameters are restrained according to the electronegativity and the
    self-Coulomb interaction corresponding to the specific atom kind.
    X = ... + w_1(E_1 + Y_1*q_1 + 0.5*J_1*q_1**2) + ...
    where w_i is the weight of the ith restraint, Y_i the
    electronegativity, J_i the self_Coulomb term,
    see "Campana et al., JCTC 5 (2009) 2866-2878".
    
    Parameters
    ----------
    atb : (N,) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and the reference data column matrix (of size M)
    ata : (N,N) float ndarray 
        dot product between the transposed model data matrix (of size N,M)
        and itself (of size M,N)
    btb : float 
        dot product between the transposed reference data column matrix and 
        itself, this is the magnitude of the reference data
    data : instance of class Settings
        contains input and check variables

    Returns
    -------
    res : (N,) float ndarray
        array of optimized parameters

    Notes
    -----
    This kind of restraint is equivalent to a parabolic restraint of
    the parameters to a charge Q_i=(-Y_i/J_i), multiplying the weight of
    each restraint by J_i.
    """

    # find the atomic number of the atom kind related to each parameter
    # check if zatoms are stored and if not retrieve them from the
    # first cube file
    find_zatom(data)
    # self-Coulomb per parameter
    joo=data.selci[data.z_per_par-1]*atautils.autokcal
    # electronegativity per parameter
    xoo=data.eneg[data.z_per_par-1]*atautils.autokcal
    # rg_weight (input in au units as REPEAT) converted to kcal/mol
    rg_weight=data.rg_weight*atautils.autokcal

    print ""
    print "This is equivalent to restraining the charges, "
    print "with a weight depending on the self Coulomb interaction J, "
    print "to values: "
    print -data.eneg[data.z_per_par-1]/data.selci[data.z_per_par-1]
    print ""

    
    r_atb =(-data.natoms_per_par*rg_weight*0.5*xoo)
    r_ata =(
    np.identity(data.npar,dtype=float)*
    data.natoms_per_par*
    rg_weight*0.5*
    joo)

    # fuse restraints and ESP
    atb = r_atb + e_atb
    ata = r_ata + e_ata

    # solve
    B,A =  add_lambda_elements(atb, ata, data)
    res_full = np.linalg.solve(A,B)
    res = res_full[:-1]

    return res

def save_ata(mat, atomic_data=None):
    """
    Write on file all data required to restart the calculations.

    Parameters
    ----------
    mat : dict
        containing:
        tail : string
            complete ata_prefix to identify where data are saved
        n_pts : int
            sum over system configurations of the
            number of valid grid points
        atb : (N,) float ndarray 
            sum over system configurations of 
            dot product between the transposed model data matrix (of size N,M)
            and the reference data column matrix (of size M) 
        ata : (N,N) float ndarray
            sum over system configurations of 
            dot product between the transposed model data matrix (of size N,M)
            and itself (of size M,N)
        btb : float
            sum over system configurations of the
            dot product between the transposed reference data column matrix and
            itself, this is the magnitude of the reference data
        e2_tot: float
            sum over system configurations of the
            sum of the squared ESP values of the full grid
    
    atomic_data : list
        can contain numbers or string. For ESP data it is the list of atomic
        numbers, for TDF data it is the list of atomic types.
    """

    directory = ata_prefix + str(mat['tail']) + "/"
    if not path.exists(directory):
        makedirs(directory)
    else:
        print "WARNING: ", directory, " overwritten."
    ataf = open(directory + "data.dat", "w")
    ataf.write("tail " + str(mat['tail']) + "\n")
    ataf.write("n_pts " + str(mat['n_pts'])+ "\n")
    ataf.write("btb {:>25.15f}\n".format(mat['btb']))
    ataf.write("e2_tot {:>18.15e}\n".format(mat['e2_tot']))
    mat['atb'].tofile(directory + "atb.bin") 
    mat['ata'].tofile(directory + "ata.bin") 
    if atomic_data != None:
        ataf.write("atomic_data: ")
        for ad in atomic_data:
            ataf.write(str(ad) + " ")
        ataf.write("\n")

def read_ata(atadir):
    """
    Read data saved with save_ata.

    Parameters
    ----------
    atadir : string
        name of the directory where files are stored

    Returns
    -------
    results : dict 
        containing:
        tail : string
            complete ata_prefix to identify where data are saved
        n_pts : int
            sum over system configurations of the
            number of valid grid points
        atb : (N,) float ndarray 
            sum over system configurations of 
            dot product between the transposed model data matrix (of size N,M)
            and the reference data column matrix (of size M) 
        ata : (N,N) float ndarray
            sum over system configurations of 
            dot product between the transposed model data matrix (of size N,M)
            and itself (of size M,N)
        btb : float 
            sum over system configurations of 
            dot product between the transposed reference data column matrix and 
            itself, this is the magnitude of the reference data
        avg_prm_traj_f : None
            this quantity is meaningless for multiple configurations computation
        e2_tot: float
            sum over system configurations of the
            sum of the squared ESP values of the full grid
        atomic_data (optional): list 
            can contain numbers or string. For ESP data it is the list of atomic
            numbers, for TDF data it is the list of atomic types.
    """

    fnames = [atadir + name for name in ["data.dat","atb.bin","ata.bin"]]
    for name in fnames:
        assert path.getsize(name) > 0, "ERROR corrupted file: %s" % name

    lines = open(fnames[0]).readlines()
    slines = [line.split() for line in lines if not "atomic_data:" in line]
    atomic_data = np.squeeze(np.array([ad.replace("atomic_data:","").split() for ad in lines if "atomic_data:" in ad]))
    if all(ad.isdigit() for ad in atomic_data):
        atomic_data = atomic_data.astype(np.int)
    results = {key: value for (key,value) in slines}

    if results['tail'].isdigit():
        results['tail'] = float(results['tail']) 
    results['n_pts'] = int(results['n_pts'])
    results['btb'] = float(results['btb'])
    results['e2_tot'] = float(results['e2_tot']) 
        
    
    results['avg_prm_traj_f'] = None 
    results['atb'] = np.fromfile(fnames[1]) 
    results['ata'] = np.fromfile(fnames[2]).reshape(len(results['atb']),-1)
    
    if atomic_data.size != 0:
        results['atomic_data'] = atomic_data
    
    return results


def find_zatom(data):
    """
    Check if atomic numbers have been already stored (either read from 
    a cube or from saved data). If not the function opens the firs cube
    file and stores them.

    Parameters
    ----------
    data : instance of class Settings
        contains input and check variables
    """
    
    if not hasattr(data, 'zatoms'):
        cube_fn = data.cube_files[0]
        atautils.init(cube_fn)
        data.zatoms = atautils.get_zatom(data.natoms)
        assert type_to_z(data), "ERROR: atomic numbers are not consistent across each atom type group"
        atautils.release()
    

def main():
    
    dt = settings()
    print ""

    # TDF data generation
    if "ESP+TDF" == dt.fit_kind:
        directory = dt.ata_tdf 
        if path.exists(directory):
            print "Reading TDF data from: ", directory
            mat=read_ata(directory)
            b_atb = mat['atb']
            b_ata = mat['ata']
            b_btb = mat['btb']
            tdf_nequations = mat['n_pts']
            # check that the atom-types used for the generation of TDF
            # atas are consistent with current .xyz or (.pdb) ones
            assert Counter(dt.names) == Counter(mat['atomic_data'])
        else:
            print "Computing TDF data"
            # read TDF data and subtract avgs
            b_atb, b_ata, b_btb, tdf_natoms, tdf_nstep =  read_full_traj_tdf_auto(dt)
            tdf_nequations = tdf_nstep*3
            mat =  dict(zip(['tail','atb', 'ata', 'btb', 'avg_prm_traj_f','n_pts', 'e2_tot'],['tdf',b_atb, b_ata, b_btb, None, tdf_nequations, 0.0]))
            save_ata(mat,dt.names)
            assert tdf_natoms == dt.natoms

        # solve and store reference TDF charges, these will not change during the entire procedure
        res, res_ber, hog_ber, hchi_b, rchi_b = solve(np.array([b_atb]), np.array([b_ata]), np.array([b_btb]),
                np.array([tdf_nequations]),np.array([1.0]), dt)
        
        print "Pure TDF charges:"
        for k,knd in enumerate(dt.atoms_per_par.keys()):
            print "{:<4s}={:>14.8f}".format(knd, res[k])

        uc_charge =  np.sum(np.multiply(dt.natoms_per_par,res))
        if abs(uc_charge) > Tiny:
            print "WARNING! Unit cell charge: ", uc_charge, "\n"
        print "TDF RRMS: ", hog_ber[0]
        print ""

        # file containing the most relevant results for all gammas
        res_h =  open("auto_weight_res.dat", "w")
        res_h.write("# gamma   rounded charges # iteration_stop num_grid_pts \n")



    # generate atas for all ESP configurations and for all gammas
    g_matrices = precompute_working_matrices_frame(dt)

    # main loop over different gammas (i.e., cube points ensembles for ESP)
    for matrices in g_matrices: 
        gamma = matrices['tail']
        print "========================================================================"
        print "Working for gamma: ", gamma
        
        # solve and store reference ESP charges, these will not change for a given gamma
        ngood_grid_pts = matrices['n_pts']
        e_atb = matrices['atb'] 
        e_ata = matrices['ata']
        e_btb = matrices['btb']
        avg_prm_traj_f = matrices['avg_prm_traj_f']
        res, res_esp, hog_esp, hchi_e, rchi_e = solve(np.array([e_atb]), 
            np.array([e_ata]), np.array([e_btb]), np.array([ngood_grid_pts]), 
                np.array([1.0]), dt)
        print "Pure ESP charges:"
        for k,knd in enumerate(dt.atoms_per_par.keys()):
            print "{:<4s}={:>14.8f}".format(knd, res[k])
        
        uc_charge =  np.sum(np.multiply(dt.natoms_per_par,res))
        if abs(uc_charge) > Tiny:
            print "WARNING! Unit cell charge: ", uc_charge, "\n"
        print "ESP RRMS: ", hog_esp[0]
        print ""


        #solve for restrained fit
        if ("RESP" == dt.fit_kind):
            res = apply_restraints(dt, e_atb, e_ata, e_btb)
        
        #solve for "Rappe-Goddard like" restrained fit
        if ("RG" == dt.fit_kind):
            res = apply_rg(dt, e_atb, e_ata, e_btb)
            uc_charge =  np.sum(np.multiply(dt.natoms_per_par,res))
            if abs(uc_charge) > Tiny:
                print "WARNING! Unrounded Unit cell charge: ", uc_charge, "\n"

        # if "non TDF" fit, computation for this gamma ends here
        if not ("ESP+TDF" == dt.fit_kind):
            # rounding resulting charges starting from fourth decimal place
            uc_rnd_ch, rnd_ch, rnd = iterative_round_charges(res, dt.natoms_per_par, rnd=4)
            
            # evaluate RRMS for rounded charges
            rnd_hchie = compute_chi(e_btb, e_atb, e_ata, rnd_ch, ngood_grid_pts)
            
            print "Final rounded charges:"
            fmt="{:<4s}={:>10."+str(rnd)+"f}"
            for k,knd in enumerate(dt.atoms_per_par.keys()):
                print fmt.format(knd, rnd_ch[k])
            
            if abs(uc_rnd_ch) > Tiny:
                print "WARNING! Unit cell charge: ", uc_rnd_ch, "\n"
            print "RRMS: ", rnd_hchie/rchi_e[0]
            
            continue



        # ESP+TDF fit only
        weight_b=dt.weight_tdf

        res_f =  open("results_" + str(gamma) + ".dat", "w")
        res_u =  open("results_unrounded_" + str(gamma) + ".dat", "w")
        res_lf = open("res_loop_" + str(gamma) + ".dat", "w")
        res_lf.write("# weight_tdf   delta%_tdf   delta%_esp    RRMS_tdf    RRMS_esp \n")
            

        # arrays to treat together TDF and ESP data
        atbs = np.array([b_atb,e_atb])
        atas = np.array([b_ata,e_ata])
        btbs = np.array([b_btb,e_btb])
        neqs = np.array([tdf_nequations,ngood_grid_pts])
        done_auto = False
        auto_res = None

        # loop over various TDF weights (increasing values)
        print "Performing weight loops..."
        for w in xrange(dt.nw_loops):
            weight_b += dt.w_incr
            # weight should fall in [0,1]
            assert  weight_b <= 1.0+Tiny and weight_b >= 0.0-Tiny, "weight_b = %r" % weight_b
            weight_e = 1.0-weight_b
            weights = np.array([weight_b,weight_e])
            res, res_full, hogs, hchis, rchis = solve(atbs, atas, btbs, neqs,
                                 weights, dt)

            # evaluate u.c. total charge after fit
            uc_charge =  np.sum(np.multiply(dt.natoms_per_par,res))

            # evaluate RMS independently for ESP and TDF
            hchie = hchis[1]
            hchib = hchis[0]
            # compute the square of the mean magnitude of the reference data
            # these have been already computed for pure TDF and ESP
            rchi_e = rchis[1]
            rchi_b = rchis[0]
            # store current RRMS for ESP
            current_hge = hogs[1]
            # store current RRMS for TDF
            current_hgb = hogs[0]
            # evaluate improvement/degradation with respect to best ESP and best TDF
            delta_ber = (current_hgb - hog_ber[0])/hog_ber[0]
            delta_esp = (current_hge - hog_esp[0])/hog_esp[0]
            res_lf.write("{:>7.4f} {:>20.16f} {:>20.16f} {:>20.16f} {:>20.16f}\n".format(weight_b,
            delta_ber*100.0, delta_esp*100.0, current_hgb, current_hge))
            
            # after first iteration, check if auto-stop criterion is met
            # i.e., [(RRMS_E - RRMS_Ebest)/RRMS_Ebest] < [(RRMS_T - RRMS_Tbest)/RRMS_Tbest]
            if w > 0 and not done_auto:
                if (delta_esp - delta_ber) > 0.0:
                    auto_res = [[rnd_ch,rnd], w-1, (delta_esp - delta_ber), rnd_hchib/rchi_b, rnd_hchie/rchi_e, uc_rnd_ch]
                    done_auto = True

            # write resulting unrounded parameters for this weight
            res_u.write(str(gamma))
            for val in res:
                res_u.write("{:>20.15f}, ".format(val))
            res_u.write(" # " + str(ngood_grid_pts) + "\n")

            # rounding resulting charges starting from fourth decimal place
            uc_rnd_ch, rnd_ch, rnd = iterative_round_charges(res, dt.natoms_per_par, rnd=4)
            
            # evaluate RMS independently for ESP and TDF for rounded charges
            rnd_hchie = compute_chi(e_btb, e_atb, e_ata, rnd_ch, ngood_grid_pts)
            rnd_hchib = compute_chi(b_btb, b_atb, b_ata, rnd_ch, tdf_nequations)

            # write results for this weight
            res_f.write(str(gamma))
            fmt = "{:>"+ str(rnd+4) + "." + str(rnd) + "f}, "
            for val in rnd_ch:
                res_f.write(fmt.format(val))
            res_f.write(" # " + str(ngood_grid_pts) + "\n")

            # after first iteration, check if minimum distance criterion is met
            # i.e., min(abs(RRMS_E-RRMS_T))
            delta = np.abs(current_hge-current_hgb)
            if w > 0:
                if (delta < delta_pre):
                    delta_res = [[rnd_ch,rnd], w, delta, rnd_hchib/rchi_b, rnd_hchie/rchi_e, uc_rnd_ch]
            else:
                delta_res = [[rnd_ch,rnd], w, delta, rnd_hchib/rchi_b, rnd_hchie/rchi_e, uc_rnd_ch]
            delta_pre = delta
            
        res_lf.close()
        res_u.close()
        res_f.close()

        # final results:
        # chose between "autoweight" or "min(abs(RRMS_E-RRMS_T))" criterion
        # NOTE: true_res[0] is a list:
        # where true_res[0][0] is a list containing the rounded charges
        # and   true_res[0][1] is the final number of decimal places used for the rounding
        true_res = delta_res
        stop_kind = "min(abs(RRMS_E-RRMS_T))"
        if auto_res and delta_res[1] > auto_res[1]:
            true_res = auto_res
            stop_kind = "RRMS_E degradation larger than RRMS_T improvement"

        # case of single imposed weight
        true_tdf_w = true_res[1]*dt.w_incr
        if dt.nw_loops == 1:
            true_tdf_w = weight_b
            stop_kind = "chosen by the user"
        
        print "Stopping due to: ", stop_kind
        print "at weights: ESP ", 1.0-true_tdf_w, "   TDF ", true_tdf_w
        print "Final rounded charges: "
        fmt = "{:<4s}={:>" + str(6+true_res[0][1]) + "." + str(true_res[0][1]) + "f}"
        for k,knd in enumerate(dt.atoms_per_par.keys()):
            print fmt.format(knd, true_res[0][0][k])
        print "ESP RRMS: ", true_res[4]
        print "TDF RRMS: ", true_res[3]
        if abs(true_res[5]) > Tiny:
            print "WARNING! Unit cell charge: ", true_res[5], "\n"
        
        # write final rounded charges
        res_h.write(str(gamma))
        fmt = "{:>" + str(4+true_res[0][1]) + "." + str(true_res[0][1]) + "f}, "
        for val in true_res[0][0]:
            res_h.write(fmt.format(val))
        res_h.write(" # " + str(true_res[1]) + " " + str(ngood_grid_pts) + "\n")

    # needed because outside the weights loop
    if ("ESP+TDF" == dt.fit_kind):
        res_h.close()
        

if __name__ == "__main__":
    print "Started on: ", datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    print "Current working directory is: ", getcwd()
    print "Pid is: ", getpid()
    print ""
    start_time = time()
    main()
    print ""
    print "Ended on: ", datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Run for seconds: ",time()-start_time
    print ""

    print "========================================================================"
    print "If you found this program useful, please cite: "
    print "M. Sant, A. Gabrieli, P. Demontis, G. B. Suffritti"
    print "Comput. Phys. Commun., accepted manuscript (2015)"
    print "http://dx.doi.org/10.1016/j.cpc.2015.10.005"
    print ""
    print "A. Gabrieli, M. Sant, P. Demontis, G. B. Suffritti"
    print "J. Chem. Theory Comput. 11 (2015) pp 3829-3843"
    print "http://dx.doi.org/10.1021/acs.jctc.5b00503"
    print "                                               "
    print "Together with the original REPEAT work:        "
    print "C. Campana, B. Mussard, T. K. Woo"
    print "J. Chem. Theory Comput. 5 (2009) 2866-2878"
    print "http://dx.doi.org/10.1021/ct9003405"
    print ""
    print "The logo has been generated with figlet: http://www.figlet.org/"
    print ""
    print "========================================================================"
    print "InfiniCharges comes with ABSOLUTELY NO WARRANTY. "
    print ""
