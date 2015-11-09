import atautils
import numpy as np
import time 
#Note:
#- to use modules inside a module use
# gfortran -c -fPIC cubecruncher.f90
# f2py -c -m atautils ata_utils.f90 cubecruncher.o
print atautils.atautils.autokcal
print atautils.atautils.init.__doc__
print atautils.atautils.filter.__doc__
print atautils.atautils.get_ata_matrix.__doc__
print atautils.atautils.release.__doc__
print atautils.atautils.tdf_smooth.__doc__
print atautils.atautils.get_zatom.__doc__
print atautils.atautils.read_xyz.__doc__
atautils.atautils.init("../tests/cubes/energy-force-v_hartree-1_0.cube")
t1 = time.clock()
n_good_grid_pts, e2_tot = atautils.atautils.filter(1.0, False)
print "ESP**2: ", e2_tot
#for i in np.linspace(0.5, 2.5, 21):
#    n_good_grid_pts = atautils.atautils.filter(i, False)
t2 = time.clock()
print "Time for filter: ", t2 - t1 
ell = np.array([16.9856, 16.9856, 16.9856])
#natoms_per_kind = [24,48,24,48,72,48,12]
#natoms_per_kind = [276]
natoms_per_kind = [24,48,24,48,48,24,48,12]
kind_idx=[i+1 for i in xrange(np.sum(natoms_per_kind))]

print atautils.atautils.get_zatom(sum(natoms_per_kind))
test_dip = False
if test_dip:
    ndip = len(open("../tests/tdf/tdf.dat").readlines()) - 4
    dpx = [i for i in range(ndip)]
    dpy = [i for i in range(ndip)]
    dpz = [i for i in range(ndip)]
    is_cp2k = False
    dpx, dpy, dpz = atautils.atautils.tdf_smooth("../tests/tdf/tdf.dat", is_cp2k, ndip, print_tdf=True)
    with open("check_sumdip", "w") as outf:
        for i,x in enumerate(dpx):
            outf.write("{:>10d}{:>20.10f}{:>20.10f}{:>20.10f}\n".format(i+1,x,dpy[i],dpz[i]))
    print "Check xyz"
    names, rx, ry, rz = atautils.atautils.read_xyz("../tests/tdf/pos.xyz", 276, ndip)
    #print names, rx, ry, rz
    names = ["".join(n).replace(" ","") for n in names]
    print names
    with open("check_pos.xyz", "w") as outf:
        for i in range(ndip):
            outf.write("276\n\n")
            for j in range(276):
                outf.write("{:>5s}{:>15.6f}{:>15.6f}{:>15.6f}\n".format(names[j],rx[j][i],ry[j][i],rz[j][i]))
        

print_pot=True 
res = atautils.atautils.get_ata_matrix(15, 1.0e-10, -1.0, n_good_grid_pts, natoms_per_kind, kind_idx, print_pot, -1.0)
print res
#for i in np.linspace(1.0, 2.0, 11):
#    n_good_grid_pts = atautils.atautils.filter(i, True)
#    res = atautils.atautils.get_ata_matrix(15, 1.0e-10, -1.0, n_good_grid_pts, natoms_per_kind, kind_idx, True, -1.0)
#    print res
atautils.atautils.release()

