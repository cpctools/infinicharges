from glob import glob
from collections import OrderedDict, Counter
from os import path, getcwd
import sys
from numpy import array
sys.path.insert(0, getcwd())
from InfiniInput import *

Tiny = 10e-9 #epsilon for comparisons
ata_prefix = "ata_"

class settings:
    'Container for all relevant variables.'
    data_dir = "./cubes/"
    tdf_data_dir = "./tdf/"
    #total charge must be 0.0 for periodic systems
    q_tot = 0.0
    #number of weight iterations
    nw_loops=101
    #implemented fit kinds
    valid_fit_kinds=["ESP+TDF","ESP","RESP","RG"]
    def __init__(self):
        self.fit_kind=fit_kind
        self.coor=coor
        self.g_start = g_start
        self.g_step = g_step
        self.n_gammas = n_gammas
        self.esp_sign = 1.0
        try:
            self.locut=locut
        except:
            self.locut=1.0e-10

        if (self.fit_kind in self.valid_fit_kinds):
            print "Performing ", self.fit_kind, " fit."
        else:
            raise Exception(self.fit_kind + " fit not recognized.")
        
        print "ESP *.cube files are expected to be in: " + self.data_dir
        self.cube_files = glob(self.data_dir + "*.cube")
        if (not self.cube_files):
            #raise Exception("*.cube file(s) not found.")
            print "WARNING: *.cube file(s) not found. The program will try to read saved data from: ./" + ata_prefix + "*"
            #pass
        if (esp_from_cp2k):
            self.esp_sign = -1.0
            print "ESP *.cube from CP2K, data are multiplied by: ", self.esp_sign
        
        if ("ESP+TDF" == self.fit_kind):
            self.tdf_from_cp2k=tdf_from_cp2k
            self.tdf_signal_fname = tdf_signal_fname
            self.tdf_coor_fname = tdf_coor_fname
            self.primitive_tdf = self.tdf_data_dir + self.tdf_signal_fname
            self.primitive_pos = self.tdf_data_dir + self.tdf_coor_fname
            self.ata_tdf = ata_prefix + "tdf/"
            print "TDF related files are expected to be in: " + self.tdf_data_dir
            if path.exists(self.ata_tdf):
                print "The program will try to read TDF saved data from: ./" + self.ata_tdf
            else:
                assert path.exists(self.primitive_tdf) and path.exists(self.primitive_pos), "ERROR: no TDF data given."
            print ""
            
            #weights iterations definitions
            assert (self.nw_loops > 1), "ERROR: wrong number of weight loops"
            try:
                self.weight_tdf=weight_tdf
                self.nw_loops=1
                self.w_incr=0.0
            except:
                self.w_incr = 1.0/(self.nw_loops-1)
                self.weight_tdf=-self.w_incr
            print "Number of weight loops: ", self.nw_loops
            print "with weight increment: ", self.w_incr
            print ""
        elif ("RESP" == self.fit_kind):
            try:
                self.restraints = OrderedDict(restraints)
            except:
                print "ERROR: Asking for RESP but no restraints given"
                raise
            try:
                self.restraint_weight=restraint_weight
            except:
                print "ERROR: Asking for RESP but no weight for restraints given"
                raise
            #Set weight for restraints:
            #if one "float" given, set this value as weight for all restrained atoms
            if type(self.restraint_weight) == type(1.0):
                self.restraint_weight = OrderedDict(
                (k,self.restraint_weight) for k in self.restraints.keys())
            else:
            #check the consistency between restraints and their weights
                self.restraint_weight = OrderedDict(restraint_weight)
                assert set(self.restraint_weight.keys()) == set(self.restraints.keys()), "ERROR: restraints and their weights are not consistent"
            
        elif ("RG" == self.fit_kind):
            try:
                self.rg_weight=rg_weight
                from rappe_goddard import eneg,selci
                self.eneg = eneg
                self.selci = selci
            except:
                print "ERROR: Asking for RG but no weight given or rappe_goddard.py file missing."
                raise
            print "RG weight: ", self.rg_weight
            print "RG restraints will be automatically generated"
            assert eneg[0] > Tiny, "ERROR: values in rappe_goddard.py file seem not to be set."

        print "Atom type and occurrence (optionally restraint and its weight):"
        self.atoms_per_par, self.kind_idx  = self.get_natoms_per_par(self.coor)
        for key, val in self.atoms_per_par.iteritems():
            if ("RESP" == self.fit_kind) and (key in self.restraints):
                print "{:<4s}-->{:>6d}  (restrained to:{:>14.8f}  with weight:{:>14.8f})".format(key, val, self.restraints[key], self.restraint_weight[key])
            else:
                print "{:<4s}-->{:>6d}".format(key, val)
        #check that restraints types actually exist in the system
        if ("RESP" in self.fit_kind):
            for key in self.restraints:
                if key not in self.atoms_per_par:
                    raise Exception ("ERROR: restraint type \"" + key + "\" does not exist (also check blank spaces)")
        self.natoms_per_par = array(self.atoms_per_par.values())
        self.npar = len(self.natoms_per_par)
        self.natoms = sum(self.natoms_per_par)
        
    def get_natoms_per_par(self, coor):
        names = open(coor,"r").readlines()
        if ".pdb" in coor:
            splitted_names=[name.split()[2] for name in names if "ATOM" in name]
        elif ".xyz" in coor: 
            splitted_names=[name.split()[0] for name in names[2:]]

        self.names = splitted_names
        a=""
        types = OrderedDict()
        for n in splitted_names:
             if a!=n:
                 types[n]=splitted_names.count(n)
                 a=n
    
        #generate mask containing the atom indexes grouped by kind
        kind_idx=[]
        for t in types:
            for i,n in enumerate(splitted_names):
                if n == t:
                    kind_idx.append(i+1)
        return types, kind_idx

    def check_types(self, types):
        cnt = Counter(types)
        for at in self.atoms_per_par:
            if self.atoms_per_par[at] != cnt[at]:
                return False
        return True



