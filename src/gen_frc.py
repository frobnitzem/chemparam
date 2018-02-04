#!/usr/bin/env python2.7

import sys
#from nwchem import *
from nw_helper import *
from numpy import *
import numpy.linalg as la
from babel import PyBabel

def start_sim(elems, x, chg):
    db = nwchem_init(1500)
    init_geom(db, elems, x)
    charge(db, "charge %d\n"%chg)
    scf(db, "scf; print none; thresh 1e-5; end\n")
    dft(db, "dft; xc b3lyp; noprint; end\n")
    mp2(db, "mp2; tight; end\n")
    init_basis(db, "6-31g*")
    #prop(db, "property\n%s\nend\n" % "\n".join(
    #    ["DIPOLE", "QUADRUPOLE", "OCTUPOLE",
    #               "vectors %s.mp2nos"%db['file_prefix']]))
    return db

Bohr = 0.52917721067 # Ang / Bohr
fac  = 2625.499638/4.184/Bohr # kcal/mol-Ang / (Har/Bohr)

# Collect N MC samples of x,f data from QM.
def main(argv):
    assert len(argv) == 5, "Usage: %s <mol.sdf> <N> <theory> <xf.npy>"%argv[1]

    mol = PyBabel(argv[1])
    mol.minimize()
    N = int(argv[2])
    theory = argv[3]
    out = argv[4]

    #db = start_sim(mol.elem, mol.x, mol.chg)

    xf = zeros((N,2,mol.atoms,3))
    i = 0
    for x, E in mc_geom(mol.x, mol.en, N, beta=1./0.6):
        xf[i,0] = x*0.1 # nm
        xf[i,1] = -mol.de(x) # kcal/Ang
        #xf[i,1] = -get_de(db, x, theory=theory)*fac
        i += 1
        if i % 100 == 0:
            print("Step %d: en = %f"%(i,E))

    save(out, xf)
    #import code
    #code.interact(local=locals())

if __name__=="__main__":
    main(sys.argv)

