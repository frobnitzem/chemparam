#!/usr/bin/env python2.7

import sys
from nwchem import *
from nw_helper import *
from numpy import *
import numpy.linalg as la
from babel import PyBabel

def start_sim(elems, x, chg, theory):
    db = nwchem_init(800)
    init_geom(db, elems, x)
    charge(db, "charge %d\n"%chg)
    scf(db, "scf; print none; thresh 1e-5; end\n")
    if theory == "dft":
        init_basis(db, "6-31g*")
        dft(db, "dft; xc b3lyp; noprint; end\n")
    elif theory == "mp2":
        init_basis(db, "cc-pvdz")
        mp2(db, "mp2; noprint; end\n")
    elif theory == "ccsd":
        init_basis(db, "cc-pvdz")
        ccsd(db, "ccsd; noprint; end\n")
    elif theory == "scf":
        init_basis(db, "6-31g*")
    else:
        raise KeyError, "Invalid Theory"
    #prop(db, "property\n%s\nend\n" % "\n".join(
    #    ["DIPOLE", "QUADRUPOLE", "OCTUPOLE",
    #               "vectors %s.mp2nos"%db['file_prefix']]))
    return db

Bohr2Ang = 0.52917721067 # Ang / Bohr
HarBohr2KcalAng = 2625.499638/4.184/Bohr2Ang # kcal/mol-Ang / (Har/Bohr)

# Collect N MC samples of x,f data from QM.
def main(argv):
    assert len(argv) == 4, "Usage: %s <in.mol> <theory> <xf.npy>"%argv[1]

    mol = PyBabel(argv[1])
    theory = argv[2]
    out = argv[3]

    db = start_sim(mol.elem, mol.x, mol.chg, theory)

    xf = load(out)
    for i in range(len(xf)):
        xf[i,1] = -get_de(db, xf[i,0], theory=theory)*HarBohr2KcalAng
    single_op(save)(out, xf)

if __name__=="__main__":
    main(sys.argv)

