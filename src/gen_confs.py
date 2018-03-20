#!/usr/bin/env python2.7

import sys
from numpy import *
from geom import *
import numpy.linalg as la
from babel import PyBabel

# Collect N MC samples of x,f data from MMFF94.
# Used for conformer generation.
def main(argv):
    assert len(argv) == 4, "Usage: %s <in.mol> <N> <xf.npy>"%argv[1]

    mol = PyBabel(argv[1])
    mol.minimize()
    N = int(argv[2])
    out = argv[3]

    xf = zeros((N,2,mol.atoms,3))
    i = 0
    skip = min(N, 100)
    for x, E in mc_geom(mol.x, mol.en, skip+N, beta=0.33/0.6):
        if i >= skip:
            xf[i-skip,0] = x # Ang
            xf[i-skip,1] = -mol.de(x) # kcal/Ang
        i += 1
        if i % 100 == 0:
            print("Step %d: en = %f"%(i,E))

    save(out, xf)

if __name__=="__main__":
    main(sys.argv)

