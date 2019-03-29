#!/usr/bin/env python
# Creates a PQR file with R and eps values hacked
# to work out to atom/OH2 pair values from 'nonbond.itp' when used
# with a specific water model (see fix_lj).

import os, sys
from chemparam import read_mol, read_top
from math import floor
path = os.path.dirname(os.path.realpath(__file__))

# Fix R_min and eps so that Lorentz-Berthelot combination
# rules with SPC water oxygen will generate the "right" answer.
# Assumes units of Angstrom and kcal/mol.
def fix_lj(R, eps, R_w = 3.55332, eps_w = 0.650167284/4.184):
    return 2*R - R_w, eps*eps/eps_w

def main(argv):
    chg = None
    if len(argv) > 3 and argv[1] == '-q':
        argv.pop(1)
        chg = map(float, open(argv.pop(1)).read().split())
    assert len(argv) == 3, "Usage: %s <in.mol> <out.pqre>"%argv[0]

    top = read_top(os.path.join(path, "nonbond.itp"))
    mol = read_mol(argv[1])
    if len(mol.name) > 3:
        res = mol.name[:3].upper()
    else:
        res = mol.name.upper()

    if chg is None:
        chg = mol.q
    else:
        assert len(chg) == mol.N, "Expected %d charges but found %d!"%(mol.N,len(chg))

    # Truncate to 10^-4 and fix rounding.
    Q = ["%7.4f"%q for q in chg]
    
    err = floor(sum(chg)+0.5) - sum(map(float,Q))
    if abs(err) >= 0.0010:
        print("Warning: Non-integer total charge!")
        print("Rounding changed %e -> %.4f"%(sum(chg), sum(map(float,Q))))
    elif abs(err) > 1e-5:
        Q[-1] = "%7.4f"%( float(Q[-1]) + err )

    print mol.t

    with open(argv[2], 'w') as f:
        for i in range(mol.N):
            R, eps = top.reps(mol.t[i], 'OH2', False)
            R, eps = fix_lj(R, eps)
            f.write("ATOM  %5d %4s %4s    1    %8.3f%8.3f%8.3f %s %7.4f %7.4f\n"
                        % (i+1, mol.names[i], res,
                           mol.x[i,0], mol.x[i,1], mol.x[i,2],
                           Q[i], 0.5*R, eps) )

if __name__=="__main__":
    main(sys.argv)
