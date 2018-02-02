#!/usr/bin/env python
# Driver program for processing Gromacs/CHARMM-specific
# information into PDB and cg_topol data structures
# used by frc_solve.

import sys, os, argparse
from numpy import load, newaxis, sum, zeros
from mol import read_mol
from gromacs.read_top import read_top
from charmm.psf import read_psf
from charmm.read_prm import read_prm
from frc_match import *
from cg_topol import *

from ewsum import ES_seed, ES_frc, dES_frc

# Re-format mol/psf data into PDB : {
#     conn : [Set Int], -- complete connection table,
#     edge : Set (Int, Int), -- unique edges,
#     x : Array Float (n, 3), -- example / reference config.
#     names : [(name : String, res : String, atype : String)],
#     mass : Array Float n,
#     atoms : Int, -- number of atoms
#     L : None | Array Float (3,3),
# }
def pdb_of_mol(mol):
    names = [ (mol.names[i], mol.res[i], mol.t[i]) \
                for i in range(mol.N) ]
    return PDB(names, mol.m, zeros((mol.N,3)),
               set([(i,j) for i,j in mol.bonds]))

# Note: LJPair (as called here) currently uses an LJ-cutoff of 11 Ang.
def topol_of_pdb(pdb, dihedrals, UB=True, LJ=True):
    # Create constraints to fix n (looked up from dihedrals).
    #   If no n, no torsion!
    def constrain_n(type):
        t = tuple(type.split("-"))
        if dihedrals.has_key(t):
            return [u[0] for u in dihedrals[t]]
        print "Warning: no dihedral n found for type:", type
        print("Assuming n = 3!")
        return [3]
    def mk_polytors(*a, **b):
        return PolyTorsion(*a, constrain_n=constrain_n, **b)

    terms = []
    terms.append(bond_terms(pdb, PolyBond))
    terms.append(angle_terms(pdb, PolyAngle))
    if UB:
	terms.append(angle_terms(pdb, PolyUB))
    if isinstance(dihedrals, dict):
        terms.append(torsion_terms(pdb, mk_polytors))
    else:
        terms.append(torsion_terms(pdb, PolyTorsion))
    terms.append(improper_terms(pdb, PolyImprop))
    if LJ == "14":
        # Achtung! Because it creates pdb.pair, this must be called first.
        terms.append(pair_terms(pdb, LJPair, n=5))
        terms.append(pair_n_terms(pdb, LJPair, n=4)) # 1-4 appends to pdb.pair
    elif LJ:
        terms.append(pair_terms(pdb, LJPair, n=4))
    else:
	pair_terms(pdb, LJPair, n=4) # still need to create pairlist
    return FFconcat(terms)

# Creates MQ by giving 1 DOF to ea. unique charge/atom type in the input.
# Then removes one by setting sum = 0.
def wrassle_es(pairs, LJ14, q, t, scale14=0.75):
    mask = []
    for i in range(len(q)-1):
	for j in range(i+1, len(q)):
	    if (i,j) not in pairs:
                if (i,j) in LJ14:
                    mask.append((i,j,scale14))
                else:
                    mask.append((i,j,0.0))

    # Create an extended type bundling charge and type name.
    ext_type = [(qi,ti) for qi,ti in zip(q,t)]
    un = list(set(ext_type))
    mult = array([ext_type.count(u) for u in un])
    MQ = zeros((len(q),len(un)-1))
    for i in range(len(q)):
	j = un.index(ext_type[i])
	if j == len(un)-1: # short straw
	    MQ[i] = -mult[:-1]/float(mult[-1])
	else:
	    MQ[i,j] = 1.0

    print "Using %d minus 1 charge group types:\n"%len(mult)
    print "\n".join("%s %f (%d)"%(ti,qi,m) \
		    for (qi,ti),m in zip(un, mult))
    #print "MQ = " + str(MQ)

    un = array([i for i,j in un[:-1]]) # strip type info.

    return un, mask, MQ

def onefour(tors):
    LJ14=set()
    for i,j,k,l in tors:
        if i<l:
           LJ14.add((i,l))
        else:
           LJ14.add((l,i))
    return LJ14

# reps : t1, t2, oneFour? -> Rmin, eps
def LJ_frc(pairs, LJ14, x, t, reps):
    #print len(tors)
    f = zeros(x.shape)
    for i,j in pairs:
        r = x[:,i] - x[:,j]
        drij = (sum(r*r, -1)**(-0.5))
        rmin, eps = reps(t[i], t[j], (i,j) in LJ14)
        sij = rmin*drij
        r *= (12*eps*(rmin**(-2))*(sij**14 - sij**8))[...,newaxis]
        f[:,i] += r
        f[:,j] -= r
    return f

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Match Them Forces')
    parser.add_argument('mol', metavar='sys.mol', type=str,
			help='System MOL/PSF file.')
    parser.add_argument('xf', metavar='xf.npy', type=str,
			help='Coordinate and force data.')
    parser.add_argument('out', metavar='out_dir', type=str,
			help='Output directory name.')
    parser.add_argument('--top', type=str,
		        help='Read LJ parameters from topol.')
    parser.add_argument('--box', metavar=('Lx', 'Ly', 'Lz'),
			type=float, nargs=3, default=None,
		        help='Box diagonal lengths.')
    parser.add_argument('--tilt', metavar=('Lxy', 'Lxz', 'Lyz'),
			type=float, nargs=3, default=None,
		        help='Box tilt factors.')
    # Toggles.
    parser.add_argument('--noUB', dest='UB', action='store_false',
		        help='Don\'t Fit Urey-Bradley Terms.')
    parser.add_argument('--noLJ', dest='LJ', action='store_false',
		        help='Don\'t Fit LJ parameters at all.')
    parser.add_argument('--chg', action='store_true',
		        help='Fit charges.')
    return parser.parse_args()

def main(args):
    # assumptions:
    top = None
    dih = None
    fudgeQQ = 1.0
    L = None

    if args.top != None:
        if args.top[-4:] == ".prm": # charmm format
            top = read_prm(args.top)
            dih = prm.dihedrals
        else: # gromacs format
            top = read_top(args.top)
            fudgeQQ = top['defaults'][4]
    else:
	top = None

    if args.box != None:
	pdb.L = diag(args.box)
	if args.tilt != None:
	    pdb.L[1,0] = args.tilt[0]
	    pdb.L[2,0] = args.tilt[1]
	    pdb.L[2,1] = args.tilt[2]
    out = args.out

    # Read input data.
    if args.mol[-4:] == '.psf':
        mol = read_psf(args.mol)
    else:
        mol = read_mol(args.mol)
    pdb = pdb_of_mol(mol)
    pdb.L = L
    xf  = load(args.xf)

    # Create topol and FM object.
    topol = topol_of_pdb(pdb, dih, args.UB, args.LJ)
    if args.LJ == False:
	if top == None:
	    raise LookupError, "A parameter file is required for subtracting LJ (using --noLJ)"
	xf[:,1] -= LJ_frc(pdb.pair, onefour(pdb.tors), xf[:,0], mol.t, top.reps)

    q, mask, MQ = wrassle_es(pdb.pair, onefour(pdb.tors), mol.q, \
                             mol.t, fudgeQQ)

    forces = frc_match(topol, pdb, 1.0, 1.0, do_nonlin=args.chg)
    forces.add_nonlin("es", q, ES_seed, ES_frc, dES_frc,
                      (mask, pdb.pair, MQ, pdb.L))
    show_index(forces.topol)

    # Do work.
    forces.append(xf[:,0], xf[:,1])
    forces.dimensionality() # double-checks well-formedness
    forces.maximize()
    forces.write_out(out)

    q = dot(MQ, forces.nonlin["es"][0])
    # Provide a minimal file from which a full gmx parameter set can
    # be written.
    mol.q = q
    mol.write_psf(os.path.join(out, "molecule.psf"))

if __name__=="__main__":
    args = parse_args(sys.argv)
    print(args)
    main(args)

