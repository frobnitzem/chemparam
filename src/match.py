#!/usr/bin/env python2.7
# Driver program for processing Gromacs/CHARMM-specific
# information into PDB and cg_topol data structures
# used by forcesolve.
#
# This program assumes units of:
# en ~ kcal/mol
# dist ~ Angstroms
# charge ~ e_0

import sys, os, argparse
from numpy import load, newaxis, sum, zeros, newaxis
from mol import read_mol
from top import read_top
from psf import read_psf
from prm import read_prm
from forcesolve import *

# Creates MQ by giving 1 DOF to ea. unique charge/atom type in the input.
# Then removes one by setting sum = 0.
def wrassle_es(excl, one4, q, t, scale14=0.75):
    mask = {}
    for i in range(len(q)-1):
	for j in range(i+1, len(q)):
	    if (i,j) in excl:
                mask[(i,j)] = 0.0
            elif (i,j) in one4:
                mask[(i,j)] = scale14

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

    print("Using %d minus 1 charge group types:\n"%len(mult))
    print("\n".join("%s %f (%d)"%(ti,qi,m) \
		    for (qi,ti),m in zip(un, mult)))
    #print("MQ = " + str(MQ))

    un = array([i for i,j in un[:-1]]) # strip type info.

    return un, mask, MQ

# excl : set( (i : int, j : int | 0 <= i < j < N) )
# LJ14 : set( (i : int,  j : int | 0 <= i < j < N) )
# x    : array(float, (N,3))
# t    : [string] | len == N
# reps : (t1 : string, t2 : string, oneFour? : bool) ->
#                                       (Rmin, eps) : (float, float)
def LJ_frc(excl, LJ14, x, t, reps):
    f = zeros(x.shape)
    for i in range(len(t)):
      for j in range(i+1, len(t)):
        if (i,j) in excl:
            continue
        print i,j, (i,j) in LJ14, t[i], t[j]
        r = x[:,i] - x[:,j]
        drij = (sum(r*r, -1)**(-0.5))
        rmin, eps = reps(t[i], t[j], (i,j) in LJ14)
        sij = rmin*drij
        r *= (12*eps*(rmin**(-2))*(sij**14 - sij**8))[...,newaxis]
        f[:,i] += r
        f[:,j] -= r
    return f

# Re-format mol/psf data into PDB : {
#     conn : [Set Int], -- complete connection table,
#     edge : Set (Int, Int), -- unique edges,
#     x : Array Float (n, 3), -- example / reference config.
#     names : [(name : String, res : String, atype : String)],
#     mass : Array Float n,
#     atoms : Int, -- number of atoms
#     L : None | Array Float (3,3),
# }
def pdb_of_mol(mol, L):
    names = [ (mol.names[i], mol.res[i], mol.t[i]) for i in range(mol.N) ]
    pdb = PDB(names, mol.m, mol.x, mol.bonds)
    pdb.L = L
    return pdb

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Match Them Forces')
    parser.add_argument('termfile', metavar='sys.py', type=str,
			help='Python file describing FF terms.')
    parser.add_argument('mol', metavar='sys.mol', type=str,
			help='MOL/PSF file describing molecule bonding.')
    parser.add_argument('xf', metavar='xf.npy', type=str,
			help='Coordinate and force data.')
    parser.add_argument('out', metavar='out_dir', type=str,
			help='Output directory name.')
    parser.add_argument('--param', type=str,
                        help='Read LJ parameters from prm / top file.')
    parser.add_argument('--box', metavar=('Lx', 'Ly', 'Lz'),
			type=float, nargs=3, default=None,
		        help='Box diagonal lengths.')
    parser.add_argument('--tilt', metavar=('Lxy', 'Lxz', 'Lyz'),
			type=float, nargs=3, default=None,
		        help='Box tilt factors.')
    # Toggles.
    parser.add_argument('--noLJ', dest='LJ', action='store_false',
		        help='Don\'t Fit LJ parameters at all.')
    parser.add_argument('--chg', action='store_true',
		        help='Fit charges.')
    return parser.parse_args()

def main(args):
    # assumptions:
    param = None
    fudgeQQ = 1.0
    L = None

    if args.param != None:
        if args.param[-4:] == ".prm": # charmm format
            param = read_prm(args.param)
        else: # gromacs format
            param = read_top(args.param)
            fudgeQQ = param.defaults['default'][4]

    if args.box != None:
	L = diag(args.box)
	if args.tilt != None:
	    L[1,0] = args.tilt[0]
	    L[2,0] = args.tilt[1]
	    L[2,1] = args.tilt[2]
    out = args.out

    # Read input data.
    if args.mol[-3:] == "psf":
        mol = read_psf(args.mol)
    elif args.mol[-3:] == "itp":
        mol = mol_of_itp( read_itp(args.mol) )
    else:
        mol = read_mol(args.mol)
    xf  = load(args.xf)

    # Create topol and FM object.
    pdb   = pdb_of_mol(mol, L)
    topol = read_terms(pdb, args.termfile)
    excl  = set(srt2(i-1,j-1) for (i,j),t in (Conn(1,2)|Conn(1,3)).run(pdb))
    one4  = set(srt2(i-1,j-1) for (i,j),t in (Conn(1,4)).run(pdb))
    if args.LJ == False:
        assert pdb.L == None, "Periodicity not supported when subtracting LJ."
        if param == None:
            raise LookupError, "A parameter file is required for subtracting LJ (using --noLJ)."
        xf[:,1] -= LJ_frc(excl, one4, xf[:,0], mol.t, param.reps)

    q, mask, MQ = wrassle_es(excl, one4, mol.q, mol.t, fudgeQQ)

    forces = frc_match(topol, pdb, 1.0, 1.0, do_nonlin=args.chg)
    # assumes kcal/mol energy units and e_0 chg. units
    forces.add_nonlin("es", q, ES_seed, ES_frc, dES_frc, (mask, MQ, L))
    show_index(forces.topol)

    # Do work.
    forces.append(xf[:,0], xf[:,1])
    #forces.dimensionality(1e-7) # double-checks well-formedness
    forces.maximize()
    forces.write_out(out)

    q = dot(MQ, forces.nonlin["es"][0])
    # Provide a minimal file from which a full gmx parameter set can
    # be written.
    mol.q = q
    mol.write_psf(os.path.join(out, "molecule.psf"))
    #mol.write_itp(os.path.join(out, "molecule.itp"), \
    #              os.path.basename(args.mol).rsplit('.',1)[0])

if __name__=="__main__":
    args = parse_args(sys.argv)
    print(args)
    main(args)

