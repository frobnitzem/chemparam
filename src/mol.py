# Methods for creating a molecule topology from a structure.
from babel import PyBabel
from numpy import array, ndarray
from datetime import datetime
from forcesolve import modprod
from itp import ITP

# assign a name to the atom in the format: %s%d % type, num
def name_atom(t, naming):
    if naming.has_key(t):
        naming[t] += 1
        i = naming[t]
    else:
        naming[t] = 1
        i = 1
    return "%s%d"%(t,i)

# name : String / filename -> IO(Mol)
def read_mol(name):
    mol = PyBabel(name) # returns generator of mols

    naming = {} # special dict to help with naming
    atoms = []
    for i in range(mol.atoms):
        t = mol.t[i] # mmff94 type name
        name = name_atom(mol.elem[i], naming)
        res = mol.res[i]
        #resn = atom.GetResidue().GetNum()
        atoms.append({ 'name': name,
                       'res': res,
                       't': t,
                       'q': mol.q[i],
                       'm': mol.mass[i],
                       'x': mol.x[i]
                    })
    return Mol(atoms, mol.bonds)

# {k:[v]}, n -> [{k:v}] (length n)
def atomize(dl, n):
    return [ dict([(k, dl[k][i]) for k,v in dl.iteritems()]) \
             for i in range(n) ]
# [{k:v}], keys -> {k,[v]}
def dvector(ld, keys):
    return dict([k, [u[k] for u in ld]] for k in keys)

def mol_of_itp(itp):
    return Mol( itp.atoms, array( itp.bonds.keys() ), \
                itp.moleculetype['name'] )

# positions are in units of Angstroms
class Mol:
    def __init__(self, atoms, bonds, name='chemicalX'):
        self.N = len(atoms)
        self.bonds = bonds
        self.name = name

        self.x     = array([a['x'] for a in atoms])
        self.q     = array([a['q'] for a in atoms])
        self.m     = array([a['m'] for a in atoms])

        self.t     = [a['t'] for a in atoms]
        self.res   = [a['res'] for a in atoms]
        self.names = [a['name'] for a in atoms]
        if len(atoms) > 0 and atoms[0].has_key('chain') \
                    and len(atoms[0]['chain'].strip()) > 0:
            self.chain = [a['chain'] for a in atoms]
        else:
            self.chain = ['A' for a in atoms]
        if len(atoms) > 0 and atoms[0].has_key('resn'):
            self.resn = [a['resn'] for a in atoms]
        else:
            self.resn = [1 for a in atoms]

    def write_itp(self, name, title=None):
        itp = self.to_itp(self.name)
        if title is not None:
            itp.moleculetype['name'] = title
        itp.write(name)

    def to_itp(self, name, title=None):
        if title is None:
            title = self.name
        G = [set() for i in range(self.N)]
        for b in self.bonds: # generate node-based connection table
            G[b[0]].add(b[1])
            G[b[1]].add(b[0])
        ang = list_angles(G)
        tor = list_dihedrals(self.bonds, G)
        imp = list_impropers(G)
        dihedrals = dict([t, 3] for t in tor)
        dihedrals.update( dict([t, 2] for t in imp) )
        pair = onefour(tor)
        return ITP( moleculetype={'name':title, 'nrexcl':3},
                    atoms = atomize({
                        'n': range(1,self.N+1),
                        'chain': self.chain,
                        'resn': self.resn,
                        'res': self.res,
                        'name': self.names,
                        't': self.t,
                        'q': self.q,
                        'm': self.m
                       }, self.N),
                    bonds=dict([(i,j), 1] for i,j in self.bonds),
                    angles=dict([t, 5] for t in ang),
                    dihedrals=dihedrals,
                    pairs=dict([(i,j), 1] for i,j in pair)
                  )

    def write_atom_psf(self, i):
        return "%8d     %-4s  %5d %4s     %-4s      %4s     %13.6f %9.5f"%(
                  i+1, self.chain[i], self.resn[i], self.res[i], self.names[i],
                  self.t[i], self.q[i], self.m[i]) \
             + "            0   0.00000     -0.301140E-02"
    def write_psf(self, name):
        d = datetime.utcnow()
        f = open(name, 'w')
        f.write("PSF EXT CMAP CHEQ XPLOR\n\n")
        f.write("         3 !NTITLE\n")
        f.write("* Generated by the ChemParameter Server\n")
        f.write("* David M. Rogers (based on ForceSolve, arXiV:1003.4741v1)\n")
        f.write("*  DATE:    %02d/%02d/%02d     %02d:%02d:%02dZ"%(
                d.month, d.day, d.year, d.hour, d.minute, d.second ) \
                + "     CREATED BY USER: davidrogers@usf.edu\n")
        f.write("\n%8d !NATOM\n"%self.N)
        [ f.write(self.write_atom_psf(i)+"\n") for i in range(self.N) ]
        G = [set() for i in range(self.N)]
        for b in self.bonds: # generate node-based connection table
            G[b[0]].add(b[1])
            G[b[1]].add(b[0])
        f.write("\n%8d !NBOND: bonds\n"%len(self.bonds))
        f.write(blocked(array(self.bonds)+1))
        ang = list_angles(G)
        f.write("\n%8d !NTHETA: angles\n"%len(ang))
        f.write(blocked(array(ang)+1, 9))
        tor = list_dihedrals(self.bonds, G)
        f.write("\n%8d !NPHI: dihedrals\n"%len(tor))
        f.write(blocked(array(tor)+1))
        imp = list_impropers(G)
        f.write("\n%8d !NIMPHI: impropers\n"%len(imp))
        f.write(blocked(array(imp)+1))
        f.write(vestigal_tail%(blocked([[0]]*self.N), \
                               blocked([[1]]*self.N)))
        f.close()

def onefour(tors):
    LJ14=set()
    for i,j,k,l in tors:
        if i<l:
           LJ14.add((i,l))
        else:
           LJ14.add((l,i))
    return LJ14

def list_angles(G):
    angles = []
    for j in range(len(G)):
        angles += [(i,j,k) for i,k in modprod(G[j], G[j]) if i < k]
    angles.sort()
    return angles

def list_dihedrals(b, G):
    tors = []
    for j,k in b:
        tors += [(i,j,k,l) for i,l in modprod(
                                       G[j]-set([k]), \
                                       G[k]-set([j]) ) if i != l]
    for i,t in enumerate(tors):
        if t[0] > t[3]:
            tors[i] = t[3], t[2], t[1], t[0]
    tors.sort()
    return tors

def list_impropers(G):
    return [(i,)+tuple(b) for i,b in enumerate(G) if len(b) == 3]

def concat(x):
    if len(x) == 0:
        return x
    if isinstance(x, ndarray):
        return x.reshape(-1)
    if isinstance(x[0], tuple):
        return sum(x, ())
    return sum(x, [])

def blocked(l, n=8):
    if len(l) == 0:
        return "\n"
    v = concat(l)
    u = (len(v)-n-1) % n + 1
    lines = [("%8d"*n)%tuple(v[i:i+n]) for i in range(0,len(v)-n,n)]
    lines.append( ("%8d"*u)%tuple(v[-u:]) )
    return "\n".join(lines) + "\n" 

vestigal_tail = '''
       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB

%s
       1       0 !NGRP NST2
       0       1       0

       1 !MOLNT
%s
       0       0 !NUMLP NUMLPH

       0 !NCRTERM: cross-terms
'''

if __name__=="__main__":
    import sys
    mol = read_mol(sys.argv[1])
    print mol.names
    print mol.t
    print mol.q
    print mol.m
    print mol.x.shape
    print len(mol.bonds) + 1 - mol.N

