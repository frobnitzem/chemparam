# Methods for creating a molecule topology from a structure.
from babel import PyBabel
from numpy import array

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

# positions are in units of Angstroms
class Mol:
    def __init__(self, atoms, bonds):
        self.N = len(atoms)
        self.bonds = bonds

        self.x     = array([a['x'] for a in atoms])
        self.q     = array([a['q'] for a in atoms])
        self.m     = array([a['m'] for a in atoms])

        self.t     = [a['t'] for a in atoms]
        self.res   = [a['res'] for a in atoms]
        self.names = [a['name'] for a in atoms]

if __name__=="__main__":
    import sys
    mol = read_mol(sys.argv[1])
    print mol.names
    print mol.t
    print mol.q
    print mol.m
    print mol.x.shape
    print len(mol.bonds) + 1 - mol.N

