# Wrapper for pybel interface to OpenBabel
# This wrapper provides a molecule class
# that is able to return atom information
# and calculate energies / forces.

import pybel
from numpy import array, zeros, floor

class PyBabel:
    def __init__(self, name):
        mol = pybel.readfile(name[-3:], name)
        for mol in mol: # take only the first conformer
            break
        xyz = mol.write('xyz').split('\n')[2:]

        self.mol = mol.OBMol
        self.atoms = len(mol.atoms)
        self.bonds = []
        for i in xrange(mol.OBMol.NumBonds()):
            bond = mol.OBMol.GetBond(i)
            self.bonds.append( (bond.GetBeginAtomIdx()-1,
                                bond.GetEndAtomIdx()-1) )

        self.elem = [line.split()[0] for line in xyz if line]
        self.x = array([atom.coords for atom in mol.atoms])
        #'acc':atom.IsHbondAcceptor()
        #'don':atom.IsHbondDonor() (heavy)
        #'donH':atom.IsHbondDonorH()

        self.mass = []
        self.res  = []
        for i in range(self.atoms):
            atom = self.mol.GetAtom(i+1)
            self.mass.append( atom.GetExactMass() )
            self.res.append( atom.GetResidue().GetName() + "" )

        self.ff = pybel._forcefields['mmff94']
        if not self.ff.Setup(self.mol):
            raise ValueError, "Error setting up MMFF94 FF."

        self.ff.GetAtomTypes(self.mol)
        self.t = [ mmff_symb[int(self.mol.GetAtom(i+1).\
                             GetData("FFAtomType").GetValue())]
                     for i in range(self.atoms) ]

        self.ff.GetPartialCharges(self.mol)
        chg = [ self.mol.GetAtom(i+1).GetData("FFPartialCharge").GetValue()+"" \
                     for i in range(self.atoms) ]
        self.q = array( map(float, chg) )
        #self.chg = int(floor(reduce(lambda x,y: x+y, self.chg, 0)+1e-8))
        self.chg = int(floor( sum(self.q) + 1e-8 ))

    def push_crd(self): # x -> mol -> ff
        for i in range(self.atoms):
            atom = self.mol.GetAtom(i+1)
            atom.SetVector(self.x[i,0], self.x[i,1], self.x[i,2])
        self.ff.SetCoordinates(self.mol)

    def pull_crd(self): # ff -> mol -> x
        self.ff.GetCoordinates(self.mol)
        self.x = zeros((self.atoms, 3))
        for i in range(self.atoms):
            atom = self.mol.GetAtom(i+1)
            self.x[i,0] = atom.GetX()
            self.x[i,1] = atom.GetY()
            self.x[i,2] = atom.GetZ()

    def minimize(self, steps=1000):
        self.ff.ConjugateGradients(steps)
        self.pull_crd()

    def en(self, x):
        self.x[:] = x
        self.push_crd()
        return self.ff.Energy(False)
    
    def de(self, x):
        return num_diff(self.en, x)

def num_diff(f, x, h = 1e-7):
    ih = 1./(2.0*h)
    h = 0.5/ih

    s = x.shape
    x = x.reshape(-1)
    df = zeros(len(x))
    for i in range(len(x)):
        x0 = x[i]
        x[i] = x0 + h
        df[i] = f(x.reshape(s))
        x[i] = x0 - h
        df[i] -= f(x.reshape(s))
        x[i] = x0

    f0 = f(x.reshape(s))

    df *= ih
    return df.reshape(s)

mmff_symb = {
    1 : 'CR',
    2 : 'C=C',
    3 : 'C=O',
    4 : 'CSP',
    5 : 'HC',
    6 : 'OR',
    7 : 'O=C',
    8 : 'NR',
    9 : 'N=C',
    10 : 'NC=O',
    11 : 'F',
    12 : 'CL',
    13 : 'BR',
    14 : 'I',
    15 : 'S',
    16 : 'S=C',
    17 : 'SO',
    18 : 'SO2',
    19 : 'SI',
    20 : 'CR4R',
    21 : 'HOR',
    22 : 'CR3R',
    23 : 'HNR',
    24 : 'HOCO',
    25 : 'PO4',
    26 : 'P',
    27 : 'HN=C',
    28 : 'HNCO',
    29 : 'HOCC',
    30 : 'CE4R',
    31 : 'HOH',
    32 : 'O2CM',
    33 : 'HOS',
    34 : 'NR+',
    35 : 'OM',
    36 : 'HNR+',
    37 : 'CB',
    38 : 'NPYD',
    39 : 'NPYL',
    40 : 'NC=C',
    41 : 'CO2M',
    42 : 'NSP',
    43 : 'NSO2',
    44 : 'STHI',
    45 : 'NO2',
    46 : 'N=O',
    47 : 'NAZT',
    48 : 'NSO',
    49 : 'O+',
    50 : 'HO+',
    51 : 'O=+',
    52 : 'HO=+',
    53 : '=N=',
    54 : 'N+=C',
    55 : 'NCN+',
    56 : 'NGD+',
    57 : 'CNN+',
    58 : 'NPD+',
    59 : 'OFUR',
    60 : 'C%',
    61 : 'NR%',
    62 : 'NM',
    63 : 'C5A',
    64 : 'C5B',
    65 : 'N5A',
    66 : 'N5B',
    67 : 'N2OX',
    68 : 'N3OX',
    69 : 'NPOX',
    70 : 'OH2',
    71 : 'HS',
    72 : 'SM',
    73 : 'SMO2',
    74 : '=S=O',
    75 : 'P=C',
    76 : 'N5M',
    77 : 'CLO4',
    78 : 'C5',
    79 : 'N5',
    80 : 'CIM+',
    81 : 'NIM+',
    82 : 'N5OX',
    87 : 'FE+2',
    88 : 'FE+3',
    89 : 'F-',
    90 : 'CL-',
    91 : 'BR-',
    92 : 'LI+',
    93 : 'NA+',
    94 : 'K+',
    95 : 'ZN+2',
    96 : 'CA+2',
    97 : 'CU+1',
    98 : 'CU+2',
    99 : 'MG+2'
}
