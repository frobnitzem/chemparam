#!/usr/bin/env python2.7
# J. Am. Chem. Soc., Vol. 114, No. 20, 1992.

from ctypes import *

calc = CDLL("compute.so")
#double compute_lr(double ai, double Ni, double Ai, double Gi,
#                  double aj, double Nj, double Aj, double Gj,
#                  int da, double *r_ij);
calc.compute_lr.argtypes = [c_double]*8 + [c_int, POINTER(c_double)]
calc.compute_lr.restype  = c_double

# type  alpha-i     N-i       A-i       G-i DA Symb   Origin
coef = """
    1     1.050     2.490     3.890     1.282 - CR     E94
    2     1.350     2.490     3.890     1.282 - C=C    E94
    3     1.100     2.490     3.890     1.282 - C=O    E94
    4     1.300     2.490     3.890     1.282 - CSP    E94
    5     0.250     0.800     4.200     1.209 - HC     C94
    6     0.70      3.150     3.890     1.282 A OR     C94
    7     0.65      3.150     3.890     1.282 A O=C    C94
    8     1.15      2.820     3.890     1.282 A NR     C94
    9     0.90      2.820     3.890     1.282 A N=C    C94
   10     1.000     2.820     3.890     1.282 A NC=O   E94
   11     0.35      3.480     3.890     1.282 A F      C94
   12     2.300     5.100     3.320     1.345 A CL     E94
   13     3.400     6.000     3.190     1.359 A BR     E94
   14     5.500     6.950     3.080     1.404 A I      E94
   15     3.00      4.800     3.320     1.345 A S      C94
   16     3.900     4.800     3.320     1.345 A S=C    E94
   17     2.700     4.800     3.320     1.345 - SO     E94
   18     2.100     4.800     3.320     1.345 - SO2    E94
   19     4.500     4.200     3.320     1.345 - SI     E94
   20     1.050     2.490     3.890     1.282 - CR4R   E94
   21     0.150     0.800     4.200     1.209 D HOR    C94
   22     1.100     2.490     3.890     1.282 - CR3R   E94
   23     0.150     0.800     4.200     1.209 D HNR    C94
   24     0.150     0.800     4.200     1.209 D HOCO   C94
   25     1.600     4.500     3.320     1.345 - PO4    E94
   26     3.600     4.500     3.320     1.345 A P      E94
   27     0.150     0.800     4.200     1.209 D HN=C   C94
   28     0.150     0.800     4.200     1.209 D HNCO   C94
   29     0.150     0.800     4.200     1.209 D HOCC   C94
   30     1.350     2.490     3.890     1.282 - CE4R   E94
   31     0.150     0.800     4.200     1.209 D HOH    C94
   32     0.75      3.150     3.890     1.282 A O2CM   C94
   33     0.150     0.800     4.200     1.209 D HOS    C94
   34     1.00      2.820     3.890     1.282 - NR+    C94
   35     1.50      3.150     3.890     1.282 A OM     X94
   36     0.150     0.800     4.200     1.209 D HNR+   C94
   37     1.350     2.490     3.890     1.282 - CB     E94
   38     0.85      2.820     3.890     1.282 A NPYD   C94
   39     1.10      2.820     3.890     1.282 - NPYL   C94
   40     1.00      2.820     3.890     1.282 A NC=C   E94
   41     1.100     2.490     3.890     1.282 - CO2M   C94
   42     1.000     2.820     3.890     1.282 A NSP    E94
   43     1.000     2.820     3.890     1.282 A NSO2   E94
   44     3.00      4.800     3.320     1.345 A STHI   C94
   45     1.150     2.820     3.890     1.282 - NO2    E94
   46     1.300     2.820     3.890     1.282 - N=O    E94
   47     1.000     2.820     3.890     1.282 A NAZT   X94
   48     1.200     2.820     3.890     1.282 A NSO    X94
   49     1.00      3.150     3.890     1.282 - O+     X94
   50     0.150     0.800     4.200     1.209 D HO+    C94
   51     0.400     3.150     3.890     1.282 - O=+    E94
   52     0.150     0.800     4.200     1.209 D HO=+   C94
   53     1.000     2.820     3.890     1.282 - =N=    X94
   54     1.30      2.820     3.890     1.282 - N+=C   C94
   55     0.80      2.820     3.890     1.282 - NCN+   E94
   56     0.80      2.820     3.890     1.282 - NGD+   E94
   57     1.000     2.490     3.890     1.282 - CNN+   E94
   58     0.80      2.820     3.890     1.282 - NPD+   E94
   59     0.65      3.150     3.890     1.282 A OFUR   C94
   60     1.800     2.490     3.890     1.282 A C%     E94
   61     0.800     2.820     3.890     1.282 A NR%    E94
   62     1.300     2.820     3.890     1.282 A NM     X94
   63     1.350     2.490     3.890     1.282 - C5A    E94
   64     1.350     2.490     3.890     1.282 - C5B    E94
   65     1.000     2.820     3.890     1.282 A N5A    E94
   66     0.75      2.820     3.890     1.282 A N5B    C94
   67     0.950     2.82      3.890     1.282 A N2OX   X94
   68     0.90      2.82      3.890     1.282 A N3OX   C94
   69     0.950     2.82      3.890     1.282 A NPOX   C94
   70     0.87      3.150     3.890     1.282 A OH2    C94
   71     0.150     0.800     4.200     1.209 D HS     C94
   72     4.000     4.800     3.320     1.345 A SM     X94
   73     3.000     4.800     3.320     1.345 - SMO2   X94
   74     3.000     4.800     3.320     1.345 - =S=O   X94
   75     4.000     4.500     3.320     1.345 A P=C    X94
   76     1.200     2.820     3.890     1.282 A N5M    X94
   77     1.500     5.100     3.320     1.345 A CLO4   X94
   78     1.350     2.490     3.890     1.282 - C5     X94
   79     1.000     2.820     3.890     1.282 A N5     X94
   80     1.000     2.490     3.890     1.282 - CIM+   C94
   81     0.80      2.820     3.890     1.282 - NIM+   C94
   82     0.950     2.82      3.890     1.282 A N5OX   X94
   87     0.45      6.        4.        1.4   - FE+2   X94
   88     0.55      6.        4.        1.4   - FE+3   X94
   89     1.4       3.48      3.890     1.282 A F-     X94
   90     4.5       5.100     3.320     1.345 A CL-    X94
   91     6.0       6.000     3.190     1.359 A BR-    X94
   92     0.15      2.        4.        1.3   - LI+    X94
   93     0.4       3.5       4.        1.3   - NA+    X94
   94     1.0       5.        4.        1.3   - K+     X94
   95     0.43      6.        4.        1.4   - ZN+2   X94
   96     0.9       5.        4.        1.4   - CA+2   X94
   97     0.35      6.        4.        1.4   - CU+1   X94
   98     0.40      6.        4.        1.4   - CU+2   X94
   99     0.35      3.5       4.        1.3   - MG+2   X94
"""

lk_da = {'-': 2, 'D': 1, 'A': -1}

# Manually translated symbols.
mmff_symb = {
        1:  ('CR',     'ALKYL CARBON, SP3'),
        2:  ('CSP2',   ' VINYLIC CARBON or GENERIC SP2 CARBON'),
        3:  ('C=',     'CARBONYL CARBON'),
        4:  ('=C=',    'ALLENIC CARBON'),
        5:  ('HC',     'H ATTACHED TO C/SI'),
        6:  ('OR',     'ALCOHOL OR ETHER OR GENERAL DIVALENT OXYGEN'),
        7:  ('O=',     'CARBONYL OXYGEN'),
        8:  ('NR',     'NITROGEN IN ALIPHATIC AMINES/IMINES'),
        9:  ('N=N',    'NITROGEN IN AZO/IMINE'),
        10: ('NC=O',   'NITROGEN IN AMIDES'),
        11: ('F',      'FLUORINE'),
        12: ('CL',     'CHLORINE'),
        13: ('BR',     'BROMINE'),
        14: ('I',      'IODINE'),
        15: ('S',      'SULFUR IN THIOETHERS AND MERCAPTANS'),
        16: ('S=C',    'TERMINAL SULFUR DOUBLY BONDED TO CARBON'),
        17: ('S=O',    'SULFUR IN SULFOXIDES'),
        18: ('SO2',    'SULFUR IN SULFONES/SULFATE'),
        19: ('Si',     'SILICON'),
        20: ('CR4R',   'CARBON IN 4-MEMBERED RINGS'),
        21: ('HOR',    'HYDROGEN IN ALCOHOL/HYDROXIDE ANION'),
        22: ('CR3R',   'CARBON IN A 3-MEMBERED RING'),
        23: ('HNR',    'H-N'),
        24: ('HOCO',   'H-O IN OXYLIC ACIDS'),
        25: ('PO',      'TETRACOORDINATE PHOSPHORUS'),
        26: ('P',      'TRICOORDINATE P, AS IN PHOSPHINES'),
        27: ('HN=N',   'AZO/IMINE HYDROGEN'),
        28: ('HSP2',   'GENERAL H ON SP2 NITROGEN'),
        29: ('HOCC',   'H-O IN ENOLS AND PHENOLS/HO-C=N'),
        30: ('CE4R',   'OLEFINIC CARBON IN 4-MEMBERED RINGS'),
        31: ('HOH',    'HYDROGEN IN H2O'),
        32: ('OX',     'TERMINAL ANIONIC OXYGEN'),
        33: ('HOS',    'H ON OXYGEN ATTACHED TO SULFUR'),
        34: ('NR+',    'QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED'),
        35: ('OM',     'OXIDE OXYGEN, NEGATIVELY CHARGED'),
        36: ('HN+',    'H ON N5+, N5A+ OR N5B+'),
        37: ('CB',     'CARBON AS IN BENZENE, PYRROLE'),
        38: ('NPYD',   'NITROGEN, AS IN PYRIDINE'),
        39: ('NPYL',   'NITROGEN, AS IN PYRROLE'),
        40: ('NC=R',   'NITROGEN ON N-C=R or C-C TRIPLE BOND'),
        41: ('CO2M',   'CARBOXYLATE ANION CARBON'),
        42: ('NSP',    'NITROGEN, TRIPLE BONDED'),
        43: ('N%',     'NITROGEN ATTACHED TO PI-GROUP'),
        44: ('STHI',   'SULFUR AS IN THIOPHENE'),
        45: ('NO',     'NITRO/NITRATE GROUP NITROGEN'),
        46: ('N=O',    'NITROSO NITROGEN'),
        47: ('NAZT',   'TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP'),
        48: ('NSO',    'DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP'),
        49: ('O+',     'POSITIVELY CHARGED OXONIUM (TRICOORDINATE) OXYGEN'),
        50: ('HO+',    'HYDROGEN ON O+ OXYGEN'),
        51: ('O=+',    'POSITIVELY CHARGED OXENIUM (DICOORDINATE) OXYGEN'),
        52: ('HO=+',   'HYDROGEN ON OXENIUM OXYGEN'),
        53: ('=N=',    'NITROGEN IN C=N=N OR -N=N=N '),
        54: ('N+=R',   'POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N'),
        55: ('NCN+',   'N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2'),
        56: ('NGD+',   'GUANIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1/3'),
        57: ('CNN+',   'C IN +N=C-N RESONANCE STRUCTURES'),
        58: ('NPD+',   'PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1'),
        59: ('OFUR',   'AROMATIC OXYGEN AS IN FURAN'),
        60: ('C%',     'ISONITRILE CARBON'),
        61: ('NR%',    'ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]'),
        62: ('NM',     'DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1'),
        63: ('C5A',    'ALPHA CARBON IN 5-MEMBERED HETEROAROMATIC RING'),
        64: ('C5B',    'BETA CARBON IN 5-MEMBERED HETEROAROMATIC RING'),
        65: ('N5A',    'ALPHA AROM HETEROCYCLIC 5-RING  NITROGEN'),
        66: ('N5B',    'BETA AROM HETEROCYCLIC 5-RING  NITROGEN'),
        67: ('N2OX',   'SP2-HYDRIDIZED N-OXIDE NITROGEN'),
        68: ('N3OX',   'SP3-HYDRIDIZED N-OXIDE NITROGEN'),
        69: ('NPOX',   'PYRIDINE N-OXIDE NITROGEN'),
        70: ('OH2',    'OXYGEN ON WATER'),
        71: ('HSP',    'H ATTACHED TO DIVALENT, DICOORDINATE S/P'),
        72: ('ST',     'TERMINAL SULFUR BONDED TO PI-GROUP'),
        73: ('SO2M',   'SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP'),
        74: ('=S=O',   'SULFINYL SULFUR, EG. IN C=S=O'),
        75: ('P=C',    'PHOSPHOROUS DOUBLY BONDED TO CARBON'),
        76: ('N5M',    'NEGATIVELY CHARGED N IN, E.G, TRI- OR TETRAZOLE ANION'),
        77: ('CLO4',   'CHLORINE IN PERCHLORATE ANION, CLO4(-)'),
        78: ('C5',     'GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING'),
        79: ('N5',     'GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING  '),
        80: ('CIM+',   'C IN N-C-N IN IMIDAZOLIUM ION'),
        81: ('N5+',    'POSITIVE N5 NITROGEN'),
        82: ('N5OX',   'N-OXIDE NITROGEN IN GENERAL 5-RING POSITION'),
        87: ('FE+2',   'IRON +2 CATION'),
        88: ('FE+3',   'IROM +3 CATION'),
        89: ('F-',     'FLUORIDE ANION'),
        90: ('CL-',    'CHLORIDE ANION'),
        91: ('BR-',    'BROMIDE ANION'),
        92: ('LI+',    'LITHIUM CATION'),
        93: ('NA+',    'SODIUM CATION'),
        94: ('K+',     'POTASSIUM CATION'),
        95: ('ZINC',   'DIPOSITIVE ZINC'),
        95: ('ZN+2',   'DIPOSITIVE ZINC'),
        96: ('CA+2',   'DIPOSITIVE CALCIUM'),
        97: ('CU+1',   'MONOPOSITIVE COPPER'),
        98: ('CU+2',   'DIPOSITIVE COPPER'),
        99: ('MG+2',   'DIPOSITIVE MAGNESIUM CATION')
    }

# Generate an LJ table using the parameters above.
def mk_table():
    lines = [line.split() for line in coef.split('\n') if len(line) >= 5]

    symb = dict([(int(tok[0]), tok[-2]) for tok in lines])
    t = [tok[-2] for tok in lines]
    # alpha-i     N-i       A-i       G-i DA Symb   Origin
    anag = [map(float,tok[1:5])+[lk_da[tok[5]]] for tok in lines]

    #print len(t) - len(set(t))

    rij = c_double(0.0)
    print("[ nonbond_params ]")
    print("; i j   func sigma (nm) eps (kJ/mol)")
    for i,n in enumerate(t):
        I = anag[i]
        for j in range(i, len(t)):
            m = t[j]
            J = anag[j]
            eps = calc.compute_lr(I[0], I[1], I[2], I[3],
                                  J[0], J[1], J[2], J[3],
                                  I[4]+J[4], byref(rij))
            #print n, m, rij.value, eps
            # convert to sigma (nm)
            r = rij.value * 0.1*2**(-1/6.0)
            print("%-6s %-6s 1 %8.6f %8.6f"%(n, m, r, eps))

    print("atom_types = {")
    for tok in lines:
        print("    %s : '%s',"%(tok[0], tok[-2]))
    print("}")

mk_table()

