# Parser/writer for gromacs top (FF) file.
# #include statements must be removed using "gromacs grompp -pp" option
# or manually pasting includes (maybe cpp)

# As implemented here, the top file contains only forcefield
# information, suitable for including in the ff-topology
# -- not molecule-specific information.
class TOP:
    def __init__(self, defaults={}, atomtypes={}, bondtypes={},
                       constrainttypes={},
                       angletypes={}, dihedraltypes={},
                       pairtypes={}, nonbond_params={}):
        # All these refer to the atomtypes, bondtypes, etc. sections
        # and contain only parameter information -- not actual atom numbers.
        self.defaults = defaults
        self.atoms = atomtypes
        self.bonds = bondtypes
        self.constraints = constrainttypes
        self.angles = angletypes
        self.dihedrals = dihedraltypes
        self.pairs = pairtypes
        self.nbs   = nonbond_params

    def lk_pair(self, t1, t2, coef):
        if coef.has_key((t2, t1)):
            t1, t2 = t2, t1
        v = coef[(t1, t2)]
        return v[0] * 2**(1/6.0), v[1]

    # Returns rmin, eps for a given type-pair.
    def reps(self, t1, t2, onefour):
        rmin, eps = self.lk_pair(t1, t2, self.nbs)
        if onefour:
            try:
                rmin, eps = self.lk_pair(t1, t2, self.pairs)
            except KeyError: # test for gen-pairs
                if self.defaults['default'][2]:
                    eps *= self.defaults['default'][3]
                else:
                    raise
        return rmin, eps

    def write(self, name):
	f = open(name, 'w')
        tf = "no"
        if self.defaults.has_key('default'):
            d = self.defaults['default']
            if d[2]:
                tf = "yes"
            f.write("[ defaults ]\n" \
               +"; nbfunc  comb-rule  gen-pairs  fudgeLJ fudgeQQ\n" \
               +"%5d     %5d      %3s       %.2f %.2f\n" % (
                 d[0], d[1], tf, d[3], d[4]) )
	f.write("\n[ atomtypes ]\n")
	writeem(write_atom, f, self.atoms)
	f.write("\n[ bondtypes ]\n")
	writeem(write_bond, f, self.bonds)
	f.write("\n[ constrainttypes ]\n")
	writeem(write_constraint, f, self.constraints)
	f.write("\n[ angletypes ]\n")
	writeem(write_angle, f, self.angles)
	f.write("\n[ dihedraltypes ]\n")
	writeem(write_dihedral, f, self.dihedrals)
        f.write("\n[ pairtypes ]\n")
	writeem(write_pair, f, self.pairs)
        f.write("\n[ nonbond_params ]\n")
	writeem(write_pair, f, self.nbs)

# ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
def read_default(p, tok):
    b = False
    if tok[2].lower() == 'yes':
        b = True
    p['default'] = (int(tok[0]), int(tok[1]), b, \
                    float(tok[3]), float(tok[4]))

# ; name  atomic_num  mass    charge   ptype          sigma      epsilon
def read_atom(p, tok):
    p[tok[0]] = (int(tok[1]), float(tok[2]), float(tok[3]), \
                 tok[4], float(tok[5]), float(tok[6]))

def write_atom(k,v): # (n, mass, charge, A/D, sigma, eps)
    return "%-6s %d %9.5f %8.3f %s %9.5f %9.5f"%((k,)+v)

# i    j  func       b0          kb
def read_bond(p, tok):
    p[(tok[0],tok[1])] = (int(tok[2]), float(tok[3]), float(tok[4]))

def write_bond(k, v): # t, r0, kb
    return "%-6s %-6s %d %8.5f %6.2f"%(k+v)

def read_constraint(p, tok):
    pass

def write_constraint(k, v): # K, r0
    return ""

# i    j    k  func       th0       cth   [ r0, kUB ]
def read_angle(p, tok):
    p[tuple(tok[:3])] = (int(tok[3]),) + tuple(map(float, tok[4:]))

def write_angle(k, v): # t, theta0, cth, [r0, kUB]
    if len(v) == 5: # UB
        return "%-6s %-6s %-6s 5 %7.3f %7.3f %7.3f %7.3f"%(k+v[1:])
    else:
        return "%-6s %-6s %-6s %d %7.3f %7.3f"%(k+v)

# i    j    k    l   func     coefficients
def read_dihedral(p, tok):
    p[tuple(tok[:4])] = (int(tok[4]),) + tuple(map(float, tok[5:]))

# either/or
# 2 chi0 K
# 3 C0 C1 ...
def write_dihedral(k, v): # [ t, k0, k1, k2, k3, k4, k5 ]
    return "%-6s %-6s %-6s %-6s %d "%(k + (v[0],)) + \
                " ".join("%9.5f"%u for u in v[1:])

# i    j   1  sigma  eps
def read_pair(p, tok):
    assert tok[2] == '1', "Only LJ nonbond_params supported."
    p[tuple(tok[:2])] = tuple(map(float, tok[3:]))

def write_pair(k, v):
    return "%-6s %-6s 1 %8.6f %8.6f"%(k + v)

def writeem(f1, f, kv):
    for key in sorted(kv.keys()):
	f.write(f1(key, kv[key]) + "\n")

# String -> IO(TOP)
def read_top(name):
    # what / how to parse
    parse = { 'defaults':  read_default,
              'atomtypes': read_atom,
              'bondtypes': read_bond,
              'constrainttypes': read_constraint,
              'angletypes': read_angle,
              'dihedraltypes': read_dihedral,
              'pairtypes': read_pair,
              'nonbond_params': read_pair,
            }
    # parse output
    out = {}

    sec = None # current section item parser
    p = {} # current parse

    # lookup parser for new section name
    def parser(sname):
        try:
            sec = parse[sname]
        except KeyError:
            print("Warning! ignoring section [ %s ]"%sname)
            return None, {}
        try:
            p = out[sname]
        except KeyError:
            p = {}
            out[sname] = p
        return sec, p

    with open(name, 'r') as f:
        no = 0 # count current line number
        line = ""
        try:
          for line in f.xreadlines():
            no += 1
            line = line.split(';', 1)[0] # remove comments
            tok = line.split()
            if len(tok) < 1:
                continue
            if tok[0] == '[': # new section
                sname = tok[1].lower()
                sec, p = parser(sname)
                continue
            if sec != None: # parse current line
                sec(p, tok)
        except ValueError:
            raise ValueError, "Parse error on line %d: %s"%(no, line)

    return TOP(**out)

