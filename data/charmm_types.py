# CHARMM forcefield with modifications to select
# UB / torsion constraints using choices from an existing prm file.

# The mk() function wrapper is required, since we can't otherwise
# use "ubs" and "dihedrals" inside of closures.
# I think the issue is how module-level variables are treated specially.
def mk():
    # Read default term choices from cgenff.
    from chemparam import read_prm

    param = read_prm("../data/cgenff.36.prm")
    dihedrals = param.dihedrals
    ubs = param.angles
    nbs = param.nonbonded

    # This code constrains dihedrals based on the 'dihedrals'
    # dict. from a prm file.
    def constrain_n(ty):
        t = tuple(ty.split("_"))[1:]
        if dihedrals.has_key(t):
            return [u[0] for u in dihedrals[t]]
        return None # no constraint
    def ConstrainTors(*a, **b):
        return PolyTorsion(*a, constrain_n=constrain_n, **b)

    def has_ub(ti,tj,tk):
        t = (ti,tj,tk)
        if not ubs.has_key(t): # default = fit the UB term
            return True
        return len(ubs[t]) == 4 # check presence in ubs

    def special14(ti,tj):
        if nbs.has_key(ti) and len(nbs[ti]) == 4:
            return True
        if nbs.has_key(tj) and len(nbs[tj]) == 4:
            return True
        return False

    return ConstrainTors, has_ub, special14

ConstrainTors, has_ub, special14 = mk()

terms = [
    Term("pbond",         PolyBond,      Conn(1,2)),
    Term("pangle",        PolyAngle,     Conn(1,2,3)),
    Term("pub",           PolyUB,        Conn(1,2,3) & TFn(has_ub)),
    Term("ptor",          ConstrainTors, Conn(1,2,3,4)),
    Term("pimproper",     PolyImprop,    OOP()),
    # Note: LJPair (as called here) currently uses an LJ-cutoff of 11 Ang.
    #PairTerm("ljpair_4+", LJPair, Conn(1,2) | Conn(1,3))
    #Term("ljpair_1,4",    LJPair,        Conn(1,4)),
    #PairTerm("ljpair_5+", LJPair, Conn(1,2) | Conn(1,3) | Conn(1,4))
    Term("ljpair_1,4",     LJPair, Conn(1,4) & TFn(special14)),
    PairTerm("ljpair_4,5+",LJPair, Conn(1,2) | Conn(1,3) \
                                   | (Conn(1,4) & TFn(special14)) )
]
