from numpy import sqrt, random, cos, sin, pi, exp

# Perform a random perturbation on the structure.
def kick(x, state, dr=0.08):
    # Select an atom.
    if len(state) == 0: # next random visitation sequence
        state.extend( random.permutation(len(x)).tolist() )
    k = state.pop(-1)

    # Select a 3D vector displacement.
    r = random.standard_normal()*dr
    z = random.random()*2 - 1.0 # z-cosine
    t = random.random()*2*pi # x/y angle
    xy = r*sqrt(1.0 - z*z)
    x[k,0] += xy*cos(t)
    x[k,1] += xy*sin(t)
    x[k,2] += r*z

# Perform 'n' Monte Carlo steps on the geometry
# 1 Hartree = 2625.499638 kJ/mol
# so kT = 2.5 kJ/mol = 1e-3 Har
def mc_geom(x0, en, n, beta=1059.1):
    x0 = x0.copy()

    rng_state = []
    E0 = en(x0) # last accepted energy = E0
    x = x0.copy()

    acc = 0
    for i in xrange(n*len(x0)):
        kick(x, rng_state)
        E = en(x) # sets geom to x
        if E <= E0 or exp(-beta*(E - E0)) > random.random():
            E0 = E
            x0[:] = x
            acc += 1
        else: # reject
            x[:] = x0

        if len(rng_state) == 0: # yield 1 cfg. every pass-through
            yield x, E0

    print("Completed MC: %d/%d accepted moves."%(acc,n*len(x0)))
    # last step was not accepted
    #if abs(E - E0) > 1e-8:
    #    set_geom(db, x0)

