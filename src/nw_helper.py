from nwchem import *
import scipy.optimize
from numpy import allclose, reshape, arange, array, ndarray, sqrt, random, cos, sin, pi, exp, zeros, diag
fmin = scipy.optimize.fmin_bfgs

def init_bq(db, X, Q):
    bq(db, """bq units angstrom
  clear
%s
end"""%'\n'.join("  %f %f %f %f"%(x[0],x[1],x[2],q) for x, q in zip(X, Q)))

def init_geom(db, names, x):
    db.names = names
    geom(db, """geometry units angstrom nocenter noautoz noprint
  symmetry c1
%s
end
"""%'\n'.join("  %3s %f %f %f"%(n,y[0],y[1],y[2]) for n,y in zip(names, x)))

def init_basis(db, bname):
    if "cc" in bname: # use spherical harmonics for correlation-consistent basis
        bas(db,  "basis spherical; * library %s; end;\n"%bname)
    else:
        bas(db,  "basis; * library %s; end;\n"%bname)

def get_en(db, x=None, theory="scf"):
    if isinstance(x, ndarray):
        set_geom(db, x)
    nwtask(db, "%s energy"%theory)
    return db["%s:energy"%theory][0]

def get_de(db, x=None, theory="scf"):
    if isinstance(x, ndarray):
        print x.shape
        set_geom(db, x)
    nwtask(db, "%s gradient"%theory)
    return db["%s:gradient"%theory].reshape((-1, 3))

def hessian(db, thr=1e-6):
    #hess(db, "hessian; thresh %e; print nucdd_cont; end\n"%thr)
    hess(db, "hessian; thresh %e; end\n"%thr)
    #with open(db['task:hessian file name']) as f:
    with open('perm/nyah.hess') as f:
        L = array( map(lambda s: float(s.replace('D','E')), f.readlines()) )

    # Parse to a full matrix
    m = -0.5 + sqrt(0.25 + len(L)*2)
    n = int(m+1e-8)
    if abs(n - m) > 1e-7:
        raise ValueError, "Invalid return len from Hessian: %d"%len(L)
    H = zeros((n, n))
    k = 0
    for i in range(n):
        j = i+1
        H[i,:j] = L[k:k+j]
        k += j
    H += H.transpose() - diag(diag(H))
    return H #.reshape((n,3,n,3))

def set_geom(db, x):
    #if allclose(db["geometry:geometry:coords"], reshape(x, -1)):
    #    return
	# cause of "cphf_solve:SCF residual greater than 1d-2"?"

    #n = len(x)
    #x = db["geometry:geometry:coords"]
    #x[:3*n] = reshape(x, -1)
    #db["geometry:geometry:coords"] = x
    init_geom(db, db.names, x)

def print_val(v):
	print "Current Value: %f"%v

# minimizes? without crashing on failure
def min_geom(db, x0, **opts):
    s = x0.shape
    x = fmin(lambda x: get_en(db, x.reshape(s), **opts), x0.reshape(-1), \
             fprime = lambda x: get_de(db, x.reshape(s), **opts).reshape(-1), \
             maxiter = 20, callback=print_val).reshape(s)
    set_geom(db, x)
    return x

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

