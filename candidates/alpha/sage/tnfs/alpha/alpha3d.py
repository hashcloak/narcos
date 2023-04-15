"""
This code is a direct translation in Python (.py) of the Magma code (.mag) of cado-nfs
available at http://cado-nfs.gforge.inria.fr, file cado-nfs/nfsd-hd/alpha3d.mag
(c) The Cado-NFS Team
Gnu Lesser General Public License (LGPL) version 2.1
"""

from sage.rings.real_mpfr import RR
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.fast_arith import prime_range
from sage.functions.log import log
from sage.arith.misc import GCD, gcd
from sage.arith.misc import valuation
from sage.misc.prandom import randint

def montecarlo_average_val_affine(f, list_p, N):
    V = [0.0 for i in range(len(list_p))]
    ZZx = f.parent()
    B = N^2
    for i in range(1, N):
        a = [randint(1,B) for i in range(3)]
        while not gcd(a) == 1:
            a = [randint(1,B) for i in range(3)]
        ax = ZZx(a)
        while not ax.is_irreducible():
            a = [randint(1,B) for i in range(3)]
            while not gcd(a) == 1:
                a = [randint(1,B) for i in range(3)]
            ax = ZZx(a)
        R = f.resultant(ax)
        for j in range(len(list_p)):
            V[j] += R.valuation(list_p[j])
    V = [ v/N for v in V]
    return V

# Assume p does not divide lc(f) nor disc(f)
def expect_val_p(f, p):
    ff = f.factor_mod(p)
    nr1 = len([r for r in ff if r[0].degree() == 1])
    nr2 = len([r for r in ff if r[0].degree() == 2])
    val_p = nr1*(p**2+p)/(p**3-1)
    val_p += 2*nr2*p**2/((p**2-1)*(p**2+p+1))
    return val_p

# Main exported function.
def alpha3d(f, pmax):
    Badp = []
    alpha = RR(0.0)
    lc = f.leading_coefficient()
    disc_f = f.discriminant()
    for p in prime_range(pmax+1):
        if (lc % p) == 0 or (disc_f % p) == 0:
            Badp.append(p)
        else:
            alpha += RR(log(p)*(1/(p-1) - expect_val_p(f, p)))

    vals = montecarlo_average_val_affine(f, Badp, 10000)
    alpha += sum([RR(log(Badp[i])*(1/(Badp[i]-1) - vals[i])) for i in range(len(Badp))])
    return float(RR(alpha))
