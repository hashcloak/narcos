import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.arith.misc import XGCD, xgcd
from sage.arith.functions import lcm
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import BrezingWeng
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_a_j1728, get_curve_generator_order_r_j1728

def polynomial_params(k):
    """Returns the parameters px,rx,cx,tx,yx,betax,lambx, and D=4 of Construction 6.2

    INPUT:
    - `k`: embedding degree, k=1 mod 2

    test:
    from tnfs.curve.fst62 import polynomial_params
    for k in [ki for ki in range(3, 30, 2) if (ki % 2) == 1]:
        px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
        x = px.variables()[0]
        print("k={} lcm(denom(px)) = {}".format(k, lcm([ci.denom() for ci in px.list()]) ))
        print("k={} lcm(denom(px(2*x))) = {}".format(k, lcm([ci.denom() for ci in (px(2*x)).list()]) ))
        print("k={} lcm(denom(px(2*x+1))) = {}".format(k, lcm([ci.denom() for ci in (px(2*x+1)).list()]) ))

    """
    if (k % 2) != 1:
        raise ValueError("Error k should be odd, but k={}=0 mod 2".format(k))
    m = 2
    u_mod_m = [1]
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    rx = QQx(cyclotomic_polynomial(4*k))
    # x is a (4*k)-th root of unity, x^2 is a (2*k)-th root of unity, and x^4 is a k-th root of unity
    # zeta_k = -x^2 because zeta_k^k = (-1)^k*x^(2*k) = -x^(2*k) = 1 because x^(2*k) = -1 (indeed, x^(4*k) = 1)
    tx = QQx(-x**2+1)
    px = QQx(x**(2*k+4) + 2*x**(2*k+2) + x**(2*k) + x**4 -2*x**2 + 1)/4
    yx = QQx(x**(k+2)+x**k)
    assert ((px+1-tx) % rx) == 0
    cx = (px+1-tx) // rx
    assert (cyclotomic_polynomial(k)(px) % rx) == 0
    assert px.is_irreducible()
    assert rx.is_irreducible()
    D = 1
    assert px == (tx**2 + D*yx**2)/4
    g,u,v = px.xgcd(yx)
    if g == 1:
        inv_yx = QQx(v)
    else:
        inv_yx = QQx(v)/QQ(g)
    # px = (tx^2 + D*yx^2) / 4
    # tx^2 = -D*yx^2
    # -D = (tx/yx)^2
    betax = QQx((inv_yx*tx) % px)
    #print("computed betax = {}".format(betax))
    lambx = QQx(x**k)
    assert ((betax**2 + D) % px) == 0
    assert ((lambx**2 + D) % rx) == 0
    return px, rx, tx, cx, yx, betax, lambx, D

def get_m_u_mod_m():
    """Return the congruence conditions m, u_mod_m on the seed x = u mod m
    so that the parameters are all integers and p(x), r(x) generate primes

    INPUT: None
    RETURN: modulus m (int), list of residues u mod m, and m=1, u_mod_m=[0] if no condition is required
    """
    m = 2
    u_mod_m = [1]
    return m, u_mod_m

def coeffs_params(k):
    """ Return the coefficients of the polynomials """
    if (k % 2) != 1:
        raise ValueError("Error k should be odd, but k={}=0 mod 2".format(k))
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    #Rx = cyclotomic_polynomial(4*k).list()
    Rx = [Integer(ri) for ri in rx.list()]
    Rx_denom = 1
    Tx = [1,0,-1] # -x^2+1
    Tx_denom = 1
    Px = [0]*(2*k+4 +1)
    Px[0] = 1
    Px[2] = -2
    Px[4] = 1
    Px[2*k] = 1
    Px[2*k+2] = 2
    Px[2*k+4] = 1
    Px_denom = 4
    C0x = cx.list()
    Cx_denom = lcm([ci.denom() for ci in C0x])
    Cx = [Integer(Cx_denom*ci) for ci in C0x]
    # tx^2 - 4*p = x^4 -2*x^2+1 - 4*px = -(x^(2*k+4) + 2*x^(2*k+2) + x^(2*k))
    # = -(x^(k+2)+x^k)^2
    Yx = [0]*(k+2 +1)
    Yx[k] = 1
    Yx[k+2] = 1
    Yx_denom = 2
    # px+1-tx = px+(4+4*x^2-4)
    # betax = a fourth root of unity modulo px, sqrt(-1) = (tx-2)/yx mod px
    # px = ((x^(k+2) + x^k)^2 + (x^2-1)^2)/4
    # computes the denominator of betax
    Bx = betax.list()
    BETAx_denom = lcm([ci.denom() for ci in Bx])
    BETAx = [Integer(BETAx_denom*ci) for ci in Bx]
    LAMBx = [0]*(k+1) ; LAMBx[k] = 1 # x^k
    LAMBx_denom = 1
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D

def poly_cofactor_gt(k):
    """Computes the co-factors for GT: Phi_k(p(x))/r(x)

    It is always irreducible:
    from tnfs.curve.fst62 import poly_cofactor_gt
    for k in [ki for ki in range(3, 50, 2) if (ki % 2) == 1]:
        print("\nFST62-k{}".format(k))
        c = poly_cofactor_gt(k)
        print(["(deg {})^{} ".format(ci.degree(), ei) if ei>1 else "(deg {})".format(ci.degree()) for (ci, ei) in c.factor()])
        x = c.variables()[0]
        print(lcm([ci.denom() for ci in (c(2*x+1)).list()]))
    """
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    cx = cyclotomic_polynomial(k)(px)
    assert (cx % rx) == 0
    cx = cx // rx
    return cx

class FST62(BrezingWeng):
    """A FST62 pairing-friendly curve of embedding degree k odd
    """
    # re-use the init function of class BrezingWeng
    def __init__(self, k, u, a=None, cofactor_r=1, verbose_init=False):
        px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, D = coeffs_params(k)
        super().__init__(k, D, u, px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, a, None, cofactor_r, verbose_init)
        # b=None
    def _repr_(self):
        return "FST62_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (pairing-friendly curve k="+str(self._k)+") with seed "+str(self._u)+"\n"+super(BrezingWeng, self)._repr_()
