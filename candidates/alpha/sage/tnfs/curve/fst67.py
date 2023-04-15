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
from sage.rings.rational_field import QQ

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import BrezingWeng

def polynomial_params(k):
    """
    Returns the parameters px,rx,cx,tx,yx,betax,lambx, and D=2 of Construction 6.7
    From Taxonomy, Freeman, Scott, Teske

    INPUT:
    - `k`: embedding degree, k=0 mod 3

    test:
    from tnfs.curve.fst67 import polynomial_params
    for k in [ki for ki in range(6, 53, 3)]:
        px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
        x = px.variables()[0]
        print("k={} lcm(denom(px)) = {}".format(k, lcm([ci.denom() for ci in px.list()]) ))
        print("k={} lcm(denom(px(2*x))) = {}".format(k, lcm([ci.denom() for ci in (px(2*x)).list()]) ))
        print("k={} lcm(denom(px(2*x+1))) = {}".format(k, lcm([ci.denom() for ci in (px(2*x+1)).list()]) ))
        if (k % 24) == 0:
            print("k={} lcm(denom(px(4*x+3))) = {}".format(k, lcm([ci.denom() for ci in (px(4*x+3)).list()]) ))
            print("k={} lcm(denom(px(4*x+1))) = {}".format(k, lcm([ci.denom() for ci in (px(4*x+1)).list()]) ))
    """
    if (k%3) != 0:
        raise ValueError("Construction 6.7, 3 | k required but k={}={} mod 3 given".format(k,k%3))
    m = 2
    u_mod_m = [1]
    if (k % 24) == 0:
        m = 4
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    l = lcm(8,k)
    rx = QQx(cyclotomic_polynomial(l))
    tx = QQx(x**(l // k) + 1)
    px = QQx((2*(x**(l//k)+1)**2 + (1-x**(l//k))**2*(x**(5*l//24) + x**(l//8) - x**(l//24))**2))/8
    yx = QQx((1-x**(l//k))*(x**(5*l//24) + x**(l//8) - x**(l//24)))/2
    tx = tx % rx
    assert ((px+1-tx) % rx) == 0
    assert (cyclotomic_polynomial(k)(px) % rx) == 0
    assert px.is_irreducible()
    assert rx.is_irreducible()
    cx = (px+1-tx) // rx
    D = 2
    assert px == (tx**2 + D*yx**2)/4
    g,u,v = px.xgcd(tx)
    if g == 1:
        inv_tx = QQx(v)
    else:
        inv_tx = QQx(v)/QQ(g)
    # px = (tx^2 + 2*yx^2) / 4
    # tx^2 = -2*yx^2
    # -1/2 = (yx/tx)^2
    # -2 = 4*(yx/tx)^2
    # -2 = (2*yx/tx)^2
    # or: tx/yx = sqrt(-2)
    betax = QQx((inv_tx*yx*2) % px)
    # x^(l//8) is a 8-th root of unity <-> (1+i)/sqrt(2)
    # x^(l//4) is a 4-th root of unity <-> i
    # (1-i)*(1+i)/sqrt(2) = 2/sqrt(2) = sqrt(2)
    # i*(1-i)*(1+i)/sqrt(2) = i*2/sqrt(2) = sqrt(-2)
    lambx = QQx(x**(l//4) * x**(l//8) * (1-x**(l//4))) % rx
    assert ((betax**2 + 2) % px) == 0
    assert ((lambx**2 + 2) % rx) == 0
    return px, rx, tx, cx, yx, betax, lambx, D

def get_m_u_mod_m(k):
    """Return the congruence conditions m, u_mod_m on the seed x = u mod m
    so that the parameters are all integers and p(x), r(x) generate primes

    INPUT: None
    RETURN: modulus m (int), list of residues u mod m, and m=1, u_mod_m=[0] if no condition is required
    """
    m = 2
    if (k % 24) == 0:
        m = 4
    u_mod_m = [1]
    return m, u_mod_m

def coeffs_params(k):
    """ Return the coefficients in ZZ of the polynomials together with the denominators """

    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    Px_denom = lcm([fi.denom() for fi in px.list()])
    Px = [Integer(fi) for fi in (Px_denom * px).list()]
    Rx_denom = lcm([fi.denom() for fi in rx.list()])
    Rx = [Integer(fi) for fi in (Rx_denom * rx).list()]
    Tx_denom = lcm([fi.denom() for fi in tx.list()])
    Tx = [Integer(fi) for fi in (Tx_denom * tx).list()]
    Cx_denom = lcm([fi.denom() for fi in cx.list()])
    Cx = [Integer(fi) for fi in (Cx_denom * cx).list()]
    Yx_denom = lcm([fi.denom() for fi in yx.list()])
    Yx = [Integer(fi) for fi in (Yx_denom * yx).list()]

    Betax_denom = lcm([fi.denom() for fi in betax.list()])
    Betax = [Integer(fi) for fi in (Betax_denom * betax).list()]

    Lambx_denom = lcm([fi.denom() for fi in lambx.list()])
    Lambx = [Integer(fi) for fi in (Lambx_denom * lambx).list()]

    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, Betax, Betax_denom, Lambx, Lambx_denom, D

def poly_cofactor_gt(k):
    """Computes the co-factors for GT: Phi_k(p(x))/r(x)

    It is always irreducible:
    from tnfs.curve.fst67 import poly_cofactor_gt
    for k in [ki for ki in range(6, 53, 3)]:
        print("\nFST67-k{}".format(k))
        c = poly_cofactor_gt(k)
        print(["(deg {})^{} ".format(ci.degree(), ei) if ei>1 else "(deg {})".format(ci.degree()) for (ci, ei) in c.factor()])
        x = c.variables()[0]
        print(lcm([ci.denom() for ci in (c(2*x+1)).list()]))
    """
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    ex = cyclotomic_polynomial(k)(px)
    assert (ex % rx) == 0
    ex = ex // rx
    return ex

class FST67(BrezingWeng):
    """A FST67 pairing-friendly curve of embedding degree k = 0 mod 3 and D=2
    """
    # re-use the init function of class BrezingWeng
    def __init__(self, k, u, a=None, b=None, cofactor_r=1, verbose_init=False):
        """
        curve parameters a,b are optional as compute_a_b() knows how to handle D=2 (j=8000): (a,b) = (-30, 56) or a quadratic twist.
        """
        px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, D = coeffs_params(k)
        super().__init__(k, D, u, px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, a, b, cofactor_r, verbose_init)
        # a*b != 0, j=8000 s D=2
    def _repr_(self):
        return "FST67_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (pairing-friendly curve k="+str(self._k)+") with seed "+str(self._u)+"\n"+super(BrezingWeng, self)._repr_()
