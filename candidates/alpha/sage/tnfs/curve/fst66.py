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
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0

# all except 0 mod 18
allowed_k = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17, 19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35, 37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53]

def polynomial_params(k):
    """Returns the parameters of Construction FST 6.6 (Freeman-Scott-Teske) where (k % 18) != 0

    INPUT:
    - `k`: embedding degree k not a multiple of 18

    RETURN: univariate polynomials px, rx, cx, tx, yx, betax, lambx, and D=3

    test:
    from tnfs.curve.fst66 import polynomial_params
    for k in [ki for ki in range(6, 50, 1) if (ki % 18) != 0]:
        px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
        x = px.variables()[0]
        print("k={} lcm(denom(px)) = {}".format(k, lcm([ci.denom() for ci in px.list()]) ))
        i = (k%2)+1
        j = 3-i
        print("k={} lcm(denom(px(3*x+{}))) = {}".format(k, i, lcm([ci.denom() for ci in (px(3*x+i)).list()]) ))
        print("k={} lcm(denom(rx(3*x+{}))) = {}".format(k, i, lcm([ci.denom() for ci in (rx(3*x+i)).list()]) ))
        print("k={} lcm(denom(px(3*x+{}))) = {}".format(k, j, lcm([ci.denom() for ci in (px(3*x+j)).list()]) ))
        print("k={} lcm(denom(rx(3*x+{}))) = {}".format(k, j, lcm([ci.denom() for ci in (rx(3*x+j)).list()]) ))

    """

    if (k % 18) == 0:
        raise ValueError("Error k cannot be a multiple of 18, but k={}={}*18 given".format(k, k/18))
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    Phi_k = QQx(cyclotomic_polynomial(k))
    D = 3
    m = 3
    u_mod_m = [(k % 2) + 1]
    if (k%6) == 1:
        rx = QQx(cyclotomic_polynomial(6*k))
        tx = QQx(-x**(k+1)+x+1)
        yx = QQx((-x**(k+1) + 2*x**k -x-1)/3)
        px = QQx((x+1)**2*(x**(2*k)-x**k+1)- 3*x**(2*k+1))/QQ(3)
        lambx = QQx(x**k-1) # the other one is -x**k
    elif (k%6) == 2:
        rx = QQx(cyclotomic_polynomial(3*k))
        tx = QQx(x**(k//2+1)-x+1)
        yx = QQx((x**(k//2+1)+2*x**(k//2)+x-1)/3)
        px = QQx((x-1)**2*(x**k-x**(k//2)+1)+3*x**(k+1))/QQ(3)
        lambx = QQx(x**(k//2)-1)
    elif (k % 18) == 3:
        rx = QQx(cyclotomic_polynomial(2*k))
        tx = QQx(x**(k//3+1)+1)
        yx = QQx((-x**(k//3+1)+2*x**(k//3)+2*x-1)/3)
        px = QQx((x**2-x+1)*(x**(2*k//3)-x**(k//3)+1)+3*x**(k//3+1))/QQ(3)
        lambx = QQx(x**(k//3)-1)
    elif (k%18) == 9 or (k%18) == 15:
        rx = QQx(cyclotomic_polynomial(2*k))
        tx = QQx(-x**(k//3+1)+x+1)
        yx = QQx((-x**(k//3+1)+2*x**(k//3)-x-1)/3)
        px = QQx((x+1)**2*(x**(2*k//3)-x**(k//3)+1)-3*x**(2*k//3+1))/QQ(3)
        lambx = QQx(x**(k//3)-1)
    elif (k%6) == 4:
        rx = QQx(cyclotomic_polynomial(3*k))
        tx = QQx(x**3+1)
        yx = QQx((x**3-1)*(2*x**(k//2)-1)/3)
        px = QQx((x**3-1)**2*(x**k-x**(k//2)+1) + 3*x**3)/QQ(3)
        lambx = QQx(x**(k//2)-1)
    elif (k%6) == 5:
        rx = QQx(cyclotomic_polynomial(6*k))
        tx = QQx(x**(k+1)+1)
        yx = QQx((-x**(k+1)+2*x**k+2*x-1)/3)
        px = QQx((x**2-x+1)*(x**(2*k)-x**k+1) + 3*x**(k+1))/QQ(3)
        lambx = QQx(x**k-1)
    elif (k%6) == 0:
        rx = QQx(cyclotomic_polynomial(k))
        tx = QQx(x+1)
        yx = QQx((x-1)*(2*x**(k//6)-1)/3)
        px = QQx((x-1)**2*(x**(k//3)-x**(k//6)+1) + 3*x)/QQ(3)
        lambx = QQx(x**(k//6)-1)
    # yx = (t-2)/sqrt(-D) % rx
    assert ((lambx**2 + lambx + 1) % rx) == 0
    if k==6:
        tx = x+1
        yx = (tx-2)*(2*x-1)/3 # no mod rx
        px = (tx**2 + 3*yx**2)/4
    #else:
    #    yx = QQx((tx-2)*(2*lambx+1) % rx)/QQ(3)
    #    tx = tx % rx
    px_denom = QQ(3) #Integer(lcm([ci.denom() for ci in px.list()]))
    assert ((px+1-tx) % rx) == 0, "Error (px+1-tx)%rx != 0, px={}*({}), tx={}, rx={}, (px+1-tx)={}".format(px_denom,px_denom*px, tx,rx, (px+1-tx).factor())
    cx = (px+1-tx) // rx
    assert (Phi_k(px) % rx) == 0, "Error k={}, Phi_k(px)%rx != 0, = {}".format(k, Phi_k(px) % rx)
    assert px.is_irreducible(), "Error px is not irreducible, px= {}".format(px.factor())
    assert rx.is_irreducible(), "Error rx is not irreducible, rx= {}".format(rx.factor())

    assert 4*px == (tx**2 + D*yx**2), "Error 4*px != tx^2+D*yx^2, {}*px = {}, tx={}, D={}, yx={}, tx^2+Dyx^2 = {}".format(px_denom, px_denom*px,tx,D,yx,tx**2+D*yx**2)
    # now computes beta s.t. (-1+sqrt(-3))/2 = beta mod px
    # computes sqrt(-3) mod px
    # 4*px = tx^2 + 3*yx^2 <=> tx^2 = -3*yx^2 <=> tx/yx = +/- sqrt(-3)
    # computes the inverse of yx
    g,u,v = px.xgcd(yx)
    if g == 1:
        inv_yx = QQx(v)
    else:
        inv_yx = QQx(v)/QQ(g)
    betax = (QQx(-1+inv_yx*tx)/2) % px
    if betax.leading_coefficient() < 0:
        betax = -betax-1
    #print("computed betax = {}".format(betax))
    
    assert ((betax**2 + betax + 1) % px) == 0
    assert ((lambx**2 + lambx + 1) % rx) == 0

    return px, rx, tx, cx, yx, betax, lambx, D

def get_m_u_mod_m(k):
    """Return the congruence conditions m, u_mod_m on the seed x = u mod m
    so that the parameters are all integers and p(x), r(x) generate primes

    INPUT: embedding degree k
    RETURN: modulus m (int), list of residues u mod m, and m=1, u_mod_m=[0] if no condition is required
    """
    m = 3
    u_mod_m = [(k % 2) + 1]
    #if (k % 2) == 0: # k = 0, 2, 4 mod 6
    #    u_mod_m = [1]
    #else:            # k = 1, 3, 5 mod 6
    #    u_mod_m = [2]
    return m, u_mod_m

def coeffs_params(k):
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    Px_denom = Integer(lcm([ci.denom() for ci in px.list()]))
    Px = [Integer(ci) for ci in (Px_denom*px).list()]
    Rx_denom = Integer(lcm([ci.denom() for ci in rx.list()]))
    Rx = [Integer(ci) for ci in (Rx_denom*rx).list()]
    Tx_denom = Integer(lcm([ci.denom() for ci in tx.list()]))
    Tx = [Integer(ci) for ci in (Tx_denom*tx).list()]
    Cx_denom = Integer(lcm([ci.denom() for ci in cx.list()]))
    Cx = [Integer(ci) for ci in (Cx_denom*cx).list()]
    Yx_denom = Integer(lcm([ci.denom() for ci in yx.list()]))
    Yx = [Integer(ci) for ci in (Yx_denom*yx).list()]
    BETAx_denom = Integer(lcm([ci.denom() for ci in betax.list()]))
    BETAx = [Integer(ci) for ci in (BETAx_denom*betax).list()]
    LAMBx_denom = Integer(lcm([ci.denom() for ci in lambx.list()]))
    LAMBx = [Integer(ci) for ci in (LAMBx_denom*lambx).list()]
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D

def poly_cofactor_gt(k):
    """Computes the co-factors for GT: Phi_k(p(x))/r(x)

    It is always irreducible:
    from tnfs.curve.fst66 import poly_cofactor_gt
    for k in [ki for ki in range(6, 53, 1) if (ki % 18) != 0]:
        print("\nFST66-k{}".format(k))
        c = poly_cofactor_gt(k)
        print(["(deg {})^{} ".format(ci.degree(), ei) if ei>1 else "(deg {})".format(ci.degree()) for (ci, ei) in c.factor()])
        x = c.variables()[0]
        i = (k%2)+1
        print(lcm([ci.denom() for ci in (c(3*x+i)).list()]))
    """
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    cx = cyclotomic_polynomial(k)(px)
    assert (cx % rx) == 0
    cx = cx // rx
    return cx

class FST66(BrezingWeng):
    """A FST66 pairing-friendly curve of embedding degree k != 0 mod 18
    """
    # re-use the init function of class BrezingWeng
    def __init__(self, k, u, b=None, cofactor_r=1, verbose_init=False):
        px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, D = coeffs_params(k)
        super().__init__(k, D, u, px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, 0, b, cofactor_r, verbose_init)
        # b=0
    def _repr_(self):
        return "FST66_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (pairing-friendly curve k="+str(self._k)+") with seed "+str(self._u)+"\n"+super(BrezingWeng, self)._repr_()
