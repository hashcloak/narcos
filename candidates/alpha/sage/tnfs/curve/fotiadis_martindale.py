import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil, sqrt
from sage.arith.misc import XGCD, xgcd, gcd
from sage.arith.functions import lcm
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.number_field.number_field import NumberField
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs
import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_a_j1728, get_curve_generator_order_r_j1728
from tnfs.curve.pairing_friendly_curve import compute_beta_lambda
from tnfs.curve.pairing_friendly_curve import BrezingWeng

allowed_code = [1, 10, 17, 19, 20, 23, 24, 25, 26, 27, 28]

def polynomial_params_from_code(code):
    
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    if code == 1:
        k=8
        D=1
        m=2
        u_mod_m = [0]
        zeta_k_mod_rx = x
        exp_tr = 1
        rx = x**4+1
        tx = x+1+rx
        yx = x**3-x**2
        px = (x**8+x**6+5*x**4+x**2+4*x+4)/4
        cx = (px+1-tx) // rx
        betax = (x**6 - x**5 - 2*x**3 + 3*x**2 - 3*x - 2)/4
        lambx = x**2
        
    elif code == 10:
        k = 9
        D = 3
        m = 3
        u_mod_m = [2]
        zeta_k_mod_rx = x
        exp_tr = 5
        rx = (x**6 + x**3 + 1)
        tx = x**5 + 1 + rx
        yx = (3*x**6 + x**5 + 5*x**3 + 2*x**2 + 4)/3
        px = (3*x**12 + 3*x**11 + x**10 + 9*x**9 + 7*x**8 + x**7 + 16*x**6 + 10*x**5 + x**4 + 13*x**3 + 4*x**2 + 7)/3
        #rx /= QQ(3) # because of u=1 mod 3
        cx = (px+1-tx) // rx
        betax = (129*x**11 - 147*x**10 - 164*x**9 + 103*x**8 - 456*x**7 - 493*x**6 + 80*x**5 - 565*x**4 - 547*x**3 - 198*x**2 - 266*x - 421)/325
        lambx = x**3

    elif code == 17:
        k=12
        D=3
        m=1
        u_mod_m = [0]
        zeta_k_mod_rx = 6*x**2
        exp_tr = 7
        # rx = (cyclotomic_polynomial(12))(6*x**2).factor()[1][0]
        # rx = rx * lcm([ri.denom() for ri in rx.list()])
        #tx = (zeta_k_mod_rx**exp_tr % rx) + 1
        #K.<w> = NumberField(rx)
        #yx = (tx-2)* ((K(-1/D).sqrt()).polynomial()) % rx
        #px = (tx**2 + D*yx**2)/4
        
        rx = 36*x**4 + 36*x**3 +18*x**2 + 6*x + 1 
        tx = (-6*x**2) + 1
        yx = 48*x**3 + 30*x**2 + 12*x + 3
        px = 1728*x**6 + 2160*x**5 + 1548*x**4 + 756*x**3 + 240*x**2 + 54*x + 7
        cx = (px+1-tx) // rx
        betax = 4320*x**5 + 3672*x**4 + 2430*x**3 + 954*x**2 + 237*x + 46
        lambx = 36*x**3 + 18*x**2 + 6*x + 1

    elif code == 19:
        k=12
        D=3
        m=30
        u_mod_m = [8,23]
        rx = (x**4 -2*x**3 -3*x**2 + 4*x + 13)/225 # because when 
        tx = (-x**3+4*x**2+5*x+6)/15
        yx = (x**3 - 4*x**2 + 5*x + 4)/15
        px = (x**6-8*x**5+21*x**4-17*x**3+13*x**2+45*x+21)/225
        cx = (px+1-tx) // rx
        zeta_k_mod_rx = tx-1
        exp_tr = 1
        lambx = (2*x**3 - 3*x**2 - 7)/15
        betax = (x**5 - 9*x**4 + 30*x**3 - 41*x**2 + 30*x + 15)/30

    elif code == 20:
        k=12
        D=3
        m=285
        u_mod_m = [209,266]
        rx = (x**4 -37*x**2 + 361)/9025 # because for u=209, 266 mod 285, p(u)  in ZZ and 9025 | r(u)
        tx = (-2*x**3 + 17*x + 95)/95
        yx = (8*x**3 + 38*x**2 - 163*x - 703)/285
        px = (x**6 + 8*x**5 - 18*x**4 - 326*x**3 - 342*x**2 + 3143*x + 6859)/1425
        cx = (px+1-tx) // rx
        zeta_k_mod_rx = tx-1
        exp_tr = 1
        betax = (5*x**5 + 61*x**4 + 173*x**3 - 832*x**2 - 5133*x - 7309)/85
        lambx = (x**2 - 21)/5

    elif code==23:
        k=16
        D=1
        m=2
        u_mod_m = [0]
        rx = x**8 + 1
        tx0 = x+1
        ht = 1
        tx = x**8 + x + 2 # x+1+rx = tx0 + ht*rx
        assert tx == tx0 + ht*rx
        yx = (x-1)*x**4
        px = (x**16 + x**10 + 5*x**8 + x**2 + 4*x + 4)/4
        cx = x**2 * (x**2 + 1) * (x**4 - x**2 + 1)/4
        zeta_k_mod_rx = (tx % rx)-1
        exp_tr = 1
        betax = (x**15 + 2*x**14 + 4*x**12 - 4*x**11 - 4*x**10 - 3*x**9 - 2*x**8 + x**7 + 10*x**6 - 8*x**5 + 12*x**4 - 12*x**3 - 12*x**2 - 11*x - 6)/16 # % px
        lambx = x**4 # % rx

    elif code==24:
        k=16
        D=1
        m=70
        u_mod_m = [10,60]
        rx = x**8 + 48*x**4 + 625
        tx0 = (2*x**5 + 41*x + 35)/35
        ht = 1
        tx = (35*x**8 + 2*x**5 + 1680*x**4 + 41*x + 21910)/35
        assert tx == tx0 + ht*rx
        yx = (x**5 + 5*x**4 + 38*x + 120)/35
        px = (245*x**16 + 28*x**13 + 23520*x**12 + x**10 + 1920*x**9 + 871225*x**8 + 48*x**6 + 45204*x**5 + 14723760*x**4 + 625*x**2 + 361148*x + 96012500)/980
        cx = (245*x**8 + 28*x**5 + 11760*x**4 + x**2 + 576*x + 152640)/980
        zeta_k_mod_rx = (tx % rx)-1
        exp_tr = 1
        betax = (41043090655*x**15 - 153792893765*x**14 + 105849347975*x**13 + 584915058007*x**12 + 3330730991539*x**11 - 10873914745505*x**10 + 5270358381374*x**9 + 45992762538358*x**8 + 91046244341433*x**7 - 264042638124995*x**6 + 75476902659297*x**5 + 1225652143320001*x**4 + 853163126296693*x**3 - 2189535628180375*x**2 + 119023447836160*x + 11207576998933556)/7817303287808 # % px
        lambx = (x**4 + 24)/7 # % rx

    elif code==25:
        k=18
        D=3
        m=3
        u_mod_m = [1]
        px = (3*x**12 - 3*x**9 + x**8 - 2*x**7 + 7*x**6 - x**5 - x**4 - 4*x**3 + x**2 - 2*x + 4)/3
        rx = x**6 - x**3 + 1
        tx0 = -x**4 + 1
        ht = 1
        tx = x**6 - x**4 - x**3 + 2
        assert tx == tx0 + ht*rx
        cx = (3*x**6 + x**2 - 2*x + 1)/3
        yx = (3*x**6 + x**4 - x**3 - 2*x + 2)/3
        zeta_k_mod_rx = x
        exp_tr = 13
        lambx = x**6
        betax = (-93*x**11 + 48*x**10 + 240*x**9 + 267*x**8 - 235*x**7 + 84*x**6 - 139*x**5 + 191*x**4 + 131*x**3 + 266*x**2 - 411*x - 112)/342

    elif code==26:
        k=18
        D=3
        m=3
        u_mod_m = [1]
        px = (21*x**12 - 6*x**10 + 1533*x**9 + x**8 - 334*x**7 + 42007*x**6 + 37*x**5 - 6199*x**4 + 512092*x**3 + 343*x**2 - 38368*x + 2343376)/21
        rx = x**6 + 37*x**3 + 343
        tx0 = (x**4 + 16*x + 7)/7
        ht = 1
        tx = (7*x**6 + x**4 + 259*x**3 + 16*x + 2408)/7
        assert tx == tx0 + ht*rx
        cx = (21*x**6 - 6*x**4 + 756*x**3 + x**2 - 112*x + 6811)/21
        yx = (21*x**6 - 5*x**4 + 763*x**3 - 94*x + 6944)/21
        zeta_k_mod_rx = (tx % rx)-1
        exp_tr = 1
        lambx = x**3 + 18
        betax = (-987945*x**11 + 5571804*x**10 + 16703934*x**9 - 54894753*x**8 + 304751239*x**7 + 893075946*x**6 - 1017840467*x**5 + 5563708427*x**4 + 15919441663*x**3 - 6297737416*x**2 + 33910970205*x + 94588671572)/22136568

    elif code==27:
        k=20
        D=1
        m=2
        u_mod_m = [1]
        rx = x**8 - x**6 + x**4 - x**2 + 1 # cyclotomic_polynomial(20)
        tx = x + 1
        px = (x**12 - 2*x**11 + x**10 + x**2 +2*x +1)/4
        cx = (x-1)**2*(x**2+1)/4
        yx = (x-1) * x**5
        zeta_k_mod_rx = tx-1
        exp_tr = 1
        betax = (x**11 - 3*x**10 + 4*x**9 - 4*x**8 + 4*x**7 - 4*x**6 + 2*x**5 + x + 1)/2 # mod px
        lambx = x**5 # mod rx

    elif code==28: # new one from Georgios Fotiadis
        k = 16
        D = 1
        m = 2
        u_mod_m = [0]
        rx = x**8+1
        tx0 = x**5+1
        ht = 1
        tx = x**8+x**5+2
        assert tx == tx0 + ht*rx
        px = (x**16 + 2*x**13 + x**10 + 5*x**8 + 6*x**5 + x**2 + 4)/4
        cx = (x*(x**3 + 1)/2)**2
        yx = x*(x**3 + 1)
        zeta_k_mod_rx = x
        exp_tr = 5
        betax = (x**12 + x**9 + 3*x**4 + x)/2 # mod px
        lambx = x**4 # mod rx

    else:
        raise ValueError("Error, FotiadisMartindale curves can have code in {} but code {} was provided".format(allowed_code, code))

    Phi_k = QQx(cyclotomic_polynomial(k))
    assert (Phi_k(tx-1) % rx) == 0, "Error k={}, Phi_k(tx-1)%rx != 0, = {}, tx=({})^{}+1 mod rx = {}, rx={}".format(k, Phi_k(tx-1) % rx, zeta_k_mod_rx, exp_tr*m, tx, rx)
    assert (Phi_k(px) % rx) == 0, "Error k={}, Phi_k(px)%rx != 0, = {}".format(k, Phi_k(px) % rx)
    assert px.is_irreducible(), "Error px is not irreducible, px= {}".format(px.factor())
    assert rx.is_irreducible(), "Error rx is not irreducible, rx= {}".format(rx.factor())
    # now computes beta s.t. (-1+sqrt(-3))/2 = beta mod px
    # computes sqrt(-3) mod px
    # 4*px = tx^2 + 3*yx^2 <=> tx^2 = -3*yx^2 <=> tx/yx = +/- sqrt(-3)
    # computes the inverse of yx
    """
    g,u,v = px.xgcd(yx)
    if g == 1:
        inv_yx = QQx(v)
    else:
        inv_yx = QQx(v)/QQ(g)
    if D==3:
        betax = (QQx(-1+inv_yx*tx)/2) % px
        lambx = QQx(zeta_k_mod_rx**(l//3)) % rx
    elif D==1:
        betax = (QQx(inv_yx*tx)) % px
        lambx = QQx(zeta_k_mod_rx**(l//4)) % rx
    if betax.leading_coefficient() < 0:
        betax = -betax-1
    if lambx.leading_coefficient() < 0:
        lambx = -lambx-1
    """
    # so that if zeta_k_mod_rx == x and 3|k, then l=k and l//3 == k//3,
    # otherwise l=3*k and x^k is a 3-rd root of unity mod rx=Phi_l(x)
    if D==3:
        assert ((betax**2 + betax + 1) % px) == 0
        assert ((lambx**2 + lambx + 1) % rx) == 0
    elif D==1:
        assert ((betax**2 + 1) % px) == 0
        assert ((lambx**2 + 1) % rx) == 0
        
    return px, rx, tx, cx, yx, betax, lambx, D, k, m, u_mod_m

def coeffs_params(code):
    px, rx, tx, cx, yx, betax, lambx, D, k, m, u_mod_m = polynomial_params_from_code(code)

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
    if betax != 0 and lambx != 0:
        BETAx_denom = Integer(lcm([ci.denom() for ci in betax.list()]))
        BETAx = [Integer(ci) for ci in (BETAx_denom*betax).list()]
        LAMBx_denom = Integer(lcm([ci.denom() for ci in lambx.list()]))
        LAMBx = [Integer(ci) for ci in (LAMBx_denom*lambx).list()]
    else:
        BETAx_denom=0; BETAx=0; LAMBx_denom=0; LAMBx=0
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D, k

class FotiadisMartindale(BrezingWeng):
    """A FotiadisMartindale pairing-friendly curve of embedding degree k from eprint 2019/555
    """
    # re-use the init function of class BrezingWeng
    def __init__(self, u, code, a=None, b=None, cofactor_r=1, verbose_init=False):
        px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, D, k = coeffs_params(code)
        super().__init__(k, D, u, px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax, betax_denom, lambx, lambx_denom, a, b, cofactor_r, verbose_init)

    def _repr_(self):
        return "FotiadisMartindale_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (pairing-friendly curve k="+str(self._k)+") with seed "+str(self._u)+"\n"+super(BrezingWeng, self)._repr_()
