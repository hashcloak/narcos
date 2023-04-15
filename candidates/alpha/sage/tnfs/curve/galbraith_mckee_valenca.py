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
from sage.rings.integer_ring import Z, ZZ
from sage.rings.rational_field import Q, QQ
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.number_field.number_field import NumberField
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs
import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import compute_beta_lambda
from tnfs.curve.pairing_friendly_curve import BrezingWeng
from tnfs.curve.pairing_friendly_curve import get_curve_generator_order_r

"""
eprint 2004/365
Galbraith McKee Valenca
Ordinary abelian varieties having small embedding degree
Published in Finite Fields and Their Applications, 13 (2007) 800-814
"""

QQx = QQ['x']; x = QQx.gen()
#ZZx = ZZ['x']; (x,) = ZZx._first_ngens(1)
params_k = {
    3: {# k=3
        2: [#h=2
            {'k': 3, 'h': 2, 'px': 8*x**2 + 2*x + 1, 'tx': -2*x},
            {'k': 3, 'h': 2, 'px': 56*x**2 + 6*x - 1, 'tx': -14*x - 2},
            {'k': 3, 'h': 2, 'px': 56*x**2 + 22*x + 1, 'tx': -14*x - 4}],
        3: [#h=3
            {'k': 3, 'h': 3, 'px': 12*x**2 + 8*x + 3, 'tx': 2*x + 1}],
        4: [#h=4
            {'k': 3, 'h': 4, 'px': 16*x**2 + 6*x + 3, 'tx': -2*x},
            {'k': 3, 'h': 4, 'px': 48*x**2 + 30*x + 5, 'tx': 6*x + 2},
            {'k': 3, 'h': 4, 'px': 112*x**2 + 26*x + 1, 'tx': -14*x - 2},
            {'k': 3, 'h': 4, 'px': 112*x**2 + 58*x + 7, 'tx': -14*x - 4}],
        5: [#h=5
            {'k': 3, 'h': 5, 'px': 20*x**2 + 12*x + 5, 'tx': 2*x + 1},
            {'k': 3, 'h': 5, 'px': 140*x**2 + 64*x + 7, 'tx': 14*x + 3},
            {'k': 3, 'h': 5, 'px': 140*x**2 + 104*x + 19, 'tx': 14*x + 5},
            {'k': 3, 'h': 5, 'px': 260*x**2 + 44*x + 1, 'tx': -26*x - 3},
            {'k': 3, 'h': 5, 'px': 260*x**2 + 164*x + 25, 'tx': -26*x - 9},
            {'k': 3, 'h': 5, 'px': 380*x**2 + 112*x + 7, 'tx': -38*x - 7},
            {'k': 3, 'h': 5, 'px': 380*x**2 + 192*x + 23, 'tx': -38*x - 11}]
    },
    4: {# k=4
        3: [# h=3
            {'k': 4, 'h': 3, 'px': 12*x**2 + 2*x + 3, 'tx': 2*x + 1},
            {'k': 4, 'h': 3, 'px': 60*x**2 + 14*x + 1, 'tx': -10*x - 1},
            {'k': 4, 'h': 3, 'px': 60*x**2 + 34*x + 5, 'tx': 10*x + 3}],
        4: [# corrected from h=2
            {'k': 4, 'h': 4, 'px': 8*x**2 + 6*x + 3, 'tx': -2*x}],
        5: [# h=5
            {'k': 4, 'h': 5, 'px': 20*x**2 + 2*x + 5, 'tx': 2*x + 1},
            {'k': 4, 'h': 5, 'px': 260*x**2 + 134*x + 17, 'tx': -26*x - 7},
            {'k': 4, 'h': 5, 'px': 260*x**2 + 186*x + 33, 'tx': 26*x + 9},
            {'k': 4, 'h': 5, 'px': 340*x**2 + 46*x + 1, 'tx': -34*x - 3},
            {'k': 4, 'h': 5, 'px': 340*x**2 + 114*x + 9, 'tx': 34*x + 5}],
        6: [# corrected from h=3
            {'k': 4, 'h': 6, 'px': 12*x**2 + 10*x + 5, 'tx': -2*x},
            {'k': 4, 'h': 6, 'px': 60*x**2 + 26*x + 3, 'tx': -10*x - 2},
            {'k': 4, 'h': 6, 'px': 60*x**2 + 46*x + 9, 'tx': 10*x + 4}],
        8: [# corrected from h=4
            {'k': 4, 'h': 8, 'px': 16*x**2 + 14*x + 7, 'tx': -2*x},
            {'k': 4, 'h': 8, 'px': 80*x**2 + 38*x + 5, 'tx': -10*x - 2},
            {'k': 4, 'h': 8, 'px': 80*x**2 + 58*x + 11, 'tx': 10*x + 4},
            {'k': 4, 'h': 8, 'px': 208*x**2 + 54*x + 3, 'tx': -26*x - 4},
            {'k': 4, 'h': 8, 'px': 208*x**2 + 106*x + 13, 'tx': 26*x + 6}],
        10:[# corrected from h=5
            {'k': 4, 'h': 10, 'px': 20*x**2 + 18*x + 9, 'tx': -2*x},
            {'k': 4, 'h': 10, 'px': 260*x**2 + 74*x + 5, 'tx': -26*x - 4},
            {'k': 4, 'h': 10, 'px': 260*x**2 + 126*x + 15, 'tx': 26*x + 6},
            {'k': 4, 'h': 10, 'px': 340*x**2 + 226*x + 37, 'tx': -34*x - 12},
            {'k': 4, 'h': 10, 'px': 340*x**2 + 294*x + 63, 'tx': 34*x + 14}],
    },
    6: {# k=6
        2: [# h=2
            {'k': 6, 'h': 2, 'px': 8*x**2 + 6*x + 3, 'tx': 2*x + 2},
            {'k': 6, 'h': 2, 'px': 24*x**2 + 6*x + 1, 'tx': -6*x}],
        3: [# h=3
            {'k': 6, 'h': 3, 'px': 12*x**2 + 4*x + 3, 'tx': -2*x + 1},
            {'k': 6, 'h': 3, 'px': 84*x**2 + 16*x + 1, 'tx': -14*x - 1},
            {'k': 6, 'h': 3, 'px': 84*x**2 + 128*x + 49, 'tx': 14*x + 11}],
        4: [# h=4
            {'k': 6, 'h': 4, 'px': 16*x**2 + 10*x + 5, 'tx': 2*x + 2},
            {'k': 6, 'h': 4, 'px': 112*x**2 + 54*x + 7, 'tx': 14*x + 4},
            {'k': 6, 'h': 4, 'px': 112*x**2 + 86*x + 17, 'tx': 14*x + 6},
            {'k': 6, 'h': 4, 'px': 208*x**2 + 30*x + 1, 'tx': -26*x - 2}, # the curve in SNARK with seed x=0x8eed757d90615e40000000
            {'k': 6, 'h': 4, 'px': 208*x**2 + 126*x + 19, 'tx': -26*x - 8}],
        5: [# h=5
            {'k': 6, 'h': 5, 'px': 20*x**2 + 8*x + 5, 'tx': -2*x + 1},
            {'k': 6, 'h': 5, 'px': 60*x**2 + 36*x + 7, 'tx': 6*x + 3},
            {'k': 6, 'h': 5, 'px': 140*x**2 + 36*x + 3, 'tx': -14*x - 1},
            {'k': 6, 'h': 5, 'px': 140*x**2 + 76*x + 11, 'tx': -14*x - 3},
            {'k': 6, 'h': 5, 'px': 260*x**2 + 96*x + 9, 'tx': 26*x + 5},
            {'k': 6, 'h': 5, 'px': 260*x**2 + 216*x + 45, 'tx': 26*x + 11},
            {'k': 6, 'h': 5, 'px': 380*x**2 + 188*x + 23, 'tx': 38*x + 9},
            {'k': 6, 'h': 5, 'px': 380*x**2 + 268*x + 47, 'tx': 38*x + 13}]
    }
}

def allowed_k():
    return list(params_k)

def is_allowed_k(k: int):
    return k in allowed_k()

def allowed_h(k: int):
    if not is_allowed_k(k):
        raise ValueError("Error k={} does not exist for Galbraith-McKee-Valenca curves, only k in {} is available".format(k, allowed_k()))
    return list(params_k[k])

allowed_k = list(params_k)
allowed_h = {k: list(params_k[k]) for k in list(params_k)}
allowed_index = {k: {h: len(params_k[k][h]) for h in params_k[k]} for k in params_k}

def polynomial_params(k: int, h: int, index: int=0):
    """
    """
    if not k in params_k:
        raise ValueError("Error k={} does not exist for Galbraith-McKee-Valenca curves, only k in {} is available".format(k, list(params_k)))
    if not h in params_k[k]:
        raise ValueError("Error h={} does not exist for Galbraith-McKee-Valenca curves with k={}, only h in {} is available".format(h, k, list(params_k[k])))
    if index < 0 or index >= len(params_k[k][h]):
        raise ValueError("Error item at index {} does not exist for Galbraith-McKee-Valenca curves, only from 0 to {} for k={}, h={}".format(index, len(params_k[k][h])-1))
    item = params_k[k][h][index]
    px = item['px']
    tx = item['tx']
    rx = (px+1-tx)/h
    # for homogeneity with other functions, do not return k not h but return c = h as a constant polynomial
    cx = QQx(item['h'])
    return px, rx, tx, cx

class GMV(EllipticCurve_finite_field):
    """
    A Galbraith-McKee-Valenca curve of embedding degree 3,4 or 6
    Ordinary abelian varieties having small embedding degree
    Published in Finite Fields and Their Applications, 13 (2007) 800-814
    eprint 2004/365

    In fact there is a close relation between these curves and
    the curves at https://eprint.iacr.org/2004/058
    Generating more MNT elliptic curves, Scott and Barreto

    """
    def __init__(self, k, h, idx_family, u, a, b, D=None, c=None):
        """
        INPUT:
        - `k`: embedding degree, this is 3, 4 or 6
        - `h`: cofactor of r such that the curve order is p+1-t = h*r, with small h
        - `u`: seed such that p=P_GMV(k,h,idx_family,u)
        - `a`: curve coefficient in y^2 = x^3 + a*x + b
        - `b`: curve coefficient in y^2 = x^3 + a*x + b
        - `D`: curve discriminant in t^2 -4*p = -D*y^2
        - `c`: additional cofactor of r: p+1-t = h*c*(r/c)

        a, b s.t. y^2 = x^3 + a*x + b has exactly order r=R_GMV(k,h,idx_family,u)
        (there are two possible choices for b: b, 1/b (a twist) but only one has order r)
        """
        if not k in [3,4,6]:
            raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))
        self._k = k # embedding degree
        self._h = h # co-factor of the polynomial rx
        self._idx = idx_family # the item number of the family
        self._a = Integer(a)
        self._b = Integer(b)
        if D != None:
            self._D = Integer(D)
        self._px, self._px_denom, self._rx, self._rx_denom, self._tx, self._tx_denom, self._cx, self._cx_denom = coeff_params(k)

        self._u = Integer(u)
        self._p = (Integer(self._px[2])*self._u + Integer(self._px[1]))*self._u + Integer(self._px[0])
        self._pbits = self._p.nbits()
        self._r = (Integer(self._rx[2])*self._u + Integer(self._rx[1]))*self._u + Integer(self._rx[0])
        self._tr = Integer(self._tx[1])*self._u + Integer(self._tx[0])
        self._c = Integer(h)
        if c != None and c != 1 and (self._r % Integer(c)) == 0:
            self._c = self._c * Integer(c)
            self._r = self._r // self._c
        if not self._c * self._r == self._p + 1 - self._tr:
            raise ValueError("Error: r*c != p+1-tr\nr={}\nc={}\np+1-tr={}\n".format(self._r,self._c,self._p+1-self._tr))
        self._y = Integer(((4*self._p - self._tr**2)//self._D).sqrt())

        try:
            self._Fp = FiniteField(self._p)
        except ValueError as err:
            print("ValueError creating Fp: {}".format(err))
            print("p= {}".format(self._p))
            raise
        except:
            print("Error creating Fp")
            raise
        if not self._r.is_prime():
            raise ValueError("Error r is not prime")

        self._Fpz = PolynomialRing(self._Fp, names=('z',))
        (self._z,) = self._Fpz._first_ngens(1)

        if a != None:
            try:
                a = Integer(a)
            except:
                raise
            self._a = a
            self._ap = self._Fp(a)
        if b != None:
            try:
                b = Integer(b)
            except:
                raise
            self._b = b
            self._bp = self._Fp(b)

        try:
            # this init function of super inherits from class EllipticCurve_generic defined in ell_generic.py
            # __init__ method inherited from ell_generic
            EllipticCurve_finite_field.__init__(self, self._Fp, [0,0,0,self._ap,self._bp])
        except ValueError as err:
            print("ValueError at EllipticCurve_finite_field.__init__: {}".format(err))
            raise
        except:
            print("An error occupred when initialising the elliptic curve")
            raise
        #print("check order")
        #self.order_checked = super(GMV,self).order()
        #if self.order_checked != (self._p+1-self._tr):
        #    raise ValueError("Wrong curve order: this one is a twist")
        #print("check order (randomized)")
        self.curve_order = self._p + Integer(1) - self._tr
        self.twist_order = self._p + Integer(1) + self._tr
        for i in range(10):
            P = self.random_element()
            if self.curve_order*P != self(0):
                if self.twist_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a twist: (p+1+tr)*P = 0\ntr={}\nr={}\np+1+tr={}\n".format(self._tr,self.curve_order,self.twist_order))
                else:
                    self.order_checked = super(GMV,self).order()
                    raise ValueError("Wrong curve order:\np+1-tr        = {}\np+1+tr        = {}\nchecked order = {}\np             = {}\nchecked trace = {}\nexpected trace= {}".format(self.curve_order,self.twist_order,self.order_checked,self._p, self._p+Integer(1)-self.order_checked, self._tr))

        # computes a generator
        self._G = get_curve_generator_order_r(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]

    def _repr_(self):
        return "GMV-k{}-h{}-i{} p{} (Galbraith-McKee-Valenca curve k={} h={} id={}) curve with seed {}\n".format(self._k,self._h,self._idx,self._pbits,self._k,self._h,self._idx,self._u)+super(GMV,self)._repr_()

    def u(self):
        return self._u
    def p(self):
        return self._p
    def r(self):
        return self._r
    def c(self):
        return self._c
    def h(self):
        return self._h
    def idx(self):
        return self._idx
    def tr(self):
        return self._tr
    def y(self):
        return self._y
    def D(self):
        return self._D
    def a(self):
        return self._a # 0
    def ap(self):
        return self._ap # 0
    def b(self):
        return self._b # Integer
    def bp(self):
        return self._bp # in Fp (finite field element)

    def px(self):
        return self._px
    def rx(self):
        return self._rx
    def tx(self):
        return self._tx
    def cx(self):
        return self._cx

    def k(self):
        return self._k
    def Fp(self):
        return self._Fp
    def Fpz(self):
        return self._Fpz, self._z
    def G(self):
        return self._G

    def print_parameters(self):
        tnfs.curve.pairing_friendly_curve.print_parameters(self)

    def print_parameters_for_RELIC(self):
        tnfs.curve.pairing_friendly_curve.print_parameters_for_RELIC(self)
