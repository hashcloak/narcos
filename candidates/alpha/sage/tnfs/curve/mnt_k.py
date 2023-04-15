# MNT curves with k=4
# the j-invariant is not special (j is not 0 or 1728)
import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.functions.other import sqrt
from sage.arith.misc import XGCD, xgcd
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_generator_order_r

def polynomial_params(k):
    if not k in [3,4,6]:
            raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))        
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    cx = QQx(1)
    if k==6:
        px = 4*x**2+1
        rx = 4*x**2 - 2*x + 1
        tx = 2*x + 1
        #D*Y^2 = 12*x^2 - 4*x + 3
        # no beta, lambda
    elif k==4:
        px = x**2 + x + 1
        rx = x**2 + 1
        tx = x + 1
        # t^2 - 4*p = (x+1)^2 - 4*x^2 -4*x -4 = -3*x^2 -2*x -3
        # an alternative representation, euivalent, is
        # px = x**2 + x + 1
        # rx = x**2 + 2*x + 2
        # tx = -x
        # t^2 - 4*p = x^2 - 4*x^2 - 4*x - 4 = -3*x^2 -4*x -4
    elif k==3:
        px = 12*x**2 - 1
        rx = 12*x**2 - 6*x + 1
        tx = 6*x - 1
    return px, rx, tx, cx

def coeff_params(k):
    if not k in [3,4,6]:
            raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))        
    Cx = [1]
    Cx_denom = 1
    Px_denom = 1
    Rx_denom = 1
    Tx_denom = 1
    if k == 6:
        Px = [1, 0, 4]
        Rx = [1, -2, 4]
        Tx = [1, 2]
    elif k == 4:
        Px = [1, 1, 1]
        Rx = [1, 0, 1]
        Tx = [1, 1]
        # alternatively one could use
        #Px = [1, 1, 1]
        #Rx = [2, 2, 1]
        #Tx = [0, -1]
    elif k == 3:
        Px = [-1, 0, 12]
        Rx = [1, -6, 12]
        Tx = [-1, 6]
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom

class MNT_k(EllipticCurve_finite_field):
    """
    A Miyaji-Nakabayashi-Takano curve of embedding degree 3,4 or 6
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.20.8113&rep=rep1&type=pdf
    """
    def __init__(self, k, u, a, b,D=None, c=None):
        """
        k is the embedding degree, this is 3, 4 or 6
        u is the seed such that p=P_MNT_k(u), and a,b are the curve coefficients, 
        s.t. y^2 = x^3 + a*x + b has exactly order r=R_MNT_k(u)
        (there are two possible choices for b: b, 1/b (a twist) but only one has order r)
        """
        if not k in [3,4,6]:
            raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))
        self._k = k # embedding degree
        self._a = Integer(a)
        self._b = Integer(b)
        if D != None:
            self._D = Integer(D)
        self._px, self._px_denom, self._rx, self._rx_denom, self._tx, self._tx_denom, self._cx, self._cx_denom = coeff_params(k)

        self._u = Integer(u)
        #self._p = sum([Integer(self._px[i])*self._u**i for i in range(len(self._px))])
        self._p = Integer(self._px[2])*self._u**2 + Integer(self._px[1])*self._u + Integer(self._px[0])
        self._pbits = self._p.nbits()
        #self._r = sum([Integer(self._rx[i])*self._u**i for i in range(len(self._rx))])
        self._r = Integer(self._rx[2])*self._u**2 + Integer(self._rx[1])*self._u + Integer(self._rx[0])
        #self._tr= sum([Integer(self._tx[i])*self._u**i for i in range(len(self._tx))])
        self._tr = Integer(self._tx[1])*self._u + Integer(self._tx[0])
        #self._c = sum([Integer(self._cx[i])*self._u**i for i in range(len(self._cx))])
        self._c = Integer(1)
        if c != None and c != 1 and (self._r % Integer(c)) == 0:
            self._c = Integer(c)
            self._r = self._r // self._c
        else:
            self._c = Integer(1)
        if not self._c * self._r == self._p + 1 - self._tr:
            raise ValueError("Error: r*c != p+1-tr\nr={}\nc={}\np+1-tr={}\n".format(self._r,self._c,self._p+1-self._tr))
        #print("computing y")
        self._y = Integer(((4*self._p - self._tr**2)//self._D).sqrt())
        #print("y done")

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
        #self.order_checked = super(MNT_k,self).order()
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
                    self.order_checked = super(MNT_k,self).order()
                    raise ValueError("Wrong curve order:\np+1-tr        = {}\np+1+tr        = {}\nchecked order = {}\np             = {}\nchecked trace = {}\nexpected trace= {}".format(self.curve_order,self.twist_order,self.order_checked,self._p, self._p+Integer(1)-self.order_checked, self._tr))
        #print("ok")
        
        # computes a generator
        self._G = get_curve_generator_order_r(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]

    def _repr_(self):
        return "MNT{} p{} (Miyaji-Nakabayashi-Takano curve k={}) curve with seed {}\n".format(self._k,self._pbits,self._k,self._u)+super(MNT_k,self)._repr_()

    def u(self):
        return self._u
    def p(self):
        return self._p
    def r(self):
        return self._r
    def c(self):
        return self._c
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

