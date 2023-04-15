import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.arith.misc import XGCD, xgcd
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_a_j1728, get_curve_generator_order_r_j1728

def polynomial_params():
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px = (x**10 + 2*x**9 + 5*x**8 + 48*x**6 + 152*x**5 + 240*x**4 + 625*x**2 + 2398*x + 3125)/980
    rx = (x**8 + 48*x**4 + 625)/61250 # 625 = 5^4, 61250 = 2*5^4*7^2
    tx = (2*x**5 + 41*x + 35)/35
    cx = 125 * (x**2 + 2*x + 5)/2 # C such that P+1-T = C*R
    yx = (x**5 + 5*x**4 + 38*x + 120)/70 # Y such that T^2 - 4*P = -4*Y^2
    # assert (px+1-tx) == cx*rx
    # assert (tx^2 + yx^2)/4 == px
    # sqrt(-1) mod P
    betax = (x**9-11*x**8-16*x**7-120*x**6-32*x**5-498*x**4-748*x**3-2740*x**2-3115*x-5651)/4018
    lambx = (x**4 + 24)/7 # sqrt(-1) mod R
    D=1
    return px, rx, tx, cx, yx, betax, lambx, D

def coeff_params():
    Px = [3125, 2398, 625, 0, 240, 152, 48, 0, 5, 2, 1]
    Px_denom = 980
    Px_Special = [1492, 628, 1378, -632, 508, -416, 230, -88, 32, -8, 1] # this is Px(x-1)
    Rx = [625, 0, 0, 0, 48, 0, 0, 0, 1]
    Rx_denom = 61250
    # when u = 25, 45 mod 70, then P(u), T(u), (P+1-T)(u) are all integers, and 61250 | R(u)
    Cx = [625, 250, 125]
    Cx_denom = 2
    Tx = [35, 41, 0, 0, 0, 2]
    Tx_denom = 35 # we need u = 25,45 mod 70 for T to be integer
    Yx = [120, 38, 0, 0, 5, 1]
    Yx_denom = 70
    BETAx = [-5651, -3115, -2740, -748, -498, -32, -120, -16, -11, 1]
    BETAx_denom = 4018 # 2 * 7^2 * 41
    LAMBx = [24, 0, 0, 0, 1]
    LAMBx_denom = 7
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom

class KSS16(EllipticCurve_finite_field):
    """
    A Kachisa-Schaefer-Scott curve of embedding degree 16
    """
    def __init__(self, u, a=None):
        """
        u is the seed such that p=P_KSS16(u), and a is the curve coefficient, s.t. y^2 = x^3 + a*x has a subgroup of order r=R_KSS16(u)
        (there are two possible choices for a: a, 1/a (a twist) but only one has order r)
        """
        self._k = 16 # embedding degree
        self._b = 0 # second curve parameter is 0 because j=1728
        self._D = 4 # 4, not 1, be careful with that
        # polynomials (will be needed later for polynomial selection)
        # P_KSS16 = (s^10 + 2*s^9 + 5*s^8 + 48*s^6 + 152*s^5 + 240*s^4 + 625*s^2 + 2398*s + 3125)/980
        # R_KSS16 = (s^8 + 48*s^4 + 625)/61250 # 625 = 5^4, 61250 = 2*5^4*7^2
        # T_KSS16 = (2*s^5 + 41*s + 35)/35
        # C_KSS16 = (125/2) * (s^2 + 2*s + 5) # C such that P+1-T = C*R
        # Y_KSS16 =  (1/70) * (s^5 + 5*s^4 + 38*s + 120) # Y such that T^2 - 4*P = -4*Y^2
        # assert (P_KSS16+1-T_KSS16) == C_KSS16*R_KSS16
        # assert (T_KSS16^2 + D_KSS16*Y_KSS16^2)/4 == P_KSS16
        # sqrt(-1) mod P
        # BETA_KSS16 = (s^9-11*s^8-16*s^7-120*s^6-32*s^5-498*s^4-748*s^3-2740*s^2-3115*s-5651)/4018
        # LAMB_KSS16 = (1/7)*(s^4 + 24) # sqrt(-1) mod R
        self.polynomial_p = [3125, 2398, 625, 0, 240, 152, 48, 0, 5, 2, 1]
        self.polynomial_p_denom = 980
        self.polynomial_r = [625, 0, 0, 0, 48, 0, 0, 0, 1]
        self.polynomial_r_denom = 61250
        # when u = 25, 45 mod 70, then P(u), T(u), (P+1-T)(u) are all integers, and 61250 | R(u)
        self.polynomial_c = [625, 250, 125]
        self.polynomial_c_denom = 2
        self.polynomial_tr= [35, 41, 0, 0, 0, 2]
        self.polynomial_tr_denom= 35 # we need u = 25,45 mod 70 for T to be integer
        self.polynomial_y = [120, 38, 0, 0, 5, 1]
        self.polynomial_y_denom = 70
        self.polynomial_beta = [-5651, -3115, -2740, -748, -498, -32, -120, -16, -11, 1]
        self.polynomial_beta_denom = 4018 # 2 * 7^2 * 41
        self.polynomial_lamb = [24, 0, 0, 0, 1]
        self.polynomial_lamb_denom = 7

        self._u = Integer(u)
        self._p = sum([Integer(self.polynomial_p[i])*self._u**i for i in range(len(self.polynomial_p))])//self.polynomial_p_denom
        self._pbits = self._p.nbits()
        self._r = sum([Integer(self.polynomial_r[i])*self._u**i for i in range(len(self.polynomial_r))])//self.polynomial_r_denom
        self._c = sum([Integer(self.polynomial_c[i])*self._u**i for i in range(len(self.polynomial_c))])//self.polynomial_c_denom
        self._tr= sum([Integer(self.polynomial_tr[i])*self._u**i for i in range(len(self.polynomial_tr))])//self.polynomial_tr_denom
        if not self._c * self._r == self._p + 1 - self._tr:
            raise ValueError("Error: r*c != p+1-tr\nr={}\nc={}\np+1-tr={}\n".format(self._r,self._c,self._p+1-self._tr))
        self._y = sum([Integer(self.polynomial_y[i])*self._u**i for i in range(len(self.polynomial_y))])//self.polynomial_y_denom
        # GLV parameter for fast scalar multiplication thanks to automorphism
        # psi: (x,y) -> (-x,i*y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+1 mod p, beta = sqrt(-1) mod p
        # and lambda is a root of x^2+1 mod r, lamb = sqrt(-1) mod r
        # there are two choices: beta, -beta, and the same for lamb: lamb, -lamb
        # arbitrarily the positive ones are chosen
        # unfortunately beta has always denominator 41, computes modulo p: / instead of //
        self._beta = sum([Integer(self.polynomial_beta[i])*self._u**i for i in range(len(self.polynomial_beta))])/Integer(self.polynomial_beta_denom)
        # LAMB is an integer when u=25,45 mod 70, which should be the case for P, T, R to be integers.
        self._lamb = sum([Integer(self.polynomial_lamb[i])*self._u**i for i in range(len(self.polynomial_lamb))])//self.polynomial_lamb_denom
        try:
            self._Fp = FiniteField(self._p)
        except ValueError as err:
            print("ValueError creating Fp: {}".format(err))
            raise
        except:
            print("Error creating Fp")
            raise
        if not self._r.is_prime():
            raise ValueError("Error r is not prime")

        # beta has denominator 41, reduces mod p:
        self._beta = Integer(self._Fp(self._beta))
        if ((self._beta**2 + 1) % self._p) != 0:
            raise ValueError("Error beta^2 + 1 != 0 mod p")
        if ((self._lamb**2 + 1) % self._r) != 0:
            raise ValueError("Error lamb^2 + 1 != 0 mod r")
        
        self._Fpz = PolynomialRing(self._Fp, names=('z',))
        (self._z,) = self._Fpz._first_ngens(1)

        self._bp = self._Fp(0) # second curve parameter is 0 because j=1728
        if a != None:
            try:
                a = Integer(a)
            except:
                raise
            self._a = a
            self._ap = self._Fp(a)
        else:
            # check that beta = U/V, where U=t/2, V = y
            self._a, self._ap = get_curve_parameter_a_j1728(self._tr, self._y, self._p, self._Fp)
        # Now self._a is such that E: y^2 = x^3 + a*x has order r
        try:
            # this init function of super inherits from class EllipticCurve_generic defined in ell_generic.py
            # __init__ method inherited from ell_generic
            EllipticCurve_finite_field.__init__(self, self._Fp, [0,0,0,self._ap,0])
        except ValueError as err:
            print("ValueError at EllipticCurve_finite_field.__init__: {}".format(err))
            raise
        except:
            print("An error occupred when initialising the elliptic curve")
            raise
        self.order_checked = super(KSS16,self).order()
        if self.order_checked != (self._p+1-self._tr):
            print("Error, wrong order")
            raise ValueError("Wrong curve order: this one is a twist")

        # computes a generator
        self._G = get_curve_generator_order_r_j1728(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]
        
        # adjust beta and lamb according to the curve
        # do we have (-x,beta*y) = lamb*(x,y)?
        if self([-self._Gx, self._Gy*self._beta]) != self._lamb*self._G:
            print("adjusting beta, lambda")
            if self([-self._Gx, self._Gy*(-self._beta)]) == self._lamb*self._G:
                self._beta = self._p-self._beta
                print("beta -> -beta")
            elif self([-self._Gx, self._Gy*self._beta]) == (-self._lamb)*self._G:
                self._lamb = -self._lamb
                print("lamb -> -lamb")
            elif self([-self._Gx, self._Gy*(-self._beta)]) == (-self._lamb)*self._G:
                self._beta = self._p-self._beta
                self._lamb = -self._lamb
                print("lamb -> -lamb")
                print("beta -> -beta")
            else:
                raise ValueError("Error while adjusting beta, lamb: compatibility not found")

    def _repr_(self):
        return "KSS16 p"+str(self._pbits)+" (Kachisa-Schaefer-Scott k=16) curve with seed "+str(self._u)+"\n"+super(KSS16,self)._repr_()

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
    def a(self):
        return self._a # Integer
    def ap(self):
        return self._ap # in Fp (finite field element)
    def b(self):
        return self._b # 0
    def bp(self):
        return self._bp # 0
    def beta(self):
        return self._beta
    def lamb(self):
        return self._lamb
    
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

