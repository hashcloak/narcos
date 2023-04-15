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
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0

def polynomial_params():
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px = (x**8 + 5*x**7 + 7*x**6 + 37*x**5 + 188*x**4 + 259*x**3 + 343*x**2 + 1763*x + 2401)/21
    rx = (x**6 + 37*x**3 + 343)/343 # 343 = 7^3
    cx = (x**2 + 5*x + 7)*49/3
    tx = (x**4 + 16*x + 7)/7
    yx = (5*x**4 + 14*x**3 + 94*x + 259)/21 # Y such that T^2 - 4*P = -3*Y^2
    betax = (x**7 + 3*x**6 + 4*x**5 + 44*x**4 + 118*x**3 + 71*x**2 + 483*x + 1118)/24
    lambx = x**3 + 18
    D=3
    return px, rx, tx, cx, yx, betax, lambx, D

def coeff_params():
    Px = [2401, 1763, 343, 259, 188, 37, 7, 5, 1]
    Px_denom = 21
    Rx = [343, 0, 0, 37, 0, 0, 1]
    Rx_denom = 343 # or 21?
    Cx = [343, 245, 49] # 49*[7, 5, 1]
    Cx_denom = 3 # 21 # or 1, and denom with r?
    Tx= [7, 16, 0, 0, 1]
    Tx_denom = 7
    Yx = [259, 94, 0, 14, 5]
    Yx_denom = 21
    BETAx = [1118, 483, 71, 118, 44, 4, 3, 1]
    BETAx_denom = 24
    LAMBx = [18, 0, 0, 1]
    LAMBx_denom = 1
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom

class KSS18(EllipticCurve_finite_field):
    """
    A Kachisa-Schaefer-Scott curve of embedding degree 18
    """
    def __init__(self, u, b=None):
        """
        u is the seed such that p=P_KSS18(u), and b is the curve coefficient, s.t. y^2 = x^3 + b has a subgroup of order r=R_KSS18(u)
        One needs u=7,14 mod 21 to ensure that P(u) and R(u) will be integers
        (there are two possible choices for b: b, 1/b (a twist) but only one has order r)
        """
        self._k = 18 # embedding degree
        self._a = 0 # first curve parameter is 0 because j=0
        self._D = 3
        # polynomials (will be needed later for polynomial selection)
        # P_KSS18 = 1/21*(s^8 + 5*s^7 + 7*s^6 + 37*s^5 + 188*s^4 + 259*s^3 + 343*s^2 + 1763*s + 2401)
        # R_KSS18 = (s^6 + 37*s^3 + 343)/343 # 343 = 7^3
        # C_KSS18 = (s^2 + 5*s + 7)*49/3
        # T_KSS18 = 1/7 * (s^4+16*s+7)
        # Y_KSS18 = (1/21) * (5*s^4 + 14*s^3 + 94*s + 259) # Y such that T^2 - 4*P = -3*Y^2
        # BETA_KSS18 = (s^7 + 3*s^6 + 4*s^5 + 44*s^4 + 118*s^3 + 71*s^2 + 483*s + 1118)/24
        # LAMB_KSS18 = s^3 + 18
        # do not divides by 21 so that it has integer coeffs.
        # pb type with 1/21 (as an int, it is 0) -> uses Integer(1)/Integer(21)
        self.polynomial_p = [2401, 1763, 343, 259, 188, 37, 7, 5, 1]
        self.polynomial_p_denom = 21
        self.polynomial_r = [343, 0, 0, 37, 0, 0, 1]
        self.polynomial_r_denom = 343 # or 21?
        self.polynomial_c = [343, 245, 49] # 49*[7, 5, 1]
        self.polynomial_c_denom = 3 # 21 # or 1, and denom with r?
        self.polynomial_tr= [7, 16, 0, 0, 1]
        self.polynomial_tr_denom= 7
        self.polynomial_y = [259, 94, 0, 14, 5]
        self.polynomial_y_denom = 21
        self.polynomial_beta = [1118, 483, 71, 118, 44, 4, 3, 1]
        self.polynomial_beta_denom = 24
        self.polynomial_lamb = [18, 0, 0, 1]
        self.polynomial_lamb_denom = 1

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
        # psi: (x,y) -> (beta*x,y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+x+1 mod p, beta = (-1 + sqrt(Fp(-3)))/2
        # and lambda is a root of x^2+x+1 mod r, lamb = (-1 + sqrt(Fr(-3)))/2
        # there are two choices: beta, -beta-1, and the same for lamb: lamb, -lamb-1
        # arbitrarily the positive ones are chosen
        if (u % 24) == 2 or (u % 24) == 10: # beta is actually an integer
            self._beta = sum([Integer(self.polynomial_beta[i])*self._u**i for i in range(len(self.polynomial_beta))])//Integer(self.polynomial_beta_denom)
        else: # computes modulo p
            self._beta = sum([Integer(self.polynomial_beta[i])*self._u**i for i in range(len(self.polynomial_beta))])/Integer(self.polynomial_beta_denom)
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

        if not ((u % 24) == 2 or (u % 24) == 10):
            # needs to be reduced mod p, otherwise has a denominator
            self._beta = Integer(self._Fp(self._beta))
        if ((self._beta**2 + self._beta + 1) % self._p) != 0:
            raise ValueError("Error beta^2 + beta + 1 != 0 mod p")
        if ((self._lamb**2 + self._lamb + 1) % self._r) != 0:
            raise ValueError("Error lamb^2 + lamb + 1 != 0 mod r")
        
        self._Fpz = PolynomialRing(self._Fp, names=('z',))
        (self._z,) = self._Fpz._first_ngens(1)

        if b != None:
            try:
                b = Integer(b)
            except:
                raise
            self._b = b
            self._bp = self._Fp(b)
        else:
            # check that beta = 2*U/(-3*V-U) before, where U=t/2, V = y/2 and 2V = 2 mod 3
            self._b, self._bp = get_curve_parameter_b_j0(self._tr, self._y, self._p, self._Fp)
        self._ap = self._Fp(0) # first curve parameter is 0 because j=0
        # Now self._b is such that E: y^2 = x^3 + b has order r
        try:
            # this init function of super inherits from class EllipticCurve_generic defined in ell_generic.py
            # __init__ method inherited from ell_generic
            EllipticCurve_finite_field.__init__(self, self._Fp, [0,0,0,0,self._bp])
        except ValueError as err:
            print("ValueError at EllipticCurve_finite_field.__init__: {}".format(err))
            raise
        except:
            print("An error occupred when initialising the elliptic curve")
            raise
        self.order_checked = super(KSS18,self).order()
        if self.order_checked != (self._p+1-self._tr):
            print("Error, wrong order")
            raise ValueError("Wrong curve order: this one is a twist")

        # computes a generator
        self._G = get_curve_generator_order_r_j0(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]
        
        # adjust beta and lamb according to the curve
        # do we have (beta*x,y) = lamb*(x,y)?
        if self([self._Gx*self._beta, self._Gy]) != self._lamb*self._G:
            print("adjusting beta, lambda")
            if self([self._Gx*(-self._beta-1), self._Gy]) == self._lamb*self._G:
                if ((u % 24) == 2 or (u % 24) == 10): # no denominator
                    self._beta = -self._beta-1
                else:
                    self._beta = self._p-self._beta-1
                print("beta -> -beta-1")
            elif self([self._Gx*self._beta, self._Gy]) == (-self._lamb-1)*self._G:
                self._lamb = -self._lamb-1
                print("lamb -> -lamb-1")
            elif self([self._Gx*(-self._beta-1), self._Gy]) == (-self._lamb-1)*self._G:
                if ((u % 24) == 2 or (u % 24) == 10): # no denominator
                    self._beta = -self._beta-1
                else:
                    self._beta = self._p-self._beta-1
                self._lamb = -self._lamb-1
                print("lamb -> -lamb-1")
                print("beta -> -beta-1")
            else:
                raise ValueError("Error while adjusting beta, lamb: compatibility not found")

    def _repr_(self):
        return "KSS18 p"+str(self._pbits)+" (Kachisa-Schaefer-Scott k=18) curve with seed "+str(self._u)+"\n"+super(KSS18,self)._repr_()

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
        return self._a # 0
    def ap(self):
        return self._ap # 0
    def b(self):
        return self._b # Integer
    def bp(self):
        return self._bp # in Fp (finite field element)
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

