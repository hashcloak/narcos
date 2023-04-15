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
    px = (x-1)**2*(x**8 - x**4 + 1)/3 + x #= (x^10-2*x^9+x^8-x^6+2*x^5-x^4+x^2+x+1)/3
    rx = x**8 - x**4 + 1
    cx = (1/3) * (x - 1)**2 # cofactor such that p+1-tr == c*r
    tx = x + 1
    yx = (x-1)*(2*x**4 - 1)/3 # = (2*x^5 - 2*x^4 - x + 1)/3
    # Y such that T^2 - 4*P = -3*Y^2
    betax = x**9 - 3*x**8 + 4*x**7 - 4*x**6 + 3*x**5 - 2*x**3 + 2*x**2 - x + 1
    lambx = x**4 -1 # or -x**4
    D = 3
    return px, rx, tx, cx, yx, betax, lambx, D

def coeff_params():
    Px = [1, 1, 1, 0, -1, 2, -1, 0, 1, -2, 1]
    Px_denom = 3
    Rx = [1, 0, 0, 0, -1, 0, 0, 0, 1]
    Rx_denom = 1
    Cx = [1, -2, 1]
    Cx_denom = 3
    Tx = [1, 1]
    Tx_denom = 1
    Yx = [1, -1, 0, 0, -2, 2]
    Yx_denom = 3
    BETAx = [1, -1, 2, -2, 0, 3, -4, 4, -3, 1]
    BETAx_denom = 1
    LAMBx = [-1, 0, 0, 0, 1]
    LAMBx_denom = 1
    D=3
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D

class BLS24(EllipticCurve_finite_field):
    """
    A Barreto-Lynn-Scott curve of embedding degree 24
    """
    def __init__(self, u, b=None):
        """
        u is the seed such that p=P_BLS24(u), and b is the curve coefficient, s.t. y^2 = x^3 + b has a subgroup of order r=R_BLS24(u)
        (there are two possible choices for b: b, 1/b (a twist) but only one has order r)
        """
        self._k = 24 # embedding degree
        self._a = 0 # first curve parameter is 0 because j=0
        self._D = 3
        # polynomials (will be needed later for polynomial selection)
        # P_BLS24 = (u-1)^2*(u^8 - u^4 + 1)/3 + u #= (s^10-2*s^9+s^8-s^6+2*s^5-s^4+s^2+s+1)/3
        # R_BLS24 = u^8 - u^4 + 1
        # C_BLS24 = (1/3) * (u - 1)^2 # cofactor such that p+1-tr == c*r
        # T_BLS24 = u + 1
        # Y_BLS24 = (u-1)*(2*u^4 - 1)/3 # = (2*s^5 - 2*s^4 - s + 1)/3
        # Y such that T^2 - 4*P = -3*Y^2
        # BETA_BLS24 = s^9 - 3*s^8 + 4*s^7 - 4*s^6 + 3*s^5 - 2*s^3 + 2*s^2 - s + 1
        # LAMB_BLS24 = s^4 -1 # or -s^4
        # do not divides by 3 so that it has integer coeffs.
        # pb type with 1/3 (as an int, it is 0) -> uses Integer(1)/Integer(3)
        self.polynomial_p = [1, 1, 1, 0, -1, 2, -1, 0, 1, -2, 1]
        self.polynomial_p_denom = 3
        self.polynomial_r = [1, 0, 0, 0, -1, 0, 0, 0, 1]
        self.polynomial_r_denom = 1
        self.polynomial_c = [1, -2, 1]
        self.polynomial_c_denom = 3
        self.polynomial_tr = [1, 1]
        self.polynomial_tr_denom = 1
        self.polynomial_y = [1, -1, 0, 0, -2, 2]
        self.polynomial_y_denom = 3
        self.polynomial_beta = [1, -1, 2, -2, 0, 3, -4, 4, -3, 1]
        self.polynomial_beta_denom = 1
        self.polynomial_lamb = [-1, 0, 0, 0, 1]
        self.polynomial_lamb_denom = 1

        self._u = Integer(u)
        self._p = (self._u-1)**2*(self._u**8 - self._u**4 + 1)//3 + self._u
        self._pbits = self._p.nbits()
        self._tr = self._u + 1
        self._r = self._u**8 - self._u**4 + 1
        self._c = (self._u-1)**2//3 # cofactor
        if not self._c * self._r == self._p + 1 - self._tr:
            raise ValueError("Error: r*c != p+1-tr\nr={}\nc={}\np+1-tr={}\n".format(self._r,self._c,self._p+1-self._tr))
        self._y =  ((self._u-1) * (2*self._u**4 - 1))//3
        # GLV parameter for fast scalar multiplication thanks to automorphism
        # psi: (x,y) -> (beta*x,y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+x+1 mod p, beta = (-1 + sqrt(Fp(-3)))/2
        # and lambda is a root of x^2+x+1 mod r, lamb = (-1 + sqrt(Fr(-3)))/2
        # there are two choices: beta, -beta-1, and the same for lamb: lamb, -lamb-1
        # arbitrarily the positive ones are chosen
        self._beta = sum([Integer(self.polynomial_beta[i])*self._u**i for i in range(len(self.polynomial_beta))])
        self._lamb = self._u**4 - 1
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
        self.order_checked = super(BLS24,self).order()
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
                self._beta = -self._beta-1
                print("beta -> -beta-1")
            elif self([self._Gx*self._beta, self._Gy]) == (-self._lamb-1)*self._G:
                self._lamb = -self._lamb-1
                print("lamb -> -lamb-1")
            elif self([self._Gx*(-self._beta-1), self._Gy]) == (-self._lamb-1)*self._G:
                self._beta = -self._beta-1
                self._lamb = -self._lamb-1
                print("lamb -> -lamb-1")
                print("beta -> -beta-1")
            else:
                raise ValueError("Error while adjusting beta, lamb: compatibility not found")

    def _repr_(self):
        return "BLS24 p"+str(self._pbits)+" (Barreto-Lynn-Scott k=24) curve with seed "+str(self._u)+"\n"+super(BLS24,self)._repr_()

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

