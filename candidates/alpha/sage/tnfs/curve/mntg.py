# Generalised MNT curves with k=6
# Modification of MNT by Scott-Barreto (DCC06, preprint https://eprint.iacr.org/2004/058)
import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.functions.other import sqrt
from sage.arith.misc import XGCD, xgcd, gcd
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.misc.functional import cyclotomic_polynomial

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_generator_order_r

def polynomial_params(k, c, d):
    """
    INPUT:
    - `k` embedding degree
    - `c`: co-factor of the curve, # E(Fp) = c*r where r is prime
    - `d`: co-factor of Phi_k(u): r*d = Phi_k(u), hence r = Phi_k(t-1)/d where t = u+1

    Returns the polynomials px, rx, tx, cx of the curve parameters
    """
    if not k in [3,4,6]:
        raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    Phi_k = QQx(cyclotomic_polynomial(k))
    rx = QQx(Phi_k/QQx(d))
    tx = QQx(x+1)
    cx = QQx(c)
    px = QQx(cx*rx + tx - 1)
    return px, rx, tx, cx

def congruence_constraints(k, d):
    """
    INPUT:
    - `k` embedding degree
    - `d`: co-factor of Phi_k(u): r*d = Phi_k(u), hence r = Phi_k(t-1)/d where t = u+1, and r in integer

    Returns the congruence constrains: m=d, u_m a list of i s.t. Phi_k(d*x+i)/d has integer coefficients
    """
    if not k in [3,4,6]:
        raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))
    if d == 1:
        u_m = [0]
        return 1, u_m
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    Phi_k = QQx(cyclotomic_polynomial(k))
    u_m = []
    for i in range(d):
        #c = gcd([ZZ(ci) for ci in (Phi_k(d*x+i)).list()])
        #print("Phi_k({}*x+{}) = {}, c = {} = {} mod {}".format(d, i, Phi_k(d*x+i), c, c%d, d))
        c = ZZ(Phi_k(i))
        if (c % d) == 0:
            u_m.append(i)
    return d, u_m

class MNTG(EllipticCurve_finite_field):
    """
    A generalised Miyaji-Nakabayashi-Takano curve of embedding degree k
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.20.8113&rep=rep1&type=pdf
    with the modification from Scott-Barreto (DCC06, preprint https://eprint.iacr.org/2004/058)
    where r = Phi_k(u)/d is prime and order = c*r = p+1-t (and t = u+1)
    The class is common to Scott-Barreto and Galbraith-McKee-Valenca with the correspondance that
    given t(x), then d = content(Phi_k(tx-1)) and rx = Phi_k(t-1)//d,
    then after that, h = content(px+1-tx), so
    rx = Phi_k(t-1)//d == (px+1-tx)//c
    """
    def __init__(self, k, u, D, c, d, a, b):
        """
        INPUT:
        - `k` embedding degree
        - `u` is the seed such that:
              t = u+1          (the trace)
              r = Phi_6(u)/d   (the curve has a subgroup of prime order r and embedding degree k)
              n = c*r          (n is the order of the curve, c is a tiny cofactor)
              p = n+u = n+t-1
        - `D`: discriminant of the curve, square-free integer s.t. t^2-4*p = -D*y^2
        - `c`: co-factor of the curve, # E(Fp) = c*r where r is prime
        - `d`: co-factor of Phi_k(u): r*d = Phi_k(u), hence r = Phi_k(t-1)/d where t = u+1
        - `a`,`b`: curve coefficients s.t. E: y^2=x^3+a*x+b has exactly order n=c*r
        """
        self._k = k # embedding degree
        self._c = Integer(c)
        self._d = Integer(d)
        self._a = Integer(a)
        self._b = Integer(b)
        self._D = Integer(D)
        # polynomials (will be needed later for polynomial selection)
        # R = Phi_k(x)/d = [1,-1,1]/d if k=6, [1,0,1]/d if k=4, [1,1,1]/d if k=3
        if self._k == 6:
            self.polynomial_r = [1, -1, 1]
        elif self._k == 4:
            self.polynomial_r = [1, 0, 1]
        elif self._k == 3:
            self.polynomial_r = [1, 1, 1]
        else:
            raise ValueError("Error k should be 3, 4 or 6 but k={} given".format(k))
        self.polynomial_r_denom = self._d
        self.polynomial_c = [self._c]
        self.polynomial_c_denom = 1
        # T = x+1
        self.polynomial_tr = [1, 1]
        self.polynomial_tr_denom = 1
        # P+1-T = c*R <=> P = c*R+T-1 = c/d*Phi_k(x) + [0,1] = (c*Phi_k + [0,d])/d
        g = gcd(c,d)
        c1 = Integer(c//g)
        d1 = Integer(d//g)
        # (c1*Phi_k + [0,d1])/d1
        self.polynomial_p = [c1*ri for ri in self.polynomial_r]
        self.polynomial_p[1] = self.polynomial_p[1] + d1
        self.polynomial_p_denom = Integer(d1)
        # no beta, lambda because de CM by sqrt(-D) is not trivial

        self._u = Integer(u)
        self._p = sum([self.polynomial_p[i]*self._u**i for i in range(len(self.polynomial_p))])//self.polynomial_p_denom
        self._pbits = self._p.nbits()
        self._tr = self._u + 1
        self._r = sum([self.polynomial_r[i]*self._u**i for i in range(len(self.polynomial_r))])//self.polynomial_r_denom
        if not self._c * self._r == self._p + 1 - self._tr:
            raise ValueError("Error: r*c != p+1-tr\nr={}\nc={}\np+1-tr={}\n".format(self._r,self._c,self._p+1-self._tr))
        #print("computing y")
        self._y = Integer(((4*self._p - self._tr**2)//self._D).sqrt())
        #print("y done")

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
        if self._pbits <= 512:
            self.order_checked = super(MNTG,self).order()
            if self.order_checked != (self._p+1-self._tr):
                raise ValueError("Wrong curve order: this one may be a twist")
        else:
            print("check order (randomized)")
            self.curve_order = self._p + Integer(1) - self._tr
            self.twist_order = self._p + Integer(1) + self._tr
            for i in range(10):
                P = self.random_element()
                if self.curve_order*P != self(0):
                    if self.twist_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a twist: (p+1+tr)*P = 0\ntr={}\nr={}\np+1+tr={}\n".format(self._tr,self.curve_order,self.twist_order))
                    else:
                        self.order_checked = super(MNTG,self).order()
                        raise ValueError("Wrong curve order:\np+1-tr        = {}\np+1+tr        = {}\nchecked order = {}\np             = {}".format(self.curve_order,self.twist_order,self.order_checked,self._p))
        #print("ok")
        
        # computes a generator
        self._G = get_curve_generator_order_r(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]

    def _repr_(self):
        return "MNT generalized p"+str(self._pbits)+" (Miyaji-Nakabayashi-Takano curve k modified Scott-Barreto and Galbraith-McKee-Valenca) curve with seed "+str(self._u)+"\n"+super(MNTG,self)._repr_()

    def u(self):
        return self._u
    def p(self):
        return self._p
    def r(self):
        return self._r
    def c(self):
        return self._c
    def d(self):
        return self._d
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
        return self.polynomial_p
    def rx(self):
        return self.polynomial_r
    def tx(self):
        return self.polynomial_tr
    def cx(self):
        return self.polynomial_c

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

