# http://doc.sagemath.org/html/en/thematic_tutorials/tutorial-objects-and-classes.html
# http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html
# http://doc.sagemath.org/html/en/thematic_tutorials/tutorial-implementing-algebraic-structures.html
# https://docs.python.org/3/library/functions.html#super
# https://rhettinger.wordpress.com/2011/05/26/super-considered-super/

import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.arith.misc import XGCD, xgcd
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
#from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
#from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0

allowed_k = [12]

def polynomial_params(k=12):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px = 36*x**4+36*x**3+24*x**2+6*x+1
    rx = 36*x**4+36*x**3+18*x**2+6*x+1
    tx = 6*x**2+1
    cx = QQx(1)
    yx=6*x**2+4*x+1
    betax=18*x**3+18*x**2+9*x+1
    lambx=36*x**3+18*x**2+6*x+1
    D=3
    return px, rx, tx, cx, yx, betax, lambx, D

def poly_cofactor_gt(k=12):
    """Computes the co-factors for GT: Phi_k(p(x))/r(x)
    It is always irreducible:
    from tnfs.curve.bn import poly_cofactor_gt
    k = 12
    print("BN")
    c = poly_cofactor_gt(k)
    print(["(deg {})^{} ".format(ci.degree(), ei) if ei>1 else "(deg {})".format(ci.degree()) for (ci, ei) in c.factor()])
    print(lcm([ci.denom() for ci in c.list()]))
    """
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    cx = cyclotomic_polynomial(k)(px)
    assert (cx % rx) == 0
    cx = cx // rx
    return cx

def coeff_params():
    Px = [1,6,24,36,36]
    Px_denom = 1
    Rx = [1,6,18,36,36]
    Rx_denom = 1
    Cx = [1] # cofactor such that p+1-tr == c*r
    Cx_denom = 1
    Tx = [1,0,6]
    Tx_denom = 1
    Yx = [1,4,6] # such that Tx^2 - 4*Px = -3*Yx^2
    Yx_denom = 1
    BETAx = [1,9,18,18]
    BETAx_denom = 1
    LAMBx = [1,6,18,36]
    LAMBx_denom = 1
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D

def poly_cofactor_twist_g1_g2(k=12):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params()
    twx = px+1+tx
    m = 1
    u_m = [0]
    # cx = 1 for BN curves
    #print("twx={}\n={}".format(twx,twx.factor()))
    # for G2, compute #E(Fp2) then compute its 6-th twist
    px2 = px**2
    tx2 = tx**2 - 2*px
    yx2 = yx*tx
    assert tx2**2 - 4*px2 == -D*yx2**2
    # now the 6-th twist that matches rx
    E2_order = px2+1-(3*yx2+tx2)/2
    assert (E2_order % rx) == 0
    g2cx = E2_order // rx # irreducible
    assert g2cx.is_irreducible()
    g2twx = px2+1+(3*yx2+tx2)/2
    #print("g2cx={}\n={}".format(g2cx,g2cx.factor()))
    #print("g2twx={}\n={}".format(g2twx,g2twx.factor()))
    g2_twxa = 36*x**4 + 36*x**3 + 18*x**2 + 1
    g2_twxb = 12*x**4 + 12*x**3 + 10*x**2 + 4*x + 1
    g2_twx0 = 3
    assert g2twx == g2_twx0 * g2_twxa * g2_twxb
    cofactor_r = [1]
    polys_cofact_twists = [twx, g2cx, g2_twxa, g2_twxb]
    label_factors = ["twx", "g2cx", "g2_twxa", "g2_twxb"]
    small_cofactors = [1, 1, 1, 1]

    return m, u_m, cofactor_r, twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors

class BN(EllipticCurve_finite_field):
    """
    A Barreto-Naehrig curve
    """
    def __init__(self, u, b=None, k=12, cofactor_r=1):
        """
        u is the seed such that p=P_BN(u), and b is the curve coefficient, s.t. y^2 = x^3 + b has exactly order r=R_BN(u)
        (there are two possible choices for b: b, 1/b (a twist) but only one has order r)
        """
        self._k = 12 # embedding degree
        self._a = 0 # first curve parameter is 0 because j=0
        self._D = 3
        # polynomials (will be needed later for polynomial selection)
        # P_BN = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1
        # R_BN = 36*u^4 + 36*u^3 + 18*u^2 + 6*u + 1
        # C_BN = 1 # cofactor such that p+1-tr == c*r
        # T_BN = 6*u^2 + 1
        # Y_BN = 6*u^2 + 4*u + 1 # such that T_BN^2 - 4*P_BN = -3*Y_BN^2
        # BETA_BN = 18*u^3 + 18*u^2 + 9*u + 1
        # LAMB_BN = 36*u^3 + 18*u^2 + 6*u + 1
        self.polynomial_p = [1,6,24,36,36]
        self.polynomial_r = [1,6,18,36,36]
        self.polynomial_c = [1]
        self.polynomial_tr= [1,0,6]
        self.polynomial_y = [1,4,6]
        self.polynomial_beta = [1,9,18,18]
        self.polynomial_lamb = [1,6,18,36]
        
        self._u = Integer(u)
        self._p = 36*self._u**4 + 36*self._u**3 + 24*self._u**2 + 6*self._u + 1
        self._pbits = self._p.nbits()
        self._tr = 6*self._u**2 + 1
        self._r = 36*self._u**4 + 36*self._u**3 + 18*self._u**2 + 6*self._u + 1
        self._c = Integer(1)
        self._y = 6*self._u**2 + 4*self._u + 1 # such that T^2 - 4*P = -3*Y^2
        # GLV parameter for fast scalar multiplication thanks to automorphism
        # psi: (x,y) -> (beta*x,y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+x+1 mod p, beta = (-1 +/- sqrt(Fp(-3)))/2
        # beta = 18*u^3 + 18*u^2 + 9*u + 1 (other one is -beta-1)
        # and lambda is a root of x^2+x+1 mod r, lamb = (-1 +/- sqrt(Fr(-3)))/2
        # lamb = 36*u^3 + 18*u^2 + 6*u + 1 (other one is -lamb-1)
        # there are two choices: beta, -beta-1, and the same for lamb: lamb, -lamb-1
        # arbitrarily the positive ones are chosen
        self._beta = 18*self._u**3 + 18*self._u**2 + 9*self._u + 1
        self._lamb = 36*self._u**3 + 18*self._u**2 + 6*self._u + 1
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
            # BN curves:
            # Pereira et al solution: b = c^4 + d^6 or b = c^6 + 4*d^4 for c, d in Fp*
            # https://eprint.iacr.org/2010/429
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
        self.order_checked = super(BN,self).order()
        if self.order_checked != self._r:
            print("Error, wrong order")
            raise ValueError("Wrong curve order: this one is a twist")

        # computes a generator
        self._G = get_curve_generator_order_r_j0(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]
        # note that there is no cofactor for BN curves because r=p+1-tr
        
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
        return "BN p"+str(self._pbits)+" (Barreto-Naehrig) curve with seed "+str(self._u)+"\n"+super(BN,self)._repr_()

    def u(self):
        return self._u
    def p(self):
        return self._p
    def r(self):
        return self._r
    def c(self):
        return self._c # cofactor is Integer(1) because r=p+1-tr and r is prime
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


