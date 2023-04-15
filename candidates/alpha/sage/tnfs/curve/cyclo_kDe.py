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
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.number_field.number_field import NumberField
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve as pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_a_j1728, get_curve_generator_order_r_j1728
from tnfs.curve.pairing_friendly_curve import get_curve_generator_order_r
import tnfs.curve.fst66
import tnfs.curve.fst64
import tnfs.curve.fst62
import tnfs.curve.fst67
import tnfs.curve.bn
import tnfs.curve.kss16
import tnfs.curve.fotiadis_martindale
#import tnfs.curve.bw6_761
#import tnfs.curve.aurifeuillean


def coeff_mult_m(k,D):
    m=1
    if (D%3) == 0:
        m=3//gcd(3,k)
    elif (D==4 or D==1):
        m=4//gcd(4,k)
    elif (D==2):
        m=8//gcd(8,k)
    # example: k=10 and D=5, with m=2, then -5 is a square mod Phi_20(x)
    elif (k==10) and (D==5):
        m=2
    else:
        m=1
    return m
    
def polynomial_params(k,D,e0,m=None):
    """Returns the parameters px,rx,cx,tx,yx,betax,lambx, and D of the Cyclotomic Construction, a generalisation of Brezing-Weng
    k: embedding degree
    D: CM discriminant (absolute value)
    e0: exponent of the trace in the Cyclotomic construction
    m: multiplicative factor so that r(x) = Phi_{k*m}(x)

    k=10, D=3, e0=1 is a special case, 
    it does not fall into the general scheme because in that particular case, 
    we have yx = rx + ((tx-2)/sqrt(-D) mod rx) and rho = 2."""

    if k==6 and D==3 and e0 == 0:
        return tnfs.curve.bw6_761.polynomial_params()
    if k==6 and D==3 and e0 == 1: # no reduction mod rx otherwise it is not irreducible
        return tnfs.curve.fst66.polynomial_params(6)
    if k==10 and D==3 and e0 == 1:
        return tnfs.curve.fst66.polynomial_params(10)

    if k==12 and (D==1 or D==4) and e0 == 1:
        return tnfs.curve.fst64.polynomial_params(12)
    if k==12 and D==3 and e0 == 0:
        return tnfs.curve.bn.polynomial_params()
    if k==16 and (D==1 or D==4) and e0 == 0:
        return tnfs.curve.kss16.polynomial_params()

    if k==15 and D==1 and e0==8:
        px, rx, tx, cx, yx, betax, lambx, D = tnfs.curve.fst62.polynomial_params(15)
        if D==4:
            D = 1
            yx = 2*yx
        return px, rx, tx, cx, yx, betax, lambx, D
    if k==15 and D==2 and e0==8:
        return tnfs.curve.fst67.polynomial_params(15)

    if k==8 and D==1 and e0 == 1:
        # the cyclo construction does not produce an irreducible px but the Fotiadis-Martindale curve with rho=2, tx = x+1+rx and yx = ((tx-2)/sqrt(-1) % rx) gives valid parameters
        return tnfs.curve.fotiadis_martindale.polynomial_params_from_code(1) # code "1" in eprint 2019/555
    if k==9 and D==3 and e0 == 5:
        # the cyclo construction does not produce an irreducible px but the Fotiadis-Martindale curve with rho=2, tx = x^5+1+rx and yx = ((tx-2)/sqrt(-3) % rx) + rx gives valid parameters
        return tnfs.curve.fotiadis_martindale.polynomial_params_from_code(10) # code "10" in eprint 2019/555
    if k==12 and D==3 and e0 == 17:
        return tnfs.curve.fotiadis_martindale.polynomial_params_from_code(17) # code "17" in eprint 2019/555
    if k==12 and D==3 and e0 == 19:
        return tnfs.curve.FotiadisMartindale.polynomial_params_from_code(19) # code "19" in eprint 2019/555
    if k==12 and D==3 and e0 == 20:
        return tnfs.curve.fotiadis_martindale.polynomial_params_from_code(20) # code "20" in eprint 2019/555
    if k==18 and D==3 and e0==0: # return Aurifeuillean family
        return tnfs.curve.aurifeuillean.polynomial_params(k=18,D=3)
        
    # example: k=10 and D=5, with m=2, then -5 is a square mod Phi_20(x)
    if (k==10) and (D==5):
        m0=2
    else:
        m0 = coeff_mult_m(k,D)
    if (m is None) or (m % m0) != 0:
        m = m0
    l = m*k
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    Phi_l = QQx(cyclotomic_polynomial(l))
    Phi_k = QQx(cyclotomic_polynomial(k))
        
    rx = Phi_l # x is a l-th root of unity, and we hope that (-D) is a square modulo rx
    # zeta_k = x
    tx = QQx(x**(e0*m)+1) % rx
    assert (Phi_k(tx-1) % rx) == 0, "Error k={}, Phi_k(tx-1)%rx != 0, = {}, tx=x^{}+1 mod rx = {}, rx={}".format(k, Phi_k(tx-1) % rx, e0*m, tx, rx)
    # there could be a default value
    # computes 1/sqrt(-D) = sqrt(-D)/D
    #K = CyclotomicField(l)
    K = NumberField(rx, names=('a',)); (a,) = K._first_ngens(1)
    inv_sqrt_D_ = (sqrt(K(-D))/K(D))
    inv_sqrt_D = QQx(inv_sqrt_D_.polynomial())
    assert ((D*inv_sqrt_D**2+1) % rx) == 0
    yx = ((tx-2)*inv_sqrt_D) % rx
    px = QQx(tx**2 + D*yx**2)/4
    px_denom = lcm([pi.denom() for pi in px.list()])
    assert ((px+1-tx) % rx) == 0, "Error (px+1-tx)%rx != 0, px={}*({}), tx={}, rx={}, (px+1-tx)={}".format(px_denom, px_denom*px, tx,rx, (px+1-tx).factor())
    cx = (px+1-tx) // rx
    assert (Phi_k(px) % rx) == 0, "Error k={}, Phi_k(px)%rx != 0, = {}".format(k, Phi_k(px) % rx)
    assert px.is_irreducible(), "Error px is not irreducible, px= {}".format(px.factor())
    assert rx.is_irreducible(), "Error rx is not irreducible, rx= {}".format(rx.factor())
    if k==9 and D==3 and (e0==1 or e0==4 or e0==7):
        rx = rx/QQ(3)
        cx = cx*3
    if k==27 and D==3 and (e0==1 or e0==19):
        rx = rx/QQ(3)
        cx = cx*3
    if (k==11) and (D==11) and ((e0 == 4) or (e0 == 8)):
        rx = rx/QQ(11)
        cx = cx*11
    if D==3:
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
        lambx = QQx(x**(l//3)) # so that if 3|k, then l=k and l//3 == k//3, otherwise l=3*k and x^k is a 3-rd root of unity mod rx=Phi_l(x)
        assert ((betax**2 + betax + 1) % px) == 0
        assert ((lambx**2 + lambx + 1) % rx) == 0
    elif (D==1) or (D==4):
        g,u,v = px.xgcd(tx)
        if g == 1:
            inv_tx = QQx(v)
        else:
            inv_tx = QQx(v)/QQ(g)
        # px = (tx^2 + 4*yx^2) / 4
        # tx^2 = -4*yx^2
        # -1 = ((2*yx)/tx)^2
        if D==4:
            betax = QQx((inv_tx*yx*2) % px)
        else:
            betax = QQx((inv_tx*yx) % px)
        #print("computed betax = {}".format(betax))
        lambx = QQx(x**(l//4))
        assert ((betax**2 + 1) % px) == 0
        assert ((lambx**2 + 1) % rx) == 0
    else:
        # computes sqrt(-D) or (-1+sqrt(-D))/2 mod r(x), p(x)
        lambx = sqrt(K(-D)).polynomial()
        #Kr.<w> = NumberField(px)
        Kp = NumberField(px, names=('ww',))#; (w,) = Kp._first_ngens(1)
        betax = sqrt(Kp(-D)).polynomial()
        assert (lambx**2 + D) % rx == 0
        assert (betax**2 + D) % px == 0
        #if (-D % 4) == 1:
        #    lambx = (-1 + lambx)/2
        #    betax = (-1 + betax)/2

    return px, rx, tx, cx, yx, betax, lambx, D

def coeffs_params(k,D,e0,m=None):
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k,D,e0,m=m)
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
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom
    
class Cyclo_kDe(EllipticCurve_finite_field):
    """
    A Cyclotomic curve of embedding degree k and discriminant D, and tx = x^(e0*m)+1, where m=3 if D=3 and k!=0 mod 3, m=4 if D=1,4 and k=1 mod 2, m=2 if D=1,4 and k=2mod 4, otherwise m=1
    """
    def __init__(self, k, D, e0, u, m=None, a=None, b=None, cofactor_r=1, verbose_init=False):
        """
        k embedding degree
        D discriminant
        u is the seed such that p=P(u), and a,b are the curve coefficients, s.t. y^2 = x^3 + a*x + b has a subgroup of order r=R(u)
        e0 in Cyclo construction (Brezing-Weng with r(x) a cyclotomic polynomial
        m such that r(x) = cyclotomic_polynomial(m*k)
        """
        self._k = k # embedding degree
        self._D = D
        m0 = coeff_mult_m(k,D)
        if (m is None) or (m % m0) != 0:
            self._m = m0
        else:
            self._m = m
        if D==3:
            self._a = 0 # first curve parameter is 0 because j=0
        elif D==1 or D==4:
            self._b = 0 # j=1728
        
        # polynomials (will be needed later for polynomial selection) in NFS and TNFS
        self._px, self._px_denom, self._rx, self._rx_denom, self._tx, self._tx_denom, self._cx, self._cx_denom, self._yx, self._yx_denom, self._betax, self._betax_denom, self._lambx, self._lambx_denom = coeffs_params(k,D,e0,m=m)
        self._u = Integer(u)
        self._p = sum([Integer(self._px[i])*self._u**i for i in range(len(self._px))])//self._px_denom
        self._pbits = self._p.nbits()
        self._r = sum([Integer(self._rx[i])*self._u**i for i in range(len(self._rx))])//self._rx_denom
        self._tr= sum([Integer(self._tx[i])*self._u**i for i in range(len(self._tx))])//self._tx_denom
        self._y = sum([Integer(self._yx[i])*self._u**i for i in range(len(self._yx))])//self._yx_denom
        self._c = sum([Integer(self._cx[i])*self._u**i for i in range(len(self._cx))])//self._cx_denom
        if not self._c * self._r == self._p + 1 - self._tr:
            if self._c * self._r == self._p + 1 + self._tr:
                raise ValueError("Error: r*c != p+1-tr but r*c = p+1+tr this is the quadratic twist\nr={}\nc={}\np+1-tr={}\n".format(self._r,self._c,self._p+1-self._tr))
            else:
                raise ValueError("Error: r*c != p+1-tr\nr={}\nc={}\np+1-tr={}\nr*c   ={}\n".format(self._r,self._c,self._p+1-self._tr, self._r*self._c))
        
        
        # GLV parameter for fast scalar multiplication thanks to automorphism
        # psi: (x,y) -> (beta*x,y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+x+1 mod p, beta = (-1 + sqrt(Fp(-3)))/2
        # and lambda is a root of x^2+x+1 mod r, lamb = (-1 + sqrt(Fr(-3)))/2
        # there are two choices: beta, -beta-1, and the same for lamb: lamb, -lamb-1
        # arbitrarily the positive ones are chosen
        if self._betax != 0 and self._lambx != 0:
            self._beta = sum([Integer(self._betax[i])*self._u**i for i in range(len(self._betax))]) # will divide by self._betax_denom later
            self._lamb = sum([Integer(self._lambx[i])*self._u**i for i in range(len(self._lambx))])
        else:
            self._beta = 0
            self._lamb = 0
        
        try:
            self._Fp = FiniteField(self._p)
        except ValueError as err:
            print("ValueError creating Fp: {}".format(err))
            print("p= {}".format(self._p))
            raise
        except:
            print("Error creating Fp")
            raise
        self._cofactor_r = Integer(cofactor_r)
        if cofactor_r > 1:
            rem_cofact_r = (self._r % self._cofactor_r)
            if (rem_cofact_r == 0):
                self._r //= self._cofactor_r
                self._c *= self._cofactor_r
            else:
                raise ValueError("Error cofactor of r given is {} but r % cofactor = {}".format(cofactor_r, rem_cofact_r))
        if not self._r.is_prime():
            raise ValueError("Error r is not prime")

        if (self._beta != 0) and (self._betax_denom != 1):
            self._betax_denom = Integer(self._betax_denom)
            if (self._beta % self._betax_denom) == 0:
                self._beta = self._beta // self._betax_denom
            else:
                self._beta = Integer(self._Fp(self._beta)/self._Fp(self._betax_denom))

        if (self._lamb != 0) and (self._lambx_denom != 1):
            self._lamb_denom = Integer(self._lambx_denom)
            if (self._lamb % self._lambx_denom) == 0:
                self._lamb = self._lamb // self._lambx_denom
            else:
                self._lamb = Integer(self._Fp(self._lamb)/self._Fp(self._lambx_denom))

        if self._D == 3:
            if ((self._beta**2 + self._beta + 1) % self._p) != 0:
                raise ValueError("Error beta^2 + beta + 1 != 0 mod p")
            if ((self._lamb**2 + self._lamb + 1) % self._r) != 0:
                raise ValueError("Error lamb^2 + lamb + 1 != 0 mod r")
        elif self._D == 1 or self._D == 4:
            if ((self._beta**2 + 1) % self._p) != 0:
                raise ValueError("Error beta^2 + 1 != 0 mod p")
            if ((self._lamb**2 + 1) % self._r) != 0:
                raise ValueError("Error lamb^2 + 1 != 0 mod r")

        self._Fpz = PolynomialRing(self._Fp, names=('z',))
        (self._z,) = self._Fpz._first_ngens(1)

        if (self._D==1 or self._D==4):
            if a is None:
                self._a, self._ap = get_curve_parameter_a_j1728(self._tr, self._y, self._p, self._Fp)
                self._bp = self._Fp(0) #  second curve parameter is 0 because j=1728
            else:
                try:
                    a = Integer(a)
                except:
                    raise
                self._a = a
                self._ap = self._Fp(a)
                self._bp = self._Fp(0) #  second curve parameter is 0 because j=1728
        elif self._D==3:
            if b is None:
                # check that beta = 2*U/(-3*V-U) before, where U=t/2, V = y/2 and 2V = 2 mod 3
                self._b, self._bp = get_curve_parameter_b_j0(self._tr, self._y, self._p, self._Fp)
                self._ap = self._Fp(0) # first curve parameter is 0 because j=0
            else:
                try:
                    b = Integer(b)
                except:
                    raise
                self._b = b
                self._bp = self._Fp(b)
                self._ap = self._Fp(0) #  first curve parameter is 0 because j=0
        else:
            if a is not None:
                try:
                    a = Integer(a)
                except:
                    raise
                self._a = a
                self._ap = self._Fp(a)
            else:
                raise ValueError("Error please provide input parameter a")
            if b is not None:
                try:
                    b = Integer(b)
                except:
                    raise
                self._b = b
                self._bp = self._Fp(b)
            else:
                raise ValueError("Error please provide input parameter b")
        if verbose_init:
            print("got a and b, now computing E([a,b])")
        # Now self._a and self._b are such that E: y^2 = x^3 + a*x + b has order r
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
        if verbose_init:
            print("curve E constructed, checking order")
        if not (D==1 or D==4 or D==3): # otherwise it is too slow
            if verbose_init:
                print("check order (randomized)")
            self.curve_order = self._p + Integer(1) - self._tr
            self.twist_order = self._p + Integer(1) + self._tr
            for i in range(10):
                P = self.random_element()
                if self.curve_order*P != self(0):
                    if self.twist_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a twist: (p+1+tr)*P = 0\ntr={}\np+1-t={}\np+1+tr={}\n".format(self._tr,self.curve_order,self.twist_order))
                    else:
                        if self._p.nbits() <= 256:
                            self.order_checked = super(Cyclo_kDe,self).order()
                        else:
                            self.order_checked = None
                        raise ValueError("Wrong curve order:\np+1-tr        = {}\np+1+tr        = {}\nchecked order = {}\np             = {}".format(self.curve_order,self.twist_order,self.order_checked,self._p))
            print("ok")
        else:
            self.order_checked = super(Cyclo_kDe,self).order()
            if self.order_checked != (self._p+1-self._tr):
                if verbose_init:
                    print("Error, wrong order")
                if self.order_checked == (self._p+1+self._tr):
                    raise ValueError("Wrong curve order: this one is the quadratic twist of order p+1+t")
                else:
                    raise ValueError("Wrong curve order: this one might be a twist (but not the quadratic twist)")

        # computes a generator
        if self._D == 3:
            self._G = get_curve_generator_order_r_j0(self)
        elif self._D == 1 or self._D == 4:
            self._G = get_curve_generator_order_r_j1728(self)
        else:
            self._G = get_curve_generator_order_r(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]
        
        # adjust beta and lamb according to the curve
        if self._D == 3:
            # do we have (beta*x,y) = lamb*(x,y)?
            if self([self._Gx*self._beta, self._Gy]) != self._lamb*self._G:
                if verbose_init:
                    print("adjusting beta, lambda")
                if self([self._Gx*(-self._beta-1), self._Gy]) == self._lamb*self._G:
                    self._beta = -self._beta-1
                    if verbose_init:
                        print("beta -> -beta-1")
                elif self([self._Gx*self._beta, self._Gy]) == (-self._lamb-1)*self._G:
                    self._lamb = -self._lamb-1
                    if verbose_init:
                        print("lamb -> -lamb-1")
                elif self([self._Gx*(-self._beta-1), self._Gy]) == (-self._lamb-1)*self._G:
                    self._beta = -self._beta-1
                    self._lamb = -self._lamb-1
                    if verbose_init:
                        print("lamb -> -lamb-1")
                        print("beta -> -beta-1")
                else:
                    raise ValueError("Error while adjusting beta, lamb: compatibility not found")
        elif self._D == 4 or self._D == 1:
            # adjust beta and lamb according to the curve
            # do we have (-x,beta*y) = lamb*(x,y)?
            if self([-self._Gx, self._Gy*self._beta]) != self._lamb*self._G:
                if verbose_init:
                    print("adjusting beta, lambda")
                if self([-self._Gx, self._Gy*(-self._beta)]) == self._lamb*self._G:
                    self._beta = self._p-self._beta
                    if verbose_init:
                        print("beta -> -beta")
                elif self([-self._Gx, self._Gy*self._beta]) == (-self._lamb)*self._G:
                    self._lamb = -self._lamb
                    if verbose_init:
                        print("lamb -> -lamb")
                elif self([-self._Gx, self._Gy*(-self._beta)]) == (-self._lamb)*self._G:
                    self._beta = self._p-self._beta
                    self._lamb = -self._lamb
                    if verbose_init:
                        print("lamb -> -lamb")
                        print("beta -> -beta")
                else:
                    raise ValueError("Error while adjusting beta, lamb: compatibility not found")


    def _repr_(self):
        return "Cyclo_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (Cyclotomic k="+str(self._k)+") curve with seed "+str(self._u)+"\n"+super(Cyclo_kDe,self)._repr_()

    def u(self):
        return self._u
    def T(self):
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
        return self._a
    def ap(self):
        return self._ap
    def b(self):
        return self._b
    def bp(self):
        return self._bp
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

    def poly_p(self):
        return self._px
    def poly_p_denom(self):
        return self._px_denom
    def poly_r(self):
        return self._rx
    def poly_r_denom(self):
        return self._rx_denom

    def miller_loop_length():
        return self._u**coeff_mult_m(self._k,self._D)
    
    def print_parameters(self):
        tnfs.curve.pairing_friendly_curve.print_parameters(self)
        
    def print_parameters_for_RELIC(self):
        tnfs.curve.pairing_friendly_curve.print_parameters_for_RELIC(self)

