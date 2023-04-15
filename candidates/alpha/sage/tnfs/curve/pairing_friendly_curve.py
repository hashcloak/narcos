import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.arith.functions import lcm
from sage.arith.misc import GCD, gcd
from sage.arith.misc import XGCD, xgcd
from sage.rings.integer import Integer
from sage.rings.integer_ring import Z, ZZ
from sage.rings.rational_field import Q, QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial

def compute_beta_lambda(px, rx, tx, yx, D):
    """
    Compute the eigenvalues beta mod px and lambda mod rx such that
    if D!=3 mod 4: beta^2 + D = 0 mod px, lambda^2 + D = 0 mod rx,
    if D=3 mod 4:  beta^2 + beta + (1+D)/4 = 0 mod px, i.e.,
    beta = (-1+sqrt(-D))/2 mod px, and
    lambda^2 + lambda + (1+D)/4 = 0 mod rx, i.e.,
    lambda = (-1+sqrt(-D))/2 mod rx.

    INPUT: polynomial parameters of a pairing-friendly elliptic curve.
    """
    QQx = QQ['x']; #(x,) = QQx._first_ngens(1)
    if D < 0:
        D = -D
    g, u, v = px.xgcd(yx) # u*px + v*yx == g
    inv_yx_modp = QQx(v)/QQ(g)
    g, u, v = rx.xgcd(yx)
    inv_yx_modr = QQx(v)/QQ(g)
    if (D % 4) != 3:
        betax = (tx * inv_yx_modp) % px
        if betax.leading_coefficient() < 0:
            betax = -betax
        assert ((betax**2 + D) % px) == 0
        lambx = ((tx-2)*inv_yx_modr) % rx
        if lambx.leading_coefficient() < 0:
            lambx = -lambx
        assert ((lambx**2 + D) % rx) == 0
    else: #  (-1+sqrt(-D))/2, minimal poly is x^2 + x + (1+D)//4
        betax = ((-1+(tx * inv_yx_modp))/2) % px
        if betax.leading_coefficient() < 0:
            betax = -betax-1
        assert ((betax**2 + betax + (1+D)//4) % px) == 0
        lambx = ((-1+(tx-2)*inv_yx_modr)/2) % rx
        if lambx.leading_coefficient() < 0:
            lambx = -lambx-1
        assert ((lambx**2 + lambx + (1+D)//4) % rx) == 0
    return betax, lambx

def poly_cofactor_twist_g1_g2(k: int, px, rx, tx, cx, yx, D):
    """Computes the curve co-factors for the twists"""
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    twx = px+1+tx
    # does it factor?
    if (D == 3) and (k % 6) == 0:
        # sextic twist
        d = 6
    elif (D == 3) and (k % 3) == 0:
        # cubic twist
        d = 3
    elif ((D == 1) or (D==4)) and (k % 4) == 0:
        # quartic twist
        d = 4
    elif (k % 2) == 0:
        # quadratic twist
        d = 2
    else:
        # no twist
        d = 1
    k1 = k // d
    # compute the curve order up to k1 before computing the d-twist order
    t1 = tx
    # (p+1-t)*(p+1+t) = p^2 + 1 + 2*p - t^2 = p^2 + 1 - (t^2 - 2*p)
    t2 = tx**2 - 2*px
    # tk = t1*t_{k-1} -p*t_{k-2}
    i = 3
    t_im1 = t2 # t_{i-1}
    t_im2 = t1 # t_{i-2}
    while i <= k1:
        t_i = t1*t_im1 - px*t_im2
        t_im2 = t_im1
        t_im1 = t_i
        i = i+1
    if k1 == 1:
        tx_k1 = t1
    elif k1 == 2:
        tx_k1 = t2
    else:
        tx_k1 = t_i
    px_k1 = px**k1
    if d==3 or d==6 or d==4:
        yx_k1_square = (tx_k1**2 - 4*px_k1)/(-D)
        lc = yx_k1_square.leading_coefficient()
        assert lc.is_square()
        yx_k1_square_monic = yx_k1_square / lc
        yx_k1_factors = yx_k1_square_monic.factor()
        yx_k1 = lc.sqrt()
        for fa, ee in yx_k1_factors:
            assert (ee % 2) == 0
            yx_k1 = yx_k1 * fa**(ee//2)
        assert yx_k1**2 == yx_k1_square
    else:
        yx_k1 = 1
    if d==3 or d==6:
        if d==6:
            E2_order = px_k1+1-(-3*yx_k1+tx_k1)/2
            E2_order_= px_k1+1-( 3*yx_k1+tx_k1)/2
            g2twx = px_k1+1+(-3*yx_k1+tx_k1)/2
            g2twx_= px_k1+1+( 3*yx_k1+tx_k1)/2
        elif d==3:
            E2_order = px_k1+1-(-3*yx_k1-tx_k1)/2
            E2_order_= px_k1+1-( 3*yx_k1-tx_k1)/2
            g2twx = px_k1+1+(-3*yx_k1-tx_k1)/2
            g2twx_= px_k1+1+( 3*yx_k1-tx_k1)/2
        if (E2_order % rx) != 0 and (E2_order_ % rx) == 0:
            E2_order = E2_order_
            g2twx = g2twx_
    elif d==4:
        if D==1:
            E2_order = px_k1 + 1 + yx_k1
            g2twx = px_k1 + 1 - yx_k1 # quadratic twist of G2
        elif D==4:
            E2_order = px_k1 + 1 + 2*yx_k1
            g2twx = px_k1 + 1 - 2*yx_k1 # quadratic twist of G2
        if (E2_order % rx) != 0 and (g2twx % rx) == 0:
            E2_order, g2twx = g2twx, E2_order
    elif d == 2:
        E2_order = px_k1 + 1 + tx_k1
        g2twx = px_k1 + 1 - tx_k1 # quadratic twist of G2
    else: # d==1
        assert d==1
        E2_order = px_k1 + 1 - tx_k1
        g2twx = px_k1 + 1 + tx_k1 # quadratic twist of G2

    assert (E2_order % rx) == 0
    g2cx = E2_order // rx # irreducible
    # do cx, twx, g2cx, g2twx factor?
    polys_cofact_twists = [cx, twx, g2cx, g2twx]
    label_factors = ["cx", "twx", "g2cx", "g2twx"]
    small_cofactors = [1, 1, 1, 1]
    return twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors

def get_curve_parameter_a_j1728(E_tr, E_y, E_p, E_Fp):
    # computes 'a' such that y^2 = x^3 + a*x / Fp has order p+1-tr
    # Rubin Silverberg, Choosing the correct elliptic curve in the CM method
    # https://eprint.iacr.org/2007/253
    # Algorithm 3.4
    # MATHEMATICS OF COMPUTATION
    # Volume 79, Number 269, January 2010, Pages 545--561
    # https://www.ams.org/journals/mcom/2010-79-269/S0025-5718-09-02266-2
    # U = tr/2, V = y/2, d=1, p = (t^2 + y^2)/4 but in KSS16.py, D=4 and y <-> y/2
    U = E_tr // 2
    dy2_4 = E_p - U**2
    E_y2 = E_y**2
    if dy2_4 == E_y2: # D=4
        V = E_y
    elif 4*dy2_4 == E_y2: # D=1
        V = E_y//2
    if (U % 2) == 1 and ((U-1) % 4  == (V % 4)):
        return Integer(-1), E_Fp(-1)
    if (U % 2) == 1 and ((U-1) % 4  != (V % 4)):
        # return any 'a' a square but not a fourth power : a^((p-1)/2) == 1, a^((p-1)/4) == -1 mod p
        # note that -1 is a square because p = 1 mod 4 (otherwise the curve would be supersingular)
        a = 1
        ap = E_Fp(1)
        exponent = (E_p-1)//4
        E_Fp_1 = E_Fp(-1)
        while not (ap**exponent == E_Fp_1):
            a += 1
            ap += 1
        return Integer(-a), -ap
    if (U % 2) == 0:
        if ((V-1) % 4 != (U % 4)):
            V = -V
        a = 1
        ap = E_Fp(1)
        #while ap**((E_p-1)/4) != (E_Fp(U)/E_Fp(V)):
        exponent = (E_p-1)//4
        val = E_Fp(U)/E_Fp(V)
        while ap**exponent != val:
            a += 1
            ap += 1
        return Integer(-a), -ap

def get_curve_parameter_b_j0(E_tr, E_y, E_p, E_Fp):
    # computes 'b' such that y^2 = x^3 + b / Fp has order p+1-tr
    # that is, find b not a square, not a cube
    # How to choose the correct twist, knowing only b, without computing
    # the order of the curve? Without even defining a curve?
    # Rubin Silverberg, Choosing the correct elliptic curve in the CM method
    # https://eprint.iacr.org/2007/253
    # Algorithm 3.5
    # MATHEMATICS OF COMPUTATION
    # Volume 79, Number 269, January 2010, Pages 545--561
    # https://www.ams.org/journals/mcom/2010-79-269/S0025-5718-09-02266-2
    # U = tr/2, V = y/2 so that p = U^2 + 3*V^2 (4*p = t^2 + d*y^2)
    UU = E_tr # 2U
    VV = E_y  # 2V
    # we should have 4*p = E_tr^2 + 3*E_y^2
    # hence p = (E_tr/2)^2 + 3*(E_y/2)^2
    if (VV % 3) == 0 and ((UU % 3) == 2):
        E_b = Integer(16)
        E_bp = E_Fp(16)
    elif (VV % 3) == 0 and ((UU % 3) == 1):
        # b should be a cube but not a square
        # if p = 3 mod 4, -1 is not a square
        # if p = 1,7 mod 8, 2 is a square
        E_b = Integer(-1)
        E_bp = E_Fp(-1)
        exponent = (E_p-1)//6
        E_Fp_1 = E_Fp(-1)
        while (E_bp**exponent != E_Fp_1):
            if E_b < 0:
                E_b = -E_b + 1
                E_bp = -E_bp + 1
            else:
                E_b = -E_b
                E_bp = -E_bp
        E_b = E_b*16
        E_bp = E_bp*16
    else: # we cannot have 3|UU and p=(UU/2)^2 + 3*(VV/2)^2 prime
        if (VV % 3) == 2:
            VV = -VV # so that now VV = 1 mod 3
        if (UU % 3) == 2:
            # computes val=2U/(3V-U) mod p and find b: b^((p-1)/6)=2U/(3V-U) mod p
            # val = 2*(2U)/(3*2V-2U) = 2*E_tr/(3*E_y-E_tr)
            val = E_Fp(UU+UU)/E_Fp(3*VV-UU)
            E_b = Integer(1)
            E_bp = E_Fp(1)
            exponent = (E_p-1)//6
            while (E_bp**exponent != val):
                if E_b < 0:
                    E_b = -E_b + 1
                    E_bp = -E_bp + 1
                else:
                    E_b = -E_b
                    E_bp = -E_bp
            E_b = E_b*16
            E_bp = E_bp*16
        else:
            # computes val=2U/(3V+U) mod p and find b: b^((p-1)/6) = 2U/(3V+U) mod p
            # val = 2*(2U)/(3*2V+2U) = 2*E_tr/(3*E_y+E_tr)
            val = E_Fp(UU+UU)/E_Fp(3*VV+UU)
            E_b = Integer(1)
            E_bp = E_Fp(1)
            exponent = (E_p-1)//6
            while (E_bp**exponent != val):
                if E_b < 0:
                    E_b = -E_b + 1
                    E_bp = -E_bp + 1
                else:
                    E_b = -E_b
                    E_bp = -E_bp
            E_b = E_b*16
            E_bp = E_bp*16
    return E_b, E_bp

def get_curve_generator_order_r(E):
    E_ap = E.ap()
    E_bp = E.bp()
    E_c = E.c()
    E_r = E.r()
    E_p = E.p()
    E_Fp = E.Fp()
    E_Fpz, z = E.Fpz()

    Gx = E_Fp(0)
    Gfound = False
    while not Gfound:
        Gx += 1
        Gy2 = Gx**3 + E_ap * Gx + E_bp
        while not Gy2.is_square():
            Gx += 1
            Gy2 = Gx**3 + E_ap * Gx + E_bp
        Gy = Gy2.sqrt()
        E_G = E([Gx,Gy])
        # needs to take into account the cofactor
        E_G = E_c*E_G # so that G has exactly order r
        Gfound = (E_G != E(0) and E_r*E_G == E(0))
    return E_G
    
def get_curve_generator_order_r_j1728(E):
    E_ap = E.ap()
    E_c = E.c()
    E_r = E.r()
    E_p = E.p()
    E_Fp = E.Fp()
    E_Fpz, z = E.Fpz()

    Gx = E_Fp(0)
    Gfound = False
    while not Gfound:
        Gx += 1
        Gy2 = Gx**3 + E_ap * Gx
        while not Gy2.is_square():
            Gx += 1
            Gy2 = Gx**3 + E_ap * Gx
        Gy = Gy2.sqrt()
        E_G = E([Gx,Gy])
        # needs to take into account the cofactor
        E_G = E_c*E_G # so that G has exactly order r
        Gfound = (E_G != E(0) and E_r*E_G == E(0))
    return E_G
    
def get_curve_generator_order_r_j0(E):
    E_bp = E.bp()
    E_c = E.c()
    E_r = E.r()
    E_p = E.p()
    E_Fp = E.Fp()
    E_Fpz, z = E.Fpz()
    
    Gy = E_Fp(0)
    Gfound = False
    while not Gfound:
        Gy += 1
        Gx3 = Gy**2 - E_bp
        #while len((z**3 - Gx3).roots()) == 0 :
        while (Gx3**((E_p-1)//3) != E_Fp(1)):
            Gy += 1
            Gx3 = E_Fp(Gy**2 - E_bp)
        Gx = (z**3 - Gx3).roots()[0][0]
        E_G = E([Gx,Gy])
        # needs to take into account the cofactor
        E_G = E_c*E_G # so that G has exactly order r
        Gfound = (E_G != E(0) and E_r*E_G == E(0))
    return E_G

def print_parameters(E):
    # https://docs.python.org/2/library/string.html#format-string-syntax
    print("u ={: d}".format(E._u))
    print("p ={: d}".format(E._p))
    print("r ={: d}".format(E._r))
    print("t ={: d} # trace".format(E._tr))
    print("y ={: d} # s.t. p=(t^2+D*y^2)/4".format(E._y))
    print("c ={: d} # cofactor".format(E._c))
    print("u ={: #x}".format(E._u))
    print("p ={: #x}".format(E._p))
    print("r ={: #x}".format(E._r))
    print("t ={: #x} # trace".format(E._tr))
    print("y ={: #x} # s.t. p=(t^2+D*y^2)/4".format(E._y))
    print("c ={: #x} # cofactor".format(E._c))
    print("log_2 p   ={0:8.2f}, p   {1:5d} bits".format(float(log(E._p,2)), E._p.nbits()))
    print("log_2 p^k ={0:8.2f}, p^k {1:5d} bits".format(float(E._k*log(E._p,2)), (E._p**E._k).nbits()))
    print("log_2 r   ={0:8.2f}, r   {1:5d} bits, {2:2d} machine words of 64 bits".format(float(log(E._r,2)), E._r.nbits(), ceil(log(E._r,2)/64)))

def print_parameters_for_RELIC(E):
    print("p in hexa, on 64-bit words, lowest significant word first.")
    i=0
    for pi in E._p.digits(2**64):
        # fill at 16 hexa digits, + 2 for the 'prefix '0x' makes 18 characters
        print("P{0} {1:0=-#018x}".format(i, pi))
        i+=1
    g,u0,v0 = xgcd(E._p,2**64)
    if abs(u0) > abs(v0):
        u0,v0 = v0,u0
    if u0 < 0:
        u0 = -u0; v0 = -v0
    print("U0 {:0=-#018x}".format(u0))
    #params for the ep module in RELIC

    print("a    ={:= #x} #({})".format(Integer(E._a), E._a))
    print("b    ={:= #x} #({})".format(Integer(E._b), E._b))
    if hasattr(E, '_beta'):
        print("BETA ={:= #x}".format(E._beta))
    if hasattr(E, '_lamb'):
        print("LAMB ={:= #x}".format(E._lamb))
    
    print("Gx ={:= #x}".format(Integer(E._Gx)))
    print("Gy ={:= #x}".format(Integer(E._Gy)))



def check_order(Fp, p, tr, a, b):
    curve_order = p+1-tr
    twist_order = p+1+tr
    E = EllipticCurve([0,0,0,Fp(a),Fp(b)])
    E0 = E(0)
    ok = True
    twist_ok = False
    for i in range(10):
        P = E.random_element()
        if curve_order*P != E0:
            ok = False
            if twist_order*P == E0:
                twist_ok = True
            break
    return ok, twist_ok

HD_dict = {
    5 : [-681472000, -1264000, 1], #H(-20)
    6 : [14670139392, -4834944, 1], #H(-24)
    10: [9103145472000, -425692800, 1], #H(-40)
    13: [-567663552000000, -6896880000, 1], #H(-52)
    14: [10064086044321563803648, 2257767342088912896, 2059647197077504, -16220384512, 1], #H(-56)
    15: [-121287375, 191025, 1], #H(-15)
    17: [-2089297506304000000000000, -318507038720000000000, -75843692160000000, -178211040000, 1], #H(-68)
    21: [-5133201653210986057826304, 88821246589810089394176, -5663679223085309952, -3196800946944, 1], #H(-84)
    22: [15798135578688000000, -6294842640000, 1], #H(-88)
    23: [12771880859375, -5151296875, 3491750, 1], #H(-23)
    26: [65437179730333545242323676123103232, -25735039642229334200564710375424, 1378339984770204584193868955648, 31013571054009020830449664, 739545196164376195072, -82028232174464, 1], #H(-104)
    29: [-100730316193548175256338136121783353344, 143376986667050616958401264069115904, -66527716583835083670963399688192, -835102260960042427461140480, -11056847669496432594944, -495202728828032, 1], #H(-116)
    30: [4934510722321469030006784000000, -2588458316335175909376000000, 26329406807264910336000, -883067971104000, 1], #H(-120)
    31: [1566028350940383, -58682638134, 39491307, 1], #H(-31)
    33: [1656636925108948992000000000000, 54984539729717250048000000000, -325211610485778048000000, -4736863498464000, 1], #H(-132)
    34: [2422829169428572504087521656832, -1834607111282472051029311488, 735960027609078992953344, -8151279336430848, 1], #H(-136)
    35: [-134217728000, 117964800, 1], #H(-35)
    37: [-7898242515936467904000000, -39660183801072000, 1], #H(-148)
    38: [472390748138731280269312000000000000000000, -1380504171426125758791680000000000000000, 2783058624787093614292992000000000000, 6854544294799483688960000000000, 17024071380555203520000000, -66246265919280000, 1], #H(-152)
    39: [20919104368024767633, 109873509788637459, -429878960946, 331531596, 1], #H(-39)
    41: [-852636173252919999445568788749874942641540406706176, 716292304882512928715138362472485709784740265984, 58876580988711431943771690012552346623541248, -82923859178811827895415538602091992842240, -72018009354152588972347870534871023616, -107971538556531472065498397540352, -161936389233870440957755392, -296853791160440320, 1], #H(-164)
    42: [496644064976895846912000000000000000, -264691184105480095991808000000000, 336511679671210230144000000, -483435712076832000, 1], #H(-168)
    46: [114574710497270997578522590458150912, 38705419208160503264676104110080, 5767007465145198439020847104, -3215890895076912384, 1], #H(-184)
    47: [16042929600623870849609375, -14982472850828613281250, 5115161850595703125, -9987963828125, 2257834125, 1], #H(-47)
    51: [6262062317568, 5541101568, 1], #H(-51)
    53: [-67450134022842979455115194007552000000000000000000, 39924086528997881772669622484992000000000000000, -11008353578715780277672803110912000000000000, -2630171369254890916959016960000000000, -628986407384453487358016000000, -73387074029381328000, 1], #H(-212)
}

def compute_a_b(D, p, tr, y=None, verbose=False):
    if D < 0:
        D = -D
    p = Integer(p)
    tr = Integer(tr)
    Fp = GF(p, proof=False)
    if D == 1 or D == 4 or D == 3:
        if y == None:
            assert ((4*p-tr**2) % D) == 0
            y = ZZ(sqrt((4*p-tr**2)//D))
        else:
            y = Integer(y)
        assert tr**2+D*y**2 == 4*p
    if D==1 or D==4:
        a, ap = get_curve_parameter_a_j1728(tr, y, p, Fp)
        b = 0
        return (a,b)
    if D==3:
        b, bp = get_curve_parameter_b_j0(tr, y, p, Fp)
        a = 0
        return (a,b)

    # a = -3*x/(x-1728) ; b = 2*x/(x-1728)
    if D in [2, 7, 11, 19, 43, 67, 163]:
        if D==2: # j=8000=20^3 hilbert_class_polynomial(-8) -> x - 8000
            (a,b) = (-30, 56) # (-15/2, 7)
        elif D==7: # j=-3375=-15^3 hilbert_class_polynomial(-7) -> x + 3375
            (a,b) = (-35,-98) # (-5/7, 2/7)
        elif D==11: #j=-32768=-2^15 hilbert_class_polynomial(-11) -> x + 32768
            (a,b) = (-264, 1694) # (-6/11, 7/44)
        elif D==19: # j=-884736 = -2^15*3^3 hilbert_class_polynomial(-19) -> x + 884736
            (a,b) = (-152, 722) # (-2/19, 1/(4*19)) = (-8*19, 2*19*19)
        elif D==43: # x + 884736000
            (a,b) = (-3440, 77658) # (-20/43, 21/(4*43))
        elif D==67: # x + 147197952000
            (a,b) = (-29480, 1948226) # (-110/67, 217/268)
        elif D==163: # x + 262537412640768000
            (a,b) = (-8697680, 9873093538) # (-53360/163, 371602/163)
        # adjust (a,b) in case this is a twist
        Fp = GF(p)
        ok, twist_ok = check_order(Fp, p, tr, a, b)
        # usually the problem comes from E being actually the quadratic twist of order (p+1+t) instead of the curve itself of order (p+1-t)
        # the quadratic twist has equation (w^4*a, w^6*b) where w is neither a square nor a cube in Fp
        # find such w
        w0 = ZZ(-1)
        w0p = Fp(-1)
        a0, b0 = a, b
        while not ok:
            while w0.is_square() or w0p.is_square():
                if w0 >= 0:
                    w0 = -w0
                    w0p = -w0p
                else:
                    w0 = -w0 + 1
                    w0p = -w0p + 1
            print("trying with a={}^2*a, b={}^3*b".format(w0,w0))
            a = a0*w0**2
            b = b0*w0**3
            ok, twist_ok = check_order(Fp, p, tr, a, b)
            if not ok:
                if w0 >= 0:
                    w0 = -w0
                    w0p = -w0p
                else:
                    w0 = -w0 + 1
                    w0p = -w0p + 1
        ab = (a,b)
        return ab
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    if D in HD_dict:
        H = Fpz(HD_dict[D])
    else:
        if ((-D) % 4) == 2 or ((-D) % 4) == 3:
            D = 4*D
        H = (hilbert_class_polynomial(-D)).change_ring(Fp)
    if verbose:
        print("#    H(-D).degree() is {}".format(H.degree()))
    rH = H.roots()
    # can we have a=-3?
    ab = None
    a_wanted = -3
    while ab is None:
        print("#looking for a={}".format(a_wanted))
        j_idx = 0
        while ab is None and j_idx < len(rH):
            # y^2 + xy = x^3 - (36/(j - 1728))*x - 1/(j - 1728)
            # (y+1/2*x)^2 - 1/4*x^2 = x^3  - (36/(j - 1728))*x - 1/(j - 1728)
            # y' = y+1/2*x
            # y'^2 = x^3 + 1/4*x^2  - (36/(j - 1728))*x - 1/(j - 1728)
            # x' = (x-1/12)
            # y'^2 = x^3 - 1/48*j/(j - 1728)*x + 1/864*j/(j - 1728)
            # (a,b) ~ (a/u^4, b/u^6)
            j = Fp(rH[j_idx][0])
            j_idx += 1
            if j == 0:
                print("unexpected j=0, a=0")
                continue
            elif j == 1728:
                print("unexpected j=1728, b=0")
                continue
            #fac_j = Fp(j)/Fp(j - 1728)
            #a = -fac_j/Fp(48)
            #b = fac_j/Fp(864)
            a = -3*j*(j - 1728)
            b = 2*j*(j - 1728)**2
            # find u s.t. (a/u^2, b/u^3) == (a_wanted, b_other)
            u4 = a/Fp(a_wanted)
            if not u4.is_square():
                continue
            u2=u4.sqrt()
            # two options here: (a_wanted, b/u2**3) and (a_wanted, -b/u2**3)
            a = a_wanted
            b = b/u2**3
            # check that this is the curve and not the twist
            ok, twist_ok = check_order(Fp, p, tr, a, b)
            if not ok:
                b = -b
                ok, twist_ok = check_order(Fp, p, tr, a, b)
                if not ok:
                    continue
            b_ZZ = ZZ(b)
            if abs(b_ZZ-p) < abs(b_ZZ):
                b_ZZ -= p
            ab = (a, b_ZZ)
        if a_wanted == -3:
            a_wanted = 1
        else:
            a_wanted += 1
    if ab is None:
        j, je = rH[0]
        a = ZZ(-3*j*(j - 1728))
        b = ZZ(2*j*(j - 1728)**2)
        if abs(a-p) < abs(a):
            a = a-p
        if abs(b-p) < abs(b):
            b = b-p
        j = ZZ(j)
        if abs(j-p) < abs(j):
            j = j-p
        print("j={}".format(j))
        print("a,b = {},{}".format(a,b))
        # simplify
        g = gcd(a,b)
        gf = g.factor()
        u = 1
        for fac,ei in gf:
            if fac > 1:
                while (a % (fac**2)) == 0 and (b%(fac**3)) == 0:
                    a = a // fac**2 # 4
                    b = b // fac**3 # 6
                    u *= fac
        print("a={}, b= {}, simplified by u = {}".format(a,b,u))
        if abs((a % p)) < abs(a):
            a = a % p
            if abs(a-p) < abs(a):
                a = a-p
        if abs((b % p)) < abs(b):
            b = b % p
            if abs(b-p) < abs(b):
                b = b-p
        # check that this is the curve and no the twist
        ok, twist_ok = check_order(Fp,p,tr,a,b)
        if not ok and not twist_ok:
            print("pb with a,b")
            return None
        ns = -2
        nsp= Fp(-2)
        while not ok and twist_ok: # this is the twist
            # find a non-square
            ns += 1
            nsp += 1
            while nsp.is_square():
                ns += 1
                nsp += 1
            a = a*ns**2
            b = b*ns**3
            ok, twist_ok = check_order(Fp,p,tr,a,b)
            print("ns = {}, try with (a,b) = ({},{})".format(ns,a,b))
        ab = (a,b)
    return ab

def eval_polynomial_Horner(px, u):
    l = len(px)
    if l == 0:
        return 0
    if l == 1:
        return px[0]
    l = l-1
    p = px[l]*u + px[l-1]
    l = l-2
    while l >= 0:
        p = u*p + px[l]
        l = l-1
    return p


class BrezingWeng(EllipticCurve_finite_field):
    """
    A Brezing--Weng pairing-friendly curve, for inheritance of families
    """

    def __init__(self, k, D, u, px, px_denom, rx, rx_denom, tx, tx_denom, cx, cx_denom, yx, yx_denom, betax=None, betax_denom=None, lambx=None, lambx_denom=None, a=None, b=None, cofactor_r=1, verbose_init=False):
        """
        Create an elliptic curve E: y^2=x^3+a*x+b / Fp

        INPUT:
        - `k` embedding degree
        - `D` discriminant (1, 3, etc)
        - `u` seed, integer
        - `a` curve coefficient
        - `b` curve coefficient
        - `cofactor_r` integer, cofactor s.t. r=rx(u)/rx_denom / cofactor_r is prime
        - `px` integer coefficients of polynomial p(x)
        - `px_denom` integer, denominator s.t. p=px(u)/px_denom
        - `rx` integer coefficients of polynomial r(x)
        - `rx_denom` integer, denominator s.t. r=rx(u)/rx_denom
        - `tx` integer coefficients of polynomial t(x)
        - `tx_denom` integer, denominator s.t. t=tx(u)/tx_denom
        - `cx` integer coefficients of polynomial c(x)
        - `cx_denom` integer, denominator s.t. c=cx(u)/cx_denom
        - `yx` integer coefficients of polynomial y(x)
        - `yx_denom` integer, denominator s.t. y=yx(u)/yx_denom
        - `betax` integer coefficients of polynomial beta(x)
        - `betax_denom` integer, denominator s.t. beta=betax(u)/betax_denom
        - `lambx` integer coefficients of polynomial lambda(x)
        - `lambx_denom` integer, denominator s.t. lambda=lambx(u)/lambx_denom

        where p is the field characteristic, r is the prime subgroup order,
        t is the curve trace, c is the cofactor s.t. p+1-t = c*r,
        y is s.t. p = (t^2+D*y^2)/4, in other terms t^2-4*p = -D*y^2 and
        -D is the square-free part,
        beta is a root mod p of either x^2+x+(D+1)/4 if D = 3 mod 4,
        or a root of x^2+D otherwise,
        lambda is a root mod r of either x^2+x+(D+1)/4 if D = 3 mod 4,
        or a root of x^2+D otherwise.
        """
        self._k = k # embedding degree
        self._D = D # curve discriminant
        if D==3:
            self._a = 0 # first curve parameter is 0 because j=0
        elif D==1 or D==4:
            self._b = 0 # j=1728
        self._px, self._px_denom = px, Integer(px_denom)
        self._rx, self._rx_denom = rx, Integer(rx_denom)
        self._tx, self._tx_denom = tx, Integer(tx_denom)
        self._cx, self._cx_denom = cx, Integer(cx_denom)
        self._yx, self._yx_denom = yx, Integer(yx_denom)

        self._u = Integer(u)
        # p
        self._p = eval_polynomial_Horner(self._px, self._u)
        if (px_denom != 1):
            if (self._p % self._px_denom) == 0:
                self._p = self._p//self._px_denom
            else:
                raise ValueError("Error px(u) % px_denom = {} != 0 with u={:#x}, px={}, px_denom={}".format(self._p % self._px_denom, self._u, self._px, self._px_denom))
        self._pbits = self._p.nbits()
        # r
        self._r = eval_polynomial_Horner(self._rx, self._u)
        if (rx_denom != 1):
            if (self._r % self._rx_denom) == 0:
                self._r = self._r//self._rx_denom
            else:
                raise ValueError("Error rx(u) % rx_denom = {} != 0 with u={:#x}, rx={}, rx_denom={}".format(self._r % self._rx_denom, self._u, self._rx, self._rx_denom))
        # tr
        self._tr = eval_polynomial_Horner(self._tx, self._u)
        if (tx_denom != 1):
            if (self._tr % self._tx_denom) == 0:
                self._tr = self._tr//self._tx_denom
            else:
                raise ValueError("Error tx(u) % tx_denom = {} != 0 with u={:#x}, tx={}, tx_denom={}".format(self._tr % self._tx_denom, self._u, self._tx, self._tx_denom))
        # c
        self._c = eval_polynomial_Horner(self._cx, self._u)
        if (cx_denom != 1):
            if (self._c % self._cx_denom) == 0:
                self._c = self._c//self._cx_denom
            else:
                raise ValueError("Error cx(u) % cx_denom = {} != 0 with u={:#x}, cx={}, cx_denom={}".format(self._c % self._cx_denom, self._u, self._cx, self._cx_denom))
        # y
        self._y = eval_polynomial_Horner(self._yx, self._u)
        if (yx_denom != 1):
            if (self._y % self._yx_denom) == 0:
                self._y = self._y//self._yx_denom
            else:
                raise ValueError("Error yx(u) % yx_denom = {} != 0 with u={:#x}, yx={}, yx_denom={}".format(self._y % self._yx_denom, self._u, self._yx, self._yx_denom))
        # beta
        if betax is None or lambx is None:
            # convert coefficient list to polynomial ring
            QQx = QQ['x']; (x,) = QQx._first_ngens(1)
            px_ = QQx(px)/Integer(px_denom)
            rx_ = QQx(rx)/Integer(rx_denom)
            tx_ = QQx(tx)/Integer(tx_denom)
            yx_ = QQx(yx)/Integer(yx_denom)
            assert 4*px_ == tx_**2 + D*yx_**2
            assert ((px_+1-tx_) % rx_) == 0
            assert (((tx_-2)**2 + D*yx_**2) % rx_) == 0
            betax_, lambx_ = compute_beta_lambda(px_, rx_, tx_, yx_, D)
            betax_denom = Integer(lcm([ci.denom() for ci in betax_.list()]))
            if betax_denom != 1:
                betax = [Integer(ci) for ci in (betax_denom*betax_).list()]
            else:
                betax = [Integer(ci) for ci in betax_.list()]
            lambx_denom = Integer(lcm([ci.denom() for ci in lambx_.list()]))
            if lambx_denom != 1:
                lambx = [Integer(ci) for ci in (lambx_denom*lambx_).list()]
            else:
                lambx = [Integer(ci) for ci in lambx_.list()]
        #
        self._betax, self._betax_denom = betax, Integer(betax_denom)
        self._lambx, self._lambx_denom = lambx, Integer(lambx_denom)
        #
        try:
            self._Fp = FiniteField(self._p)
        except ValueError as err:
            print("ValueError creating Fp: {}".format(err))
            print("self._p = {} of type {}".format(self._p, type(self._p)))
            print("self._px = {} self._px_denom = {}".format(self._px, self._px_denom))
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
            raise ValueError("Error r = {} is not prime".format(self._r))
        # beta
        self._beta = eval_polynomial_Horner(self._betax, self._u)
        if (betax_denom != 1):
            if (self._beta % self._betax_denom) == 0:
                self._beta = self._beta//self._betax_denom
            else: # compute it mod p
                # res_p = self._p % self._betax_denom
                # residue + res_p*x = 0 mod betax_denom
                self._beta = Integer(self._Fp(self._beta)/self._Fp(self._betax_denom))
        else:
            self._beta = Integer(self._beta)
        # lambda
        self._lamb = eval_polynomial_Horner(self._lambx, self._u)
        if (lambx_denom != 1):
            if (self._lamb % self._lambx_denom) == 0:
                self._lamb = self._lamb//self._lambx_denom
            else:
                self._lamb = Integer(self._Fp(self._lamb)/self._Fp(self._lambx_denom))
        else:
            self._lamb = Integer(self._lamb)
        #
        if (self._D % 4) == 3:
            if ((self._beta**2 + self._beta + (self._D+1)//4) % self._p) != 0:
                raise ValueError("Error beta^2 + beta + (D+1)/4 != 0 mod p")
            if ((self._lamb**2 + self._lamb + (self._D+1)//4) % self._r) != 0:
                raise ValueError("Error lamb^2 + lamb + (D+1)/4 != 0 mod r")
        else:
            if ((self._beta**2 + self._D) % self._p) != 0:
                raise ValueError("Error beta^2 + D != 0 mod p")
            if ((self._lamb**2 + self._D) % self._r) != 0:
                raise ValueError("Error lamb^2 + D != 0 mod r")
        # Fp
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
            if a is not None and b is not None:
                try:
                    a = Integer(a)
                except:
                    raise
                self._a = a
                self._ap = self._Fp(a)
                try:
                    b = Integer(b)
                except:
                    raise
                self._b = b
                self._bp = self._Fp(b)
            else:
                a, b = compute_a_b(self._D, self._p, self._tr, self._y)
                self._a = Integer(a)
                self._ap = self._Fp(a)
                self._b = Integer(b)
                self._bp = self._Fp(b)
        if verbose_init:
            print("got a and b, now computing E([a,b])")
        # Now self._a and self._b are such that E: y^2 = x^3 + a*x + b has order r
        try:
            # this init function of super inherits from class EllipticCurve_generic defined in ell_generic.py
            # __init__ method inherited from ell_generic
            #EllipticCurve_finite_field.__init__(self, self._Fp, [0,0,0,self._ap,self._bp])
            super().__init__(self._Fp, [0,0,0,self._ap,self._bp])
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
                            #self.order_checked = super(BrezingWeng, self).order()
                            self.order_checked = super().order()
                        else:
                            self.order_checked = None
                        raise ValueError("Wrong curve order:\np+1-tr        = {}\np+1+tr        = {}\nchecked order = {}\np             = {}".format(self.curve_order,self.twist_order,self.order_checked,self._p))
            if verbose_init:
                print("ok")
        else:
            #self.order_checked = super(BrezingWeng,self).order()
            self.order_checked = super().order()
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
        #return "BrezingWeng_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (pairing-friendly curve k="+str(self._k)+") with seed "+str(self._u)+"\n"+super(BrezingWeng,self)._repr_()
        return "BrezingWeng_k"+str(self._k)+"_D"+str(self._D)+" p"+str(self._pbits)+" (pairing-friendly curve k="+str(self._k)+") with seed "+str(self._u)+"\n"+super()._repr_()

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

    def miller_loop_length(self):
        return self._u**coeff_mult_m(self._k,self._D)

    def print_parameters(self):
        print_parameters(self)

    def print_parameters_for_RELIC(self):
        print_parameters_for_RELIC(self)
