import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.arith.misc import XGCD, xgcd, gcd, valuation
from sage.arith.functions import lcm
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0

# it depends on the implementation of poly_cofactor_twist_g1_g2
allowed_k = [6,9,12,15,21,24,27,48] # TODO: 30,33,39,42,45

def polynomial_params(k):
    if (k % 3) != 0 or (k % 18 == 0):
        raise ValueError("Error with BLS k=0 mod 3 and k != 0 mod 18 but k={}={} mod 3 = {} mod 18 given".format(k, k%3, k%18))
    #return tnfs.curve.cyclo_kDe.polynomial_params(k,D=3,e0=1)
    # BLS
    D = 3
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    tx = QQx(x+1)
    Phi_k = QQx(cyclotomic_polynomial(k))
    rx = Phi_k
    # 3 | k -> sqrt(-3) does exit mod rx
    # if k=0 mod 6, x^(k//6) is a sixth root of unity (1+sqrt(-3))/2 -> 2*x^(k//6)-1 = sqrt(-3)
    # if k=3 mod 6, x^(k//3) is a cube root of unity (-1+sqrt(-3))/2 -> 2*x^(k//3)+1 = sqrt(-3)
    if k%6 == 0:
        yx = (tx-2)*(2*x**(k//6)-1)/3
    else:
        yx = (tx-2)*(2*x**(k//3)+1)/3
    if k != 6:
        yx = yx % rx
    px = (tx**2 + 3*yx**2)/4
    lambx = QQx(x**(k//3))
    px_denom = Integer(lcm([ci.denom() for ci in px.list()]))
    assert ((px+1-tx) % rx) == 0, "Error (px+1-tx)%rx != 0, px={}*({}), tx={}, rx={}, (px+1-tx)={}".format(px_denom,px_denom*px, tx,rx, (px+1-tx).factor())
    cx = (px+1-tx) // rx
    if 3**valuation(k,3) == k: # k is a power of 3
        rx = rx/3
        cx = cx*3
    assert (Phi_k(px) % rx) == 0, "Error k={}, Phi_k(px)%rx != 0, = {}".format(k, Phi_k(px) % rx)
    assert px.is_irreducible(), "Error px is not irreducible, px= {}".format(px.factor())
    assert rx.is_irreducible(), "Error rx is not irreducible, rx= {}".format(rx.factor())

    assert 4*px == (tx**2 + D*yx**2), "Error 4*px != tx^2+D*yx^2, {}*px = {}, tx={}, D={}, yx={}, tx^2+Dyx^2 = {}".format(px_denom, px_denom*px,tx,D,yx,tx**2+D*yx**2)
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
    if betax.leading_coefficient() > 0:
        betax = -betax-1
    #print("computed betax = {}".format(betax))
    
    assert ((betax**2 + betax + 1) % px) == 0
    assert ((lambx**2 + lambx + 1) % rx) == 0

    return px, rx, tx, cx, yx, betax, lambx, D

def coeffs_params(k):
    #return tnfs.curve.cyclo_kDe.coeff_params(k,D=3,e0=1)
    if (k % 3) != 0 or (k % 18 == 0):
        raise ValueError("Error with BLS k=0 mod 3 and k != 0 mod 18 but k={} = {} mod 3 = {} mod 18 given".format(k, k%3, k%18))
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
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
    BETAx_denom = Integer(lcm([ci.denom() for ci in betax.list()]))
    BETAx = [Integer(ci) for ci in (BETAx_denom*betax).list()]
    LAMBx_denom = Integer(lcm([ci.denom() for ci in lambx.list()]))
    LAMBx = [Integer(ci) for ci in (LAMBx_denom*lambx).list()]
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D

def poly_cofactor_gt(k):
    """Computes the co-factors for GT: Phi_k(p(x))/r(x)

    It is always irreducible:
    from tnfs.curve.bls import poly_cofactor_gt
    for k in [ki for ki in range(6, 49, 3) if (ki % 18) != 0]:
        print("\nBLS-{}".format(k))
        c = poly_cofactor_gt(k)
        print(["(deg {})^{} ".format(ci.degree(), ei) if ei>1 else "(deg {})".format(ci.degree()) for (ci, ei) in c.factor()])
        x = c.variables()[0]
        print(lcm([ci.denom() for ci in (c(3*x+1)).list()]))
    """
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    cx = cyclotomic_polynomial(k)(px)
    assert (cx % rx) == 0
    cx = cx // rx
    return cx

def poly_cofactor_twist_g1_g2(k):
    """Computes the curve co-factors for the twists"""
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    twx = px+1+tx
    m = 6
    u_m = [1, 4]
    # for specific needs: k=21 we want p=1 mod 7, k=17 we want p=1 mod9
    if k==6:
        # the cofactor of the order of E(Fp) has factors cxa and cxb
        cxa = (x-1)/3
        cxb = x-1
        assert (px+1-tx) == cxa*cxb*rx
        # the quadratic twist of E(Fp) has order twx = twxa*twxb
        # twx = (1/3) * (x^2 - 4*x + 7) * (x^2 + x + 1)
        twxa = x**2 - 4*x + 7
        twxb = (x**2 + x + 1)/3
        assert twx == twxa * twxb
        # the sextic twist over Fp
        E2_order = px+1-(-3*yx+tx)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        assert g2cx.is_irreducible()
        g2twx = px+1+(-3*yx+tx)/2
        assert g2twx.is_irreducible()
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, twxa, twxb, g2cx, g2twx]
        label_factors = ["cxa", "twxa", "twxb", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1, 1]

    elif k==9:
        # the cofactor of the order of E(Fp) has factors cxa and cxb
        # note that for k a power of 3, then the division by 3 is done on r, not c
        cxa = (x-1)
        cxb = (x-1)
        assert (px+1-tx) == cxa*cxb*rx
        # the quadratic twist of E(Fp) has order twx which is irreducible
        assert twx.is_irreducible()
        # E(Fp3) and its cubic twist
        px2 = px**2
        tx2 = tx**2 - 2*px
        px3 = px**3
        tx3 = tx*tx2 - px*tx # = tx**3 - 2*px*tx - px*tx = tx**3 - 3*px*tx
        yx3 = yx * (tx**2 - px)
        assert tx3**2 - 4*px3 == -D*yx3**2
        E2_order = px3+1-(-3*yx3-tx3)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx
        assert g2cx.is_irreducible()
        g2twx = px3+1+(-3*yx3-tx3)/2
        assert g2twx.is_irreducible()
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, twx, g2cx, g2twx]
        label_factors = ["cxa", "twx", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1]

    elif k==12:
        # the cofactor of the order of E(Fp) has factors cxa and cxb
        cxa = (x-1)/3
        cxb = (x-1)
        assert (px+1-tx) == cxa*cxb*rx
        # the quadratic twist of E(Fp) has order twx, irreducible
        assert twx.is_irreducible()
        # for G2, compute #E(Fp2) then compute its 6-th twist
        px2 = px**2
        tx2 = tx**2 - 2*px
        yx2 = yx*tx
        assert tx2**2 - 4*px2 == -D*yx2**2
        # now the 6-th twist that matches rx
        E2_order = px2+1-(-3*yx2+tx2)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        assert g2cx.is_irreducible()
        g2twx = px2+1+(-3*yx2+tx2)/2
        g2_twxa = (x**2 + x + 1)/3 # and is multiple of 3
        g2_twxb = (x**4 - 3*x**3 + 2*x**2 + 1)
        g2_twxc = (x**6 - 2*x**5 + 5*x**3 - 3*x**2 + x + 7)/3
        assert g2twx == g2_twxa * g2_twxb * g2_twxc
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, twx, g2cx, g2_twxa, g2_twxb, g2_twxc]
        label_factors = ["cxa", "twx", "g2cx", "g2twxa", "g2twxb", "g2twxc"]
        small_cofactors = [1, 1, 1, 1, 1, 1]

    elif k==15:
        # the cofactor of the order of E(Fp) has factors cxa and cxb
        cxa = (x-1)/3
        cxb = (x-1)
        cxc = (x**2 + x + 1)
        assert (px+1-tx) == cxa*cxb*cxc*rx
        # the quadratic twist of E(Fp) has order twx, irreducible
        assert twx.is_irreducible()
        # for G2, compute #E(Fp2) then compute its 6-th twist
        px5 = px**5
        #tx2 = tx**2 - 2*px
        #tx3 = tx*tx2 - px*tx
        #tx4 = tx*tx3 - px*tx2
        #tx5 = tx*tx4 - px*tx3
        tx5 = tx**5 - 5*px*tx**3 + 5*px**2*tx
        yx5 = yx * (tx**4 - 3*px*tx**2 + px**2)
        assert tx5**2 - 4*px5 == -D*yx5**2
        # now the 3-rd twist
        E2_order = px5+1+( 3*yx5+tx5)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # 3 factors of degrees 2, 10, 40
        g2cxa = (x**2 + x + 1)/3
        g2cxb = (x**10 - 3*x**9 + 3*x**8 - 3*x**6 + 4*x**5 - 6*x**4 + 6*x**3 - 6*x + 7)/3
        g2cxc = (x**40 - 7*x**39 + 21*x**38 - 36*x**37 + 42*x**36 - 39*x**35 + 24*x**34 + 9*x**33 - 42*x**32 + 56*x**31 - 56*x**30 + 39*x**29 - 12*x**28 - 3*x**27 + 21*x**26 - 32*x**25 + 8*x**24 - 6*x**23 + 27*x**22 - 12*x**20 - 9*x**19 - 15*x**18 - 3*x**17 + 44*x**16 + 16*x**15 - 12*x**14 - 51*x**13 - 33*x**12 + 6*x**11 + 103*x**10 + 26*x**9 + 3*x**8 - 108*x**7 - 42*x**6 + 15*x**5 + 66*x**4 + 24*x**3 + 9*x**2 + 8*x + 31)/81
        assert g2cx == 3 * g2cxa * g2cxb * g2cxc

        g2twx = px5+1-( 3*yx5+tx5)/2
        g2twxa = (x**12 - 2*x**11 + x**10 + x**7 + x**6 - 2*x**5 + x**2 + x + 1)/3
        g2twxb = (x**48 - 8*x**47 + 28*x**46 - 56*x**45 + 70*x**44 - 52*x**43 - 7*x**42 + 125*x**41 - 286*x**40 + 385*x**39 - 319*x**38 + 98*x**37 + 227*x**36 - 616*x**35 + 910*x**34 - 859*x**33 + 422*x**32 + 164*x**31 - 706*x**30 + 1165*x**29 - 1327*x**28 + 914*x**27 - 160*x**26 - 427*x**25 + 814*x**24 - 1189*x**23 + 1208*x**22 - 625*x**21 - 28*x**20 + 259*x**19 - 454*x**18 + 980*x**17 - 979*x**16 + 329*x**15 + 34*x**14 + 83*x**13 + 281*x**12 - 964*x**11 + 551*x**10 - 11*x**9 + 131*x**8 - 118*x**7 - 520*x**6 + 563*x**5 + 19*x**4 + 16*x**3 - 125*x**2 - 131*x + 271)/81
        assert g2twx == g2twxa * g2twxb
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, cxc, twx, g2cxa, g2cxb, g2cxc, g2twxa, g2twxb]
        label_factors = ["cxa", "cxc", "twx", "g2cxa", "g2cxb", "g2cxc", "g2twxa", "g2twxb"]
        small_cofactors = [1, 3, 1, 1, 1, 1, 1, 1]

    elif k==21:
        m = 6*7 # 42
        u_m = [1, 22]
        # the cofactor of the order of E(Fp) has factors cxa and cxb
        cxa = (x-1)/3
        cxb = (x-1)
        cxc = (x**2 + x + 1)
        assert (px+1-tx) == cxa*cxb*cxc*rx
        # the quadratic twist of E(Fp) has order twx, irreducible
        assert twx.is_irreducible()
        # for G2, compute #E(Fp2) then compute its 6-th twist
        px7 = px**7
        #tx2 = tx**2 - 2*px
        #tx3 = tx*tx2 - px*tx
        #tx4 = tx*tx3 - px*tx2
        #tx5 = tx*tx4 - px*tx3
        #tx6 = tx*tx5 - px*tx4
        #tx7 = tx*tx6 - px*tx5
        tx7 = tx**7 - 7*px*tx**5 + 14*px**2*tx**3 - 7*px**3*tx
        yx7 = yx * (tx**6 - 5*px*tx**4 + 6*px**2*tx**2 - px**3)
        assert tx7**2 - 4*px7 == -D*yx7**2
        # now the 3-rd twist
        E2_order = px7+1+( 3*yx7+tx7)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # 3 factors of degrees 2, 14, 84
        g2cxa = (x**2 + x + 1)/3
        g2cxb = (x**14 - 3*x**13 + 3*x**12 - 3*x**10 + 3*x**9 - 2*x**7 + 3*x**6 - 3*x**5 + 3*x**3 - 3*x**2 + 4)/3
        g2cxc = (x**84 - 11*x**83 + 55*x**82 - 166*x**81 + 341*x**80 - 517*x**79 + 627*x**78 - 655*x**77 + 569*x**76 - 212*x**75 - 624*x**74 + 1902*x**73 - 3190*x**72 + 4015*x**71 - 4274*x**70 + 3963*x**69 - 2660*x**68 - 179*x**67 + 4315*x**66 - 8470*x**65 + 11264*x**64 - 12274*x**63 + 11643*x**62 - 8868*x**61 + 3192*x**60 + 4704*x**59 - 12552*x**58 + 18108*x**57 - 20406*x**56 + 19605*x**55 - 15753*x**54 + 8427*x**53 + 1605*x**52 - 11589*x**51 + 18690*x**50 - 21771*x**49 + 21141*x**48 - 17418*x**47 + 10791*x**46 - 1674*x**45 - 7719*x**44 + 14025*x**43 - 16029*x**42 + 15303*x**41 - 12858*x**40 + 8100*x**39 - 1413*x**38 - 5151*x**37 + 9234*x**36 - 9414*x**35 + 7386*x**34 - 5928*x**33 + 3648*x**32 + 939*x**31 - 4836*x**30 + 6195*x**29 - 5196*x**28 + 2421*x**27 - 951*x**26 + 771*x**25 + 1581*x**24 - 4296*x**23 + 4659*x**22 - 3187*x**21 + 200*x**20 + 449*x**19 + 1105*x**18 + 823*x**17 - 3002*x**16 + 2082*x**15 - 1682*x**14 + 547*x**13 + 671*x**12 + 1509*x**11 - 720*x**10 - 1607*x**9 - 307*x**8 - 187*x**7 + 582*x**6 + 536*x**5 + 848*x**4 - 496*x**3 - 8*x**2 - 140*x + 547)/729
        assert g2cx == 3 * g2cxa * g2cxb * g2cxc

        g2twx = px7+1-( 3*yx7+tx7)/2
        g2twxa = (x**16 - 2*x**15 + x**14 + x**9 - 5*x**8 + 4*x**7 + x**2 - 2*x + 4)/3
        g2twxb = (x**96 - 12*x**95 + 66*x**94 - 220*x**93 + 495*x**92 - 792*x**91 + 924*x**90 - 786*x**89 + 426*x**88 + 143*x**87 - 1089*x**86 + 2463*x**85 - 3761*x**84 + 4158*x**83 - 3345*x**82 + 1764*x**81 + 192*x**80 - 2679*x**79 + 5631*x**78 - 7989*x**77 + 8442*x**76 - 6754*x**75 + 3735*x**74 - 81*x**73 - 3947*x**72 + 7731*x**71 - 10176*x**70 + 10542*x**69 - 8766*x**68 + 5193*x**67 - 596*x**66 - 3807*x**65 + 7017*x**64 - 8785*x**63 + 9072*x**62 - 7650*x**61 + 4653*x**60 - 1008*x**59 - 2097*x**58 + 4248*x**57 - 5517*x**56 + 5670*x**55 - 4773*x**54 + 3240*x**53 - 1368*x**52 - 510*x**51 + 1998*x**50 - 2655*x**49 + 2682*x**48 - 2331*x**47 + 1395*x**46 - 528*x**45 - 36*x**44 + 540*x**43 - 732*x**42 + 810*x**41 - 423*x**40 + 63*x**39 - 99*x**38 - 405*x**37 + 531*x**36 - 315*x**35 + 270*x**34 - 4*x**33 - 627*x**32 + 2067*x**31 - 1082*x**30 - 1062*x**29 + 450*x**28 + 300*x**27 - 483*x**26 - 1443*x**25 + 3001*x**24 + 297*x**23 - 960*x**22 - 148*x**21 - 432*x**20 - 1704*x**19 - 891*x**18 + 3831*x**17 + 1050*x**16 - 1662*x**15 + 213*x**14 - 774*x**13 - 2321*x**12 - 1476*x**11 + 2232*x**10 + 2357*x**9 + 315*x**8 - 435*x**7 - 111*x**6 - 756*x**5 - 711*x**4 + 167*x**3 + 504*x**2 + 417*x + 547)/729
        assert g2twx == g2twxa * g2twxb
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, cxc, twx, g2cxa, g2cxb, g2cxc, g2twxa, g2twxb]
        label_factors = ["cxa", "cxc", "twx", "g2cxa", "g2cxb", "g2cxc", "g2twxa", "g2twxb"]
        small_cofactors = [1, 3, 1, 1, 1, 1, 1, 1]

    elif k==24:
        cxa = (x-1)/3
        cxb = (x-1)
        assert (px+1-tx) == cxa*cxb*rx
        # twx = 1/3 * (x^2 + x + 1) * (x^8 - 3*x^7 + 3*x^6 - 4*x^4 + 6*x^3 - 3*x^2 - 3*x + 7)
        twx1 = (x**2 + x + 1)/3
        twx2 = (x**8 - 3*x**7 + 3*x**6 - 4*x**4 + 6*x**3 - 3*x**2 - 3*x + 7)
        assert twx == twx1 * twx2
        # for G2, compute #E(Fp4) then compute its 6-th twist
        px2 = px**2
        tx2 = tx**2 - 2*px
        px4 = px2**2
        tx4 = tx2**2 - 2*px2
        assert (px4 + 1 - tx4) == (px+1-tx)*(px+1+tx)*(px2+1+tx2)
        yx4 = yx*tx*(tx**2-2*px)
        assert tx4**2 - 4*px4 == -D*yx4**2
        # now the 6-th twist
        E2_order = px4+1-(-3*yx4+tx4)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        assert g2cx.is_irreducible()
        #print("g2cx={}\n=({})/81".format(g2cx,81*g2cx))
        g2twx = px4+1+(-3*yx4+tx4)/2
        #print("g2twx = {}".format(g2twx.factor()))
        g2_twxa = (x**10 - 2*x**9 + x**8 - x**6 - x**5 + 2*x**4 + x**2 + 4*x + 4) / 3
        g2_twxb = (x**10 - 2*x**9 + x**8 - x**6 + 5*x**5 - 4*x**4 + x**2 - 2*x + 4) / 3
        g2_twxc = (x**20 - 4*x**19 + 6*x**18 - 4*x**17 - x**16 + 8*x**15 - 12*x**14 + 8*x**13 + x**12 - 6*x**11 + 9*x**10 - 12*x**9 + 4*x**8 + 2*x**7 - 12*x**6 + 8*x**5 + 5*x**4 + 2*x**3 + 6*x**2 - 4*x + 4) / 9
        assert g2twx == g2_twxa * g2_twxb * g2_twxc
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, twx, g2cx, g2_twxa, g2_twxb, g2_twxc]
        label_factors = ["cxa", "twx", "g2cx", "g2twxa", "g2twxb", "g2twxc"]
        small_cofactors = [1, 1, 1, 1, 1, 1]

    elif k==27:
        m = 18
        u_m = [1, 10]
        # k is a power of 3, the denominator 3 is divided from r not c
        cxa = (x-1)
        cxb = (x-1)
        assert (px+1-tx) == cxa*cxb*rx
        assert twx.is_irreducible()
        # for G2: cubic twist over Fp9
        px2 = px**2
        px3 = px*px2
        px6 = px3**2
        px9 = px6*px3
        tx2 = tx**2 - 2*px
        tx3 = tx*tx2 - px*tx
        tx6 = tx3**2 - 2*px3
        tx9 = tx3*tx6 - px3*tx3
        yx9 = yx * (tx**2 - px) * (tx6 + px3)
        assert px9 == (tx9**2 + 3*yx9**2)/4
        # now the 3-th twist
        E2_order = px9+1+( 3*yx9+tx9)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        assert g2cx.is_irreducible()
        # quadratic twist of E2
        g2twx = px9+1-( 3*yx9+tx9)/2
        assert g2twx.is_irreducible()
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, twx, g2cx, g2twx]
        label_factors = ["cxa", "twx", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1]

    elif k == 48:
        #print("cx = {} = {}".format(cx, cx.factor()))
        #print("twx = {} = {}".format(twx, twx.factor()))
        cxa = (x-1)/3
        cxb = (x-1)
        assert (px+1-tx) == cxa*cxb*rx
        assert twx.is_irreducible()
        # for G2, sextic twist over Fp8
        px2 = px**2
        px4 = (px2)**2
        px8 = (px4)**2
        tx2 = tx**2 - 2*px
        tx4 = tx2**2 - 2*px2
        tx8 = tx4**2 - 2*px4
        assert (px8 + 1 - tx8) == (px+1-tx)*(px+1+tx)*(px2+1+tx2)*(px4+1+tx4)
        yx8 = yx * tx * tx2 * tx4
        assert tx8**2 - 4*px8 == -D*yx8**2
        # now the 6-th twist
        E2_order = px8+1-(-3*yx8+tx8)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        #print("g2cx={}\n={}".format(g2cx,g2cx.factor()))
        g2twx = px8+1+(-3*yx8+tx8)/2
        #print("g2twx = {}\n   = {}".format(g2twx, g2twx.factor()))
        assert g2cx.is_irreducible()
        g2_twxa = (x**2 + x + 1)
        g2_twxb = (x**16 - 3*x**15 + 3*x**14 - 3*x**12 + 3*x**11 - 3*x**9 + 2*x**8 + 1)
        g2_twxc = (x**18 - 2*x**17 + x**16 - x**10 + 5*x**9 - 4*x**8 + x**2 + x + 7)/9
        g2_twxd = (x**36 - 4*x**35 + 6*x**34 - 4*x**33 + x**32 - 2*x**28 + 8*x**27 - 12*x**26 + 8*x**25 - 2*x**24 + 3*x**20 - 6*x**19 + 9*x**18 - 12*x**17 + 6*x**16 - 2*x**12 + 2*x**11 + 6*x**10 + 8*x**9 - 14*x**8 + x**4 + 2*x**3 - 3*x**2 - 4*x + 13)/9
        g2_twxe = (x**72 - 8*x**71 + 28*x**70 - 56*x**69 + 70*x**68 - 56*x**67 + 28*x**66 - 8*x**65 - 3*x**64 + 32*x**63 - 112*x**62 + 224*x**61 - 280*x**60 + 224*x**59 - 112*x**58 + 32*x**57 + 6*x**56 - 68*x**55 + 208*x**54 - 380*x**53 + 460*x**52 - 380*x**51 + 208*x**50 - 68*x**49 - 6*x**48 + 92*x**47 - 232*x**46 + 356*x**45 - 400*x**44 + 356*x**43 - 232*x**42 + 92*x**41 + 3*x**40 - 80*x**39 + 154*x**38 - 200*x**37 + 205*x**36 - 164*x**35 + 100*x**34 - 44*x**33 - 6*x**32 + 44*x**31 - 52*x**30 + 68*x**29 - 16*x**28 - 112*x**27 + 56*x**26 + 80*x**25 - 42*x**24 - 8*x**23 + 10*x**22 - 20*x**21 - 38*x**20 + 196*x**19 - 152*x**18 - 116*x**17 + 114*x**16 - 4*x**15 - 4*x**14 + 8*x**13 - 28*x**12 - 172*x**11 + 104*x**10 + 140*x**9 - 39*x**8 + 4*x**7 + 10*x**6 + 16*x**5 + 37*x**4 + 52*x**3 - 44*x**2 - 68*x + 73)/81
        assert g2twx == g2_twxa * g2_twxb * g2_twxc * g2_twxd * g2_twxe
        cofactor_r = [1, 1]
        polys_cofact_twists = [cxa, twx, g2cx, g2_twxa, g2_twxb, g2_twxc, g2_twxd, g2_twxe]
        label_factors = ["cxa", "twx", "g2cx", "g2_twxa", "g2_twxb", "g2_twxc", "g2_twxd", "g2_twxe"]
        small_cofactors = [1, 1, 1, 3, 1, 1, 1, 1]

    else:
        raise ValueError("Error poly_cofactor_twist_g1_g2 in bls.py: k={} is not implemented, allowed_k={}".format(k, allowed_k))

    return m, u_m, cofactor_r, twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors


class BLS(EllipticCurve_finite_field):
    """
    A Barreto-Lynn-Scott curve of embedding degree k
    """
    def __init__(self, k, u, b=None, cofactor_r=1,verbose_init=False):
        """
        Create a BLS curve

        INPUT:
        - `k`: embedding degree
        - `u`: seed such that p=P_BLS_k(u) and r = R_BLS_k(u)
        - `b`: curve coefficient, y^2 = x^3 + b has a subgroup of order r
        - `cofactor_r`: a cofactor in r=P_BLS_k(r) (if not given and r not prime, raise error)

        u = 1 mod 3, so that P(u) is integer and not a multiple of 3,
        for k=12, u=1 mod 6 gives b=1
        (there are two possible choices for b: b, 1/b (a twist) but only one has order r)
        """
        self._k = k # embedding degree
        self._a = 0 # first curve parameter is 0 because j=0
        self._D = 3
        
        self._px, self._px_denom, self._rx, self._rx_denom, self._tx, self._tx_denom, self._cx, self._cx_denom, self._yx, self._yx_denom, self._betax, self._betax_denom, self._lambx, self._lambx_denom, _D = coeffs_params(k)
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

        if ((self._beta**2 + self._beta + 1) % self._p) != 0:
            raise ValueError("Error beta^2 + beta + 1 != 0 mod p")
        if ((self._lamb**2 + self._lamb + 1) % self._r) != 0:
            raise ValueError("Error lamb^2 + lamb + 1 != 0 mod r")

        self._Fpz = PolynomialRing(self._Fp, names=('z',))
        (self._z,) = self._Fpz._first_ngens(1)

        if b == None:
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
        # check the curve order.
        # With D=3 and j=0, there are 6 possible orders
        # p+1-t
        # p+1+t # quadratic twist
        # p+1-(3*y-t)/2  # cubic twist
        # p+1-(-3*y-t)/2 # other cubic twist
        # p+1-(3*y+t)/2  # sextic twist
        # p+1-(-3*y+t)/2 # other sextic twist

        self.curve_order = self._p + 1 - self._tr
        self.twist_order = self._p + 1 + self._tr
        cubic_twist1_order = self._p + 1 - (3*self._y-self._tr)//2
        cubic_twist2_order = self._p + 1 - (-3*self._y-self._tr)//2
        sextic_twist1_order = self._p + 1 - (3*self._y+self._tr)//2
        sextic_twist2_order = self._p + 1 - (-3*self._y+self._tr)//2
        for i in range(10):
            P = self.random_element()
            if self.curve_order*P != self(0):
                if self.twist_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a twist: (p+1+tr)*P = 0\ntr={}\nr={}\np+1+tr={}\n".format(self._tr,self.curve_order,self.twist_order))
                elif cubic_twist1_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a cubic twist: (p+1-(3*y-t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(3*y-t)/2)={}\n".format(self._tr,self._y,self.curve_order,cubic_twist1_order))
                elif cubic_twist2_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a cubic twist: (p+1-(-3*y-t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(-3*y-t)/2)={}\n".format(self._tr,self._y,self.curve_order,cubic_twist2_order))
                elif sextic_twist1_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a sextic twist: (p+1-(3*y+t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(3*y+t)/2)={}\n".format(self._tr,self._y,self.curve_order,sextic_twist1_order))
                elif sextic_twist2_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a sextic twist: (p+1-(-3*y+t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(-3*y+t)/2)={}\n".format(self._tr,self._y,self.curve_order,sextic_twist2_order))

        # computes a generator
        self._G = get_curve_generator_order_r_j0(self)
        self._Gx = self._G[0]
        self._Gy = self._G[1]

        # adjust beta and lamb according to the curve
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

    def _repr_(self):
        return "BLS{}-{} curve (Barreto-Lynn-Scott k={}) with seed {} ({:#x})\n".format(self._k, self._pbits, self._k, self._u, self._u)+super(BLS,self)._repr_()

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

    def poly_p(self):
        return self._px
    def poly_p_denom(self):
        return self._px_denom
    def poly_r(self):
        return self._rx
    def poly_r_denom(self):
        return self._rx_denom

    def miller_loop_length():
        return self._u
    
    def print_parameters(self):
        tnfs.curve.pairing_friendly_curve.print_parameters(self)
        
    def print_parameters_for_RELIC(self):
        tnfs.curve.pairing_friendly_curve.print_parameters_for_RELIC(self)
