"""
KSS pairing-friendly curves from the paper
Constructing Brezing-Weng pairing friendly elliptic curves using elements in the cyclotomic field
Ezekiel J. Kachisa and Edward F. Schaefer and Michael Scott
https://eprint.iacr.org/2007/452

"""

import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from sage.functions.log import log
from sage.functions.other import ceil
from sage.arith.misc import XGCD, xgcd
from sage.arith.functions import lcm
from sage.rings.integer import Integer
from sage.rings.rational_field import Q, QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_b_j0, get_curve_generator_order_r_j0
from tnfs.curve.pairing_friendly_curve import get_curve_parameter_a_j1728, get_curve_generator_order_r_j1728

allowed_k = [8,16,18,32,36,40,54]
possible_k = [8,16,18,32,36,40,54]

def polynomial_params(k):
    if not k in possible_k:
        raise ValueError("Error with KSS k={} does not exist".format(k))
    if not k in allowed_k:
        raise ValueError("Error with KSS k={} is not implemented".format(k))

    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    if k == 8:
        px=(x**6+2*x**5-3*x**4+8*x**3-15*x**2-82*x+125)/180
        rx= (x**4-8*x**2+25)/(18*25)
        tx = (2*x**3-11*x+15)/15
        cx = 5*(x**2 + 2*x + 5)/2
        yx = (x**3 + 5*x**2 + 2*x - 20)/15
        betax = (3*x**5 + 7*x**4 + 8*x**3 + 34*x**2 - 41*x - 113)/66
        lambx = (x**2 - 4)/3
        D=1
    elif k == 16:
        px = (x**10 + 2*x**9 + 5*x**8 + 48*x**6 + 152*x**5 + 240*x**4 + 625*x**2 + 2398*x + 3125)/980
        rx = (x**8 + 48*x**4 + 625)/61250 # 625 = 5^4, 61250 = 2*5^4*7^2
        tx = (2*x**5 + 41*x + 35)/35
        cx = 125 * (x**2 + 2*x + 5)/2 # C such that P+1-T = C*R
        yx = (x**5 + 5*x**4 + 38*x + 120)/70 # Y such that T^2 - 4*P = -4*Y^2
        betax = (x**9-11*x**8-16*x**7-120*x**6-32*x**5-498*x**4-748*x**3-2740*x**2-3115*x-5651)/4018
        lambx = (x**4 + 24)/7 # sqrt(-1) mod R
        D=4
    elif k == 18:
        px = (x**8 + 5*x**7 + 7*x**6 + 37*x**5 + 188*x**4 + 259*x**3 + 343*x**2 + 1763*x + 2401)/21
        rx = (x**6 + 37*x**3 + 343)/343 # 343 = 7^3
        cx = (x**2 + 5*x + 7)*49/3
        tx = (x**4 + 16*x + 7)/7
        yx = (5*x**4 + 14*x**3 + 94*x + 259)/21 # Y such that T^2 - 4*P = -3*Y^2
        betax = (x**7 + 3*x**6 + 4*x**5 + 44*x**4 + 118*x**3 + 71*x**2 + 483*x + 1118)/24
        lambx = x**3 + 18
        D=3
    elif k == 32:
        px = (x**18-6*x**17+13*x**16+57120*x**10-344632*x**9 + 742560*x**8 + 815730721*x**2 - 4948305594*x + 10604499373)/2970292
        rx = (x**16+57120*x**8+815730721)/93190709028482 # den = 2 * 13^8 * 239^2
        tx = (-2*x**9 - 56403*x + 3107)/3107
        yx = (3*x**9 - 13*x**8 + 86158*x - 371280)/(13*239)
        cx = 13**7 * (x**2 - 6*x + 13)/2
        betax = (x**17 + 469*x**16 - 2824*x**15 + 12272*x**14 - 36712*x**13 + 159536*x**12 - 477256*x**11 + 2073968*x**10 - 6147208*x**9 + 26788318*x**8 - 81107540*x**7 + 350475892*x**6 - 1054398020*x**5 + 4556186596*x**4 - 13707174260*x**3 + 59230425748*x**2 - 177377534659*x + 382523614565)/6443591526
        lambx = (x**8 + 28560)/239
        D = 1
    elif k == 36:
        px = (x**14-4*x**13+7*x**12+683*x**8-2510*x**7 + 4781*x**6 + 117649*x**2 - 386569*x + 823543)/28749
        rx = (x**12 + 683*x**6 + 117649)/161061481
        cx = 16807 * (x**2 - 4*x + 7)/3
        #rx = (x**12 + 683*x**6 + 117649)/9583
        #cx = (x**2 - 4*x + 7)/3
        tx = (2*x**7+757*x+259)/259
        yx = (4*x**7 - 14*x**6 + 1255*x - 4781)/(3*7*37)
        betax = (x**13 - 43*x**12 + 170*x**11 - 574*x**10 + 1190*x**9 - 4018*x**8 + 9013*x**7 - 29264*x**6 + 53726*x**5 - 195244*x**4 + 376082*x**3 - 1366708*x**2 + 2750223*x - 5494929)/1036333
        lambx = (x**6 + 323)/37
        D = 3
    elif k == 40:
        px = (x**22 - 2*x**21 + 5*x**20 + 6232*x**12 - 10568*x**11 + 31160*x**10 + 9765625*x**2 - 13398638*x + 48828125)/1123380
        rx = (x**16 + 8*x**14 + 39*x**12 + 112*x**10 - 79*x**8 + 2800*x**6 + 24375*x**4 + 125000*x**2 + 390625)/2437890625
        cx = 78125 * (x**2 - 2*x + 5) * (x**4 - 8*x**2 + 25) / 36
        #rx = (x**16 + 8*x**14 + 39*x**12 + 112*x**10 - 79*x**8 + 2800*x**6 + 24375*x**4 + 125000*x**2 + 390625)/31205
        #cx = (x**2 - 2*x + 5) * (x**4 - 8*x**2 + 25) / 36
        tx = (2*x**11 + 6469*x + 1185)/1185
        yx = (x**11 - 5*x**10 + 2642*x - 15580)/(3*5*79)
        betax = (3*x**21 - 167*x**20 + 352*x**19 - 1640*x**18 + 1760*x**17 - 8200*x**16 + 8800*x**15 - 41000*x**14 + 44000*x**13 - 205000*x**12 + 238696*x**11 - 1037954*x**10 + 944204*x**9 - 5096020*x**8 + 4721020*x**7 - 25480100*x**6 + 23605100*x**5 - 127400500*x**4 + 118025500*x**3 - 637002500*x**2 + 619424375*x - 1612604207)/ (2 * 3 * 79**2 * 6469)
        lambx = (x**10 + 3116)/237
        D = 1
    elif k == 54:
        px = 1 + 3*x + 3*x**2 + 3**5 * x**9 + 3**5 * x**10 + 3**6 * x**10 + 3**6 * x**11 + 3**9 * x**18 + 3**10 * x**19 + 3**10 * x**20
        rx = 1 + 3**5 * x**9 + 3**9 * x**18
        cx = 1 + 3*x + 3*x**2
        tx = 1 + 3**5 * x**10
        yx = 3**5 * x**10 + 2*3**4 * x**9 + 2*x + 1
        betax = (3**10*x**19 + 3**10*x**18 + 2*3**9*x**17 + 3**9*x**16 + 2*3**8*x**15 + 3**8*x**14 + 2*3**7*x**13 + 3**7*x**12 + 2*3**6*x**11 + 2*3**6*x**10 + 2**2*3**5*x**9 + 3**5*x**8 + 2*3**4*x**7 + 3**4*x**6 + 2*3**3*x**5 + 3**3*x**4 + 2*3**2*x**3 + 3**2*x**2 + 3**2*x + 2)/2
        lambx = 243*x**9 + 1
        D=3

    assert (px+1-tx) == cx*rx
    assert (tx**2 + D*yx**2) == 4*px
    if D == 1 or D == 4:
        assert ((betax**2 + 1) % px) == 0
        assert ((lambx**2 + 1) % rx) == 0
    elif D == 3:
        assert ((betax**2 + betax + 1) % px) == 0
        assert ((lambx**2 + lambx + 1) % rx) == 0
    else:
        assert ((betax**2 + D) % px) == 0
        assert ((lambx**2 + D) % rx) == 0

    return px, rx, tx, cx, yx, betax, lambx, D

def coeffs_params(k):
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
    if betax != 0 and lambx != 0:
        BETAx_denom = Integer(lcm([ci.denom() for ci in betax.list()]))
        BETAx = [Integer(ci) for ci in (BETAx_denom*betax).list()]
        LAMBx_denom = Integer(lcm([ci.denom() for ci in lambx.list()]))
        LAMBx = [Integer(ci) for ci in (LAMBx_denom*lambx).list()]
    else:
        BETAx_denom=0; BETAx=0; LAMBx_denom=0; LAMBx=0
    return Px, Px_denom, Rx, Rx_denom, Tx, Tx_denom, Cx, Cx_denom, Yx, Yx_denom, BETAx, BETAx_denom, LAMBx, LAMBx_denom, D

def poly_cofactor_gt(k):
    """Computes the co-factors for GT: Phi_k(p(x))/r(x)
    It is always irreducible:
    from tnfs.curve.kss import poly_cofactor_gt
    for k in [8,16,18,32,36,40,54]:
        print("\nKSS-{}".format(k))
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

def poly_cofactor_twist_g1_g2(k):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k)
    twx = px+1+tx
    if k == 8:
        # for G2: quartic twist over Fp2
        tx2 = tx**2 - 2*px
        px2 = px**2
        assert (px2 + 1 - tx2) == (px+1-tx)*(px+1+tx)
        yx2 = tx*yx
        assert tx2**2 - 4*px2 == -D*yx2**2
        # now the 4-th twist that matches rx
        if D == 1:
            E2_order = px2 + 1 + yx2
            g2twx = px2 + 1 - yx2
        elif D == 4:
            E2_order = px2 + 1 + 2*yx2
            g2twx = px2 + 1 - 2*yx2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        m = 30
        #u_m = [5,7,17,19,25,29]
        #cofactor_r = [450, 18, 18, 18, 450, 18] # cofactor of r according to u_m
        #cofactor_r = [25, 1, 1, 1, 25, 1] # cofactor of r according to u_m
        u_m = [5, 25]
        cofactor_r = [1, 1]
        polys_cofact_twists = [cx, twx, g2cx, g2twx]
        label_factors = ["cx", "twx", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1]

    elif k == 16:
        # for G2: quartic twist over Fp4
        px2 = px**2
        px4 = (px2)**2
        tx2 = tx**2 - 2*px
        tx4 = tx2**2 - 2*px2
        assert (px4 + 1 - tx4) == (px+1-tx)*(px+1+tx)*(px2+1+tx2)
        yx4 = yx*tx*(tx**2-2*px)
        assert tx4**2 - 4*px4 == -D*yx4**2
        # now the 4-th twist that matches rx
        if D == 4:
            E2_order = px4 + 1 + 2*yx4
            g2twx = px4 + 1 - 2*yx4
        elif D ==1:
            E2_order = px4 + 1 + yx4
            g2twx = px4 + 1 - yx4
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        m = 70 # we should have u = 25,45 mod 70 <=> +/- 25 mod 70
        u_m = [25, 45] # 25 or 45
        cofactor_r = [1, 1]
        polys_cofact_twists = [cx, twx, g2cx, g2twx]
        label_factors = ["cx", "twx", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1]

    elif k == 18:
        # for G2: sextic twist over Fp3
        px2 = px**2
        tx2 = tx**2 - 2*px
        tx3 = tx*tx2 - px*tx
        px3 = px*px2
        yx3 = yx * (px - tx**2)
        assert (tx3**2 - 4*px3) == -D*yx3**2
        # now the 6-th twist that matches rx
        E2_order = px3+1-( 3*yx3+tx3)/2
        g2twx = px3+1+( 3*yx3+tx3)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        # for px(u) to be an integer we need u = 5,7,14,19 mod 21 but for u_m = 7, 19 then 3 | px(21*x+u_m)
        # for rx(u) to be an integer we need u = 0,7,14 mod 21 (u=0 mod 7) and 3 | rx(21*x + 7)
        m = 21
        u_m = [14] # 14 mod 21 only. With 7 mod 21, px is always multiple of 3.
        cofactor_r = [1]
        polys_cofact_twists = [cx, twx, g2cx, g2twx]
        label_factors = ["cx", "twx", "g2cx", "g2twx"]
        small_cofactors = [343, 3, 3, 1]

    elif k == 32:
        # for G2: quartic twist over Fp8
        px2 = px**2
        px4 = (px2)**2
        px8 = (px4)**2
        tx2 = tx**2 - 2*px
        tx4 = tx2**2 - 2*px2
        tx8 = tx4**2 - 2*px4
        assert (px8 + 1 - tx8) == (px+1-tx)*(px+1+tx)*(px2+1+tx2)*(px4+1+tx4)
        #yx4 = yx*tx*(tx**2-2*px)
        #yx8 = yx * tx * (tx**2 - 2*px) * (tx**4 - 4*px*tx**2 + 2*px**2)
        yx8 = yx * tx * tx2 * tx4
        assert tx8**2 - 4*px8 == -D*yx8**2
        # now the 4-th twist that matches rx
        if D == 4:
            E2_order = px8 + 1 - 2*yx8
            g2twx = px8 + 1 + 2*yx8
        elif D == 1:
            E2_order = px8 + 1 - yx8
            g2twx = px8 + 1 + yx8
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        m = 6214
        u_m = [325, 5889] # 325 or 5889 = -325 mod 6214
        cofactor_r = [1, 1]
        polys_cofact_twists = [cx, twx, g2cx, g2twx]
        label_factors = ["cx", "twx", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1]

    elif k == 36:
        # for G2: sextic twist over Fp6
        px2 = px**2
        px4 = px2**2
        px6 = px4*px2
        tx2 = tx**2 - 2*px
        tx4 = tx2**2 - 2*px2
        tx6 = tx2*tx4 - px2*tx2
        assert (px6 + 1 - tx6) == (px+1-tx)*(px+1+tx)*(px2+1+tx2 + px + tx*(px+1))*(px2+1+tx2 + px - tx*(px+1))
        yx6 = tx * yx * (tx**2 - 3*px) * (tx**2 - px)
        assert (tx6**2 - 4*px6) == -D*yx6**2
        # now the 6-th twist that matches rx
        E2_order = px6 + 1 - (-3*yx6 + tx6)/2
        g2twx = px6 + 1 +(-3*yx6+tx6)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible
        m = 777
        u_m = [287, 308, 497, 539, 728, 749]
        # the denominator of tx is 7*37, the denominator of yx is 3*7*37, the denominator of px is 3*7*37^2
        # indeed for certain values of x mod 3*7*37^2, then px is not irreducible: x = 2093,
        # but it simplify a lot the description to consider the congruence equation modulo 3*7*37 = 777 only.
        # thanks to Jean Gasnier for discussions about it.
        cofactor_r = [1, 1, 1, 1, 1, 1]
        polys_cofact_twists = [cx, twx, g2cx, g2twx]
        label_factors = ["cx", "twx", "g2cx", "g2twx"]
        small_cofactors = [1]*len(u_m)

    elif k == 40:
        # cx = (78125/36) * (x^2 - 2*x + 5) * (x^4 - 8*x^2 + 25)
        cxa = (x**2 - 2*x + 5)/4
        cxb = (x**4 - 8*x**2 + 25)/18
        cx0 = 78125*2
        assert cx == cxa*cxb*cx0
        # for G2: quartic twist over Fp10
        px2 = px**2
        px4 = (px2)**2
        px5 = px4*px
        px10 = px5**2
        tx2 = tx**2 - 2*px
        tx3 = tx*tx2 - px*tx
        tx4 = tx*tx3 - px*tx2
        tx5 = tx*tx4 - px*tx3
        tx10 = tx5**2 - 2*px5
        assert (px10+1-tx10) == (px5+1-tx5) * (px5+1+tx5)
        #t^2 * (t^2 - 4*p) * (t^4 - 5*p*t^2 + 5*p^2)^2 * (t^4 - 3*p*t^2 + p^2)^2
        yx10 = yx * tx5 * (tx*tx3 + px2)
        assert tx10**2 - 4*px10 == -D*yx10**2
        # now the 4-th twist that matches rx
        if D == 4:
            E2_order = px10 + 1 + 2*yx10
            g2twx = px10 + 1 - 2*yx10
        elif D == 1:
            E2_order = px10 + 1 + yx10
            g2twx = px10 + 1 - yx10
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx
        # there are 3 parts
        g2cx0 = 2
        g2cxa = (x**4 - 8*x**2 + 25)/(25*18) # this is the same as cxb
        g2cxb = (x**40 - 4*x**39 + 22*x**38 - 52*x**37 + 176*x**36 - 316*x**35 + 858*x**34 - 1228*x**33 + 2464*x**32 - 1924*x**31 + 10726*x**30 - 30756*x**29 + 191120*x**28 - 428268*x**27 + 1572410*x**26 - 2657244*x**25 + 7801280*x**24 - 10551252*x**23 + 23099990*x**22 - 17978916*x**21 + 48136994*x**20 - 77629356*x**19 + 558361526*x**18 - 1168445628*x**17 + 4722694208*x**16 - 7406831124*x**15 + 23822515514*x**14 - 30043508292*x**13 + 72512768912*x**12 - 55177288236*x**11 + 106258013446*x**10 - 63737472620*x**9 + 548188465648*x**8 - 1053509178500*x**7 + 4766701317834*x**6 - 6834636612500*x**5 + 24428898901472*x**4 - 28339363437500*x**3 + 76263658265926*x**2 - 55848992187500*x + 94754225231233)/5608811664
        # arange for denominators
        g2cxc = g2cx // (g2cx0*g2cxa*g2cxb)

        g2twx0 = 2
        g2twxa = (x**44 - 4*x**43 + 14*x**42 - 20*x**41 + 25*x**40 + 12464*x**34 - 46064*x**33 + 166912*x**32 - 230320*x**31 + 311600*x**30 + 58369074*x**24 - 197579328*x**23 + 747170508*x**22 - 978909600*x**21 + 1459226850*x**20 + 121718750000*x**14 - 373406874032*x**13 + 1489818644656*x**12 - 1811027136880*x**11 + 3048293571200*x**10 + 95367431640625*x**4 - 261692148437500*x**3 + 1117837978524302*x**2 - 1220696679687500*x + 2402039916499225)/2523965248800
        g2twxb = g2twx // (g2twx0*g2twxa)

        m = 2370
        u_m = [415, 1165, 1205, 1955] # 1955 = 2370-415, 1205 = 2370-1165
        cofactor_r = [1, 1, 1, 1]
        polys_cofact_twists = [cx, twx, g2cxa, g2cxb, g2cxc, g2twxa, g2twxb]
        label_factors = ["cx", "twx", "g2cxa", "g2cxb", "g2cxc", "g2twxa", "g2twxb"]
        small_cofactors = [1, 1, 1, 1, 1, 1, 1]

    elif k == 54:
        # for G2: sextic twist over Fp9
        px2 = px**2
        px3 = px*px2
        px6 = px3**2
        px9 = px6*px3
        tx2 = tx**2 - 2*px
        tx3 = tx*tx2 - px*tx
        tx6 = tx3**2 - 2*px3
        tx9 = tx3*tx6 - px3*tx3
        yx9 = yx * (tx**2 - px) * (tx6 + px3)

        # now the 6-th twist that matches rx
        E2_order = px9+1-( 3*yx9+tx9)/2
        g2twx = px9+1+( 3*yx9+tx9)/2
        assert (E2_order % rx) == 0
        g2cx = E2_order // rx # irreducible

        m = 1
        u_m = [0]
        cofactor_r = [1]
        polys_cofact_twists = [cx, twx, g2cx, g2twx]
        label_factors = ["cx", "twx", "g2cx", "g2twx"]
        small_cofactors = [1, 1, 1, 1]

    return m, u_m, cofactor_r, twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors

def congruence_constraints(k):
    if k == 8:
        m = 30
        u_m = [5, 25]
    elif k == 16:
        m = 70
        u_m = [25, 45]
    elif k == 18:
        m = 21
        u_m = [14] # 14 mod 21 only. With 7 mod 21, px is always multiple of 3.
    elif k == 32:
        m = 6214
        u_m = [325, 5889] # 325 or 5889 = -325 mod 6214
    elif k == 36:
        m = 777
        u_m = [287, 308, 497, 539, 728, 749] # it means 2 mod 3, 0 mod 7, and one of 9, 12, 16, 21, 25, 28 mod 37
    elif k == 40:
        m = 2370
        u_m = [415, 1165, 1205, 1955] # 1955 = 2370-415, 1205 = 2370-1165
    elif k == 54:
        m = 1
        u_m = [0]
    return m, u_m

class KSS(EllipticCurve_finite_field):
    """
    A Kachisa-Schaefer-Scott curve of embedding degree k in (8,16,18,32,36,40,54)
    """
    def __init__(self, k, u, a=None, b=None, cofactor_r=None, verbose_init=False):
        """
        Create a KSS curve

        INPUT:
        - `k`: embedding degree
        - `u`: seed such that p=P_KSS_k(u) and r = R_KSS_k(u)
        - `a`: curve coefficient, if D=1: y^2 = x^3 + a*x has a subgroup of order r
        - `b`: curve coefficient, if D=3: y^2 = x^3 + b has a subgroup of order r
        - `cofactor_r`: a cofactor in r=P_KSS_k(r) (if not given and r not prime, raise error)

        """
        self._k = k # embedding degree
        self._px, self._px_denom, self._rx, self._rx_denom, self._tx, self._tx_denom, self._cx, self._cx_denom, self._yx, self._yx_denom, self._betax, self._betax_denom, self._lambx, self._lambx_denom, self._D = coeffs_params(k)

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
        # if D=3:
        # psi: (x,y) -> (beta*x,y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+x+1 mod p, beta = (-1 + sqrt(Fp(-3)))/2
        # and lambda is a root of x^2+x+1 mod r, lamb = (-1 + sqrt(Fr(-3)))/2
        # there are two choices: beta, -beta-1, and the same for lamb: lamb, -lamb-1
        # if D=1:
        # psi: (x,y) -> (-x,i*y) = [lambda mod r]*(x,y) if (x,y) is of order r
        # where beta is a root of x^2+1 mod p, beta = sqrt(-1) mod p
        # and lambda is a root of x^2+1 mod r, lamb = sqrt(-1) mod r
        # there are two choices: beta, -beta, and the same for lamb: lamb, -lamb
        # arbitrarily the positive ones are chosen in the polynomials, adjust the choice
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
        if cofactor_r is not None:
            self._cofactor_r = Integer(cofactor_r)
        else:
            self._cofactor_r = Integer(1)
        if self._cofactor_r > 1:
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

        if self._D == 3:
            self._a = 0
            self._ap = self._Fp(0) # first curve parameter is 0 because j=0
            if b == None:
                # check that beta = 2*U/(-3*V-U) before, where U=t/2, V = y/2 and 2V = 2 mod 3
                self._b, self._bp = get_curve_parameter_b_j0(self._tr, self._y, self._p, self._Fp)
            else:
                try:
                    b = Integer(b)
                except:
                    raise
                self._b = b
                self._bp = self._Fp(b)
        elif self._D == 1 or self._D == 4:
            self._b = 0
            self._bp = self._Fp(0) # second curve parameter is 0 because j=1728
            if a == None:
                # check that beta = U/V, where U=t/2, V = y
                self._a, self._ap = get_curve_parameter_a_j1728(self._tr, self._y, self._p, self._Fp)
            else:
                try:
                    a = Integer(a)
                except:
                    raise
                self._a = a
                self._ap = self._Fp(a)
        elif a == None and b == None:
            raise ValueError("Only D=1 or D=3 is implemented for guessing a,b, please provide a,b")
        else:
            try:
                a = Integer(a)
                b = Integer(b)
            except:
                raise
            self._a = a
            self._ap = self._Fp(a)
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
        # check the curve order.
        # With D=3 and j=0, there are 6 possible orders
        # p+1-t
        # p+1+t # quadratic twist
        # p+1-(3*y-t)/2  # cubic twist
        # p+1-(-3*y-t)/2 # other cubic twist
        # p+1-(3*y+t)/2  # sextic twist
        # p+1-(-3*y+t)/2 # other sextic twist
        # With D=1 and j=1728, there are 4 possible orders
        # p+1-t
        # p+1+t # quadratic twist
        # p+1-y # quartic twist
        # p+1+y # other quartic twist
        # with the convention D=4, this is p+1-2*y and p+1+2*y

        self.curve_order = self._p + 1 - self._tr
        self.twist_order = self._p + 1 + self._tr
        if self._D == 3:
            cubic_twist1_order = self._p + 1 - (3*self._y-self._tr)//2
            cubic_twist2_order = self._p + 1 - (-3*self._y-self._tr)//2
            sextic_twist1_order = self._p + 1 - (3*self._y+self._tr)//2
            sextic_twist2_order = self._p + 1 - (-3*self._y+self._tr)//2
        elif self._D == 1:
            quartic_twist1_order = self._p + 1 - self._y
            quartic_twist2_order = self._p + 1 + self._y
        elif self._D == 4:
            quartic_twist1_order = self._p + 1 - 2*self._y
            quartic_twist2_order = self._p + 1 + 2*self._y
        for i in range(10):
            P = self.random_element()
            if self.curve_order*P != self(0):
                if self.twist_order*P == self(0):
                    raise ValueError("Wrong curve order: this one is a twist: (p+1+tr)*P = 0\ntr={}\nr={}\np+1+tr={}\n".format(self._tr,self.curve_order,self.twist_order))
                elif self._D == 3:
                    if cubic_twist1_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a cubic twist: (p+1-(3*y-t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(3*y-t)/2)={}\n".format(self._tr,self._y,self.curve_order,cubic_twist1_order))
                    elif cubic_twist2_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a cubic twist: (p+1-(-3*y-t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(-3*y-t)/2)={}\n".format(self._tr,self._y,self.curve_order,cubic_twist2_order))
                    elif sextic_twist1_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a sextic twist: (p+1-(3*y+t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(3*y+t)/2)={}\n".format(self._tr,self._y,self.curve_order,sextic_twist1_order))
                    elif sextic_twist2_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a sextic twist: (p+1-(-3*y+t)/2)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-(-3*y+t)/2)={}\n".format(self._tr,self._y,self.curve_order,sextic_twist2_order))
                elif self._D == 1:
                    if quartic_twist1_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a quartic twist: (p+1-y)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-y)={}\n".format(self._tr,self._y,self.curve_order,quartic_twist1_order))
                    elif quartic_twist2_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a quartic twist: (p+1+y)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1+y)={}\n".format(self._tr,self._y,self.curve_order,quartic_twist2_order))
                elif self._D == 4:
                    if quartic_twist1_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a quartic twist: (p+1-2*y)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1-2*y)={}\n".format(self._tr,self._y,self.curve_order,quartic_twist1_order))
                    elif quartic_twist2_order*P == self(0):
                        raise ValueError("Wrong curve order: this one is a quartic twist: (p+1+2*y)*P = 0\ntr={}\ny={}\np+1-t={}\n(p+1+2*y)={}\n".format(self._tr,self._y,self.curve_order,quartic_twist2_order))

        # computes a generator
        if self._D == 3:
            self._G = get_curve_generator_order_r_j0(self)
        elif self._D == 1 or self._D == 4:
            self._G = get_curve_generator_order_r_j1728(self)
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
        elif self._D == 1 or self._D == 4:
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
        return "KSS{}-{} curve (Kachisa-Schaefer-Scott k={}) with seed {} ({:#x})\n".format(self._k, self._pbits, self._k, self._u, self._u)+super(KSS,self)._repr_()

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
        return self._u

    def print_parameters(self):
        tnfs.curve.pairing_friendly_curve.print_parameters(self)

    def print_parameters_for_RELIC(self):
        tnfs.curve.pairing_friendly_curve.print_parameters_for_RELIC(self)
