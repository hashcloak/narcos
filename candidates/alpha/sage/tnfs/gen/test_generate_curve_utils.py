from sage.all_cmdline import *   # import sage library

from sage.misc.functional import cyclotomic_polynomial
from sage.rings.integer_ring import ZZ

import tnfs
import tnfs.curve.mnt_k
import tnfs.curve.mntg
import tnfs.curve.galbraith_mckee_valenca
import tnfs.curve.kss
import tnfs.curve.bls
import tnfs.curve.bn
from tnfs.gen.generate_curve_utils import find_seed_congruence_for_adicity_QQpoly, find_seed_congruence_for_adicity_ZZpoly

# sage -python -m tnfs.gen.test_generate_curve_utils

def test_find_seed_congruence_for_adicity_vx_ox(vx, v_str, other_poly_str, u_mod_m, m, val_l, l):
    print(v_str)
    x = vx.parent().gen(0)
    dict_mj_uj = find_seed_congruence_for_adicity_QQpoly(vx, u_mod_m, m, val_l, l)
    print("dict_mj_uj = {}".format(dict_mj_uj))
    for mj in dict_mj_uj:
        for uj in dict_mj_uj[mj]:
            v_uj = ZZ(vx(uj))
            print("uj={}={:#x} mod mj={}={:#x}={} v(uj)={:#x} {} bits {}-valuation = {}, expected {}".format(uj, uj, mj, mj, mj.factor(), v_uj, v_uj.nbits(), l, v_uj.valuation(l), val_l))
            for other_poly, other_str in other_poly_str:
                ox = other_poly(mj*x + uj)
                cont = gcd([ZZ(oi) for oi in ox.list()])
                print("    content({}(mj*x+uj)) = {}, {}-valuation {}".format(other_str, cont, l, cont.valuation(l)))
    
def test_find_seed_congruence_for_adicity_ZZ_vx_ox(vx, v_str, other_poly_str, val_l, l):
    print(v_str)
    x = vx.parent().gen(0)
    dict_mj_uj = find_seed_congruence_for_adicity_ZZpoly(vx, val_l, l)
    print("dict_mj_uj = {}".format(dict_mj_uj))
    for mj in dict_mj_uj:
        for uj in dict_mj_uj[mj]:
            v_uj = ZZ(vx(uj))
            #v_ujx = vx(x*mj + uj)
            print("uj={}={:#x} mod mj={}={:#x}={} v(uj)={:#x} {} bits {}-valuation = {}, expected {}".format(uj, uj, mj, mj, mj.factor(), v_uj, v_uj.nbits(), l, v_uj.valuation(l), val_l))
            for other_poly, other_str in other_poly_str:
                ox = other_poly(mj*x + uj)
                cont = gcd([ZZ(oi) for oi in ox.list()])
                print("    content({}(mj*x+uj)) = {}, {}-valuation {}".format(other_str, cont, l, cont.valuation(l)))

def test_find_seed_congruence_for_adicity_QQpoly():
    ## KSS polynomial

    valuations = [64,64,64,64,64,64]
    i = 0
    for k in [8,16,18,32,40,54]: #36
        px,rx,tx,cx,yx,betax,lambx, D = tnfs.curve.kss.polynomial_params(k)
        m, u_mod_m = tnfs.curve.kss.congruence_constraints(k)
        val_l = valuations[i]
        l = 2
        x = px.parent().gen(0)
        print("\nKSS-{} rx = {}".format(k, rx.factor()))
        for vx, v_str, other_poly_str in [(px-1, "px-1", [(rx-1, "(rx-1)"), (px, "px"), (rx, "rx")]), (px(-x)-1, "p(-x)-1", [(rx(-x)-1, "(r(-x)-1)"), (px(-x), "p(-x)"), (rx(-x), "r(-x)")]), (rx-1, "rx-1", [(px-1, "(px-1)"), (rx, "rx"), (px, "px")]), (rx(-x)-1, "r(-x)-1", [(px(-x)-1, "(p(-x)-1)"), (rx(-x), "r(-x)"), (px(-x), "p(-x)")])]:
            test_find_seed_congruence_for_adicity_vx_ox(vx, v_str, other_poly_str, u_mod_m, m, val_l, l)
        i = i+1

    # BLS
    i = 0
    valuations = [64,64,64,64]
    for k in [12,24,15,27]:
        px,rx,tx,cx,yx,betax,lambx, D = tnfs.curve.bls.polynomial_params(k)
        m = 3
        u_mod_m = [1]
        val_l = valuations[i]
        l = 2
        x = px.parent().gen(0)
        print("\nBLS-{} rx = {}".format(k, rx.factor()))
        for vx, v_str, other_poly_str in [(px-1, "px-1", [(rx-1, "(rx-1)"), (px, "px"), (rx, "rx")]), (rx-1, "rx-1", [(px-1, "(px-1)"), (rx, "rx"), (px, "px")])]:
            test_find_seed_congruence_for_adicity_vx_ox(vx, v_str, other_poly_str, u_mod_m, m, val_l, l)
        _u_mod_m = [2]
        for vx, v_str, other_poly_str in [(px(-x)-1, "p(-x)-1", [(rx(-x)-1, "(r(-x)-1)"), (px(-x), "p(-x)"), (rx(-x), "r(-x)")]), (rx(-x)-1, "r(-x)-1", [(px(-x)-1, "(p(-x)-1)"), (rx(-x), "r(-x)"), (px(-x), "p(-x)")])]:
            test_find_seed_congruence_for_adicity_vx_ox(vx, v_str, other_poly_str, _u_mod_m, m, val_l, l)
        i = i+1

def test_find_seed_congruence_for_adicity_ZZpoly():
    # BN polynomials
    px,rx,tx,cx,yx,betax,lambx, D = tnfs.curve.bn.polynomial_params()
    val_l = 64
    l = 2
    x = px.parent().gen(0)
    print("\nBN rx = {} px = {}".format(rx, px))
    for vx, v_str, other_poly_str in [(px-1, "px-1", [(rx-1, "(rx-1)"), (px, "px"), (rx, "rx")]), (px(-x)-1, "p(-x)-1", [(rx(-x)-1, "(r(-x)-1)"), (px(-x), "p(-x)"), (rx(-x), "r(-x)")]), (rx-1, "rx-1", [(px-1, "(px-1)"), (rx, "rx"), (px, "px")]), (rx(-x)-1, "r(-x)-1", [(px(-x)-1, "(p(-x)-1)"), (rx(-x), "r(-x)"), (px(-x), "p(-x)")])]:
        test_find_seed_congruence_for_adicity_ZZ_vx_ox(vx, v_str, other_poly_str, val_l, l)

def test_find_seed_congruence_for_adicity_MNT():
    # MNT polynomials
    for k in [3,4,6]:
        px,rx,tx,cx = tnfs.curve.mnt_k.polynomial_params(k)
        val_l = 64
        l = 2
        x = px.parent().gen(0)
        print("\nMNT{} rx = {} px = {}".format(k, rx, px))
        for vx, v_str, other_poly_str in [(px-1, "px-1", [(rx-1, "(rx-1)"), (px, "px"), (rx, "rx")]), (px(-x)-1, "p(-x)-1", [(rx(-x)-1, "(r(-x)-1)"), (px(-x), "p(-x)"), (rx(-x), "r(-x)")]), (rx-1, "rx-1", [(px-1, "(px-1)"), (rx, "rx"), (px, "px")]), (rx(-x)-1, "r(-x)-1", [(px(-x)-1, "(p(-x)-1)"), (rx(-x), "r(-x)"), (px(-x), "p(-x)")])]:
            test_find_seed_congruence_for_adicity_ZZ_vx_ox(vx, v_str, other_poly_str, val_l, l)

def test_find_seed_congruence_for_adicity_GMV():
    # Galbraith McKee Valenca polynomials
    val_l = 64
    l = 2
    params_gmv = tnfs.curve.galbraith_mckee_valenca.allowed_index
    for k in params_gmv:
        for h in params_gmv[k]:
            for i in range(params_gmv[k][h]):
                px,rx,tx,cx = tnfs.curve.galbraith_mckee_valenca.polynomial_params(k,h,i)
                x = px.parent().gen(0)
                print("\nGMV-k{}-h{}-{} rx = {} px = {} tx = {}".format(k, h, i, rx, px, tx))
                for vx, v_str, other_poly_str in [(px-1, "px-1", [(rx-1, "(rx-1)"), (px, "px"), (rx, "rx")]), (px(-x)-1, "p(-x)-1", [(rx(-x)-1, "(r(-x)-1)"), (px(-x), "p(-x)"), (rx(-x), "r(-x)")]), (rx-1, "rx-1", [(px-1, "(px-1)"), (rx, "rx"), (px, "px")]), (rx(-x)-1, "r(-x)-1", [(px(-x)-1, "(p(-x)-1)"), (rx(-x), "r(-x)"), (px(-x), "p(-x)")])]:
                    test_find_seed_congruence_for_adicity_ZZ_vx_ox(vx, v_str, other_poly_str, val_l, l)

def test_find_seed_congruence_for_adicity_MNTG():
    ## MNT-Generalized polynomials
    val_l = 64
    l = 2
    params_gmv = tnfs.curve.galbraith_mckee_valenca.allowed_index
    for k in params_gmv:
        for h in params_gmv[k]:
            for i in range(params_gmv[k][h]):
                tx_gmv = tnfs.curve.galbraith_mckee_valenca.params_k[k][h][i]['tx']
                phik = cyclotomic_polynomial(k)(tx_gmv-1)
                d = gcd([ZZ(ci) for ci in phik.list()])
                if d > 1:
                    px,rx,tx,cx = tnfs.curve.mntg.polynomial_params(k,h,d)
                    m, u_mod_m = tnfs.curve.mntg.congruence_constraints(k,d)
                    assert m == d
                    ux = tx_gmv-1
                    cu = gcd([ZZ(ci) for ci in ux.list()])
                    lux = ux.leading_coefficient()
                    print("tx-1 = {}, leading_coeff = {} = {}".format(ux, lux, lux.factor()))
                    assert (lux) % d == 0
                    print("m, u_mod_m = {}, {}, tx-1 = {}, cu = {}, (tx-1) % {} = {}".format(m, u_mod_m, tx_gmv-1, cu, m, ZZ((tx_gmv-1).constant_coefficient()) % m))
                    assert ZZ(((tx_gmv-1).constant_coefficient()) % d) in u_mod_m
                    x = px.parent().gen(0)
                    print("\nMNTG-{} c={} d={} rx = {} px = {} tx = {}".format(k, h, d, rx, px, tx))
                    print("")
                    for vx, v_str, other_poly_str in [(px-1, "px-1", [(rx-1, "(rx-1)"), (px, "px"), (rx, "rx")]), (px(-x)-1, "p(-x)-1", [(rx(-x)-1, "(r(-x)-1)"), (px(-x), "p(-x)"), (rx(-x), "r(-x)")]), (rx-1, "rx-1", [(px-1, "(px-1)"), (rx, "rx"), (px, "px")]), (rx(-x)-1, "r(-x)-1", [(px(-x)-1, "(p(-x)-1)"), (rx(-x), "r(-x)"), (px(-x), "p(-x)")])]:
                        # check if QQ is needed or only ZZ?
                        test_find_seed_congruence_for_adicity_vx_ox(vx, v_str, other_poly_str, u_mod_m, m, val_l, l)
            i = i+1

if __name__ == "__main__":
    test_find_seed_congruence_for_adicity_QQpoly()
    test_find_seed_congruence_for_adicity_ZZpoly()
    test_find_seed_congruence_for_adicity_MNT()
    #test_find_seed_congruence_for_adicity_GMV()
    #test_find_seed_congruence_for_adicity_MNTG()
