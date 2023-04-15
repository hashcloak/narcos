from sage.all_cmdline import *   # import sage library
from sage.rings.rational_field import Q, QQ

import tnfs
import tnfs.curve.galbraith_mckee_valenca
from tnfs.curve.galbraith_mckee_valenca import allowed_k, allowed_h, allowed_index, params_k
from tnfs.curve.galbraith_mckee_valenca import polynomial_params
from tnfs.gen.generate_curve_utils import find_seed_congruence_for_adicity_ZZpoly
# sage -python -m tnfs.curve.test_galbraith_mckee_valenca

def test_polynomials():
    for k in allowed_k:
        for h in allowed_h[k]:
            for index in range(allowed_index[k][h]):
                item = params_k[k][h][index]
                px = item['px']
                tx = item['tx']
                h1 = item['h'] # corrected h
                order = px+1-tx
                dyx = tx**2-4*px
                #d = lcm([di.denominator() for di in order.list()])
                #c = gcd([ci.numerator() for ci in (d*order).list()])
                c = gcd([ZZ(ci) for ci in order.list()])
                if (c % h1) == 0:
                    rx = order // h1
                    if c != h1:
                        print("k={} h={:2d} h1={:2d} px = {:20s} tx = {:10s} px+1-tx={:24s} rx = {}*({:17s}) tx^2-4*px = {} = {}".format(k, h, h1, str(px), str(tx), str(order), c//h, str((rx//(c//h))), dyx, dyx.factor()))
                    else:
                        print("k={} h={:2d} h1={:2d} px = {:20s} tx = {:10s} px+1-tx={:24s} rx = {:21s} tx^2-4*px = {} = {}".format(k, h, h1, str(px), str(tx), str(order), str(rx), dyx, dyx.factor()))
                        
                else:
                    rx = order/h1
                    print("problem h1={:2d} does not divide px+1-tx".format(h1))
                    print("k={} h={:2d} h1={:2d} px = {:20s} tx = {:10s} px+1-tx={:24s} rx = {:21s} tx^2-4*px = {} = {}".format(k, h, h1, str(px), str(tx), str(order), str(rx), dyx, dyx.factor()))

def test_congruence(k, h, idx, val):
    item = params_k[k][h][idx]
    px = item['px']
    tx = item['tx']
    h1 = item['h'] # corrected h
    order = px+1-tx
    rx = order // h1
    dyx = tx**2-4*px
    print("k={} h={:2d} px = {:20s} tx = {:10s} px+1-tx={:24s} rx = {:21s} tx^2-4*px = {} = {}".format(k, h1, str(px), str(tx), str(order), str(rx), dyx, dyx.factor()))
    dict_mj_uj_p = find_seed_congruence_for_adicity_ZZpoly(px-1, val, 2)
    dict_mj_uj_r = find_seed_congruence_for_adicity_ZZpoly(rx-1, val, 2)
    print("2^{} | p-1: {}".format(val, dict_mj_uj_p))
    print("2^{} | r-1: {}".format(val, dict_mj_uj_r))

def test_all_congruences(val):
    for k in allowed_k:
        for h in allowed_h[k]:
            for index in range(allowed_index[k][h]):
                test_congruence(k, h, index, val)

if __name__ == "__main__":
    test_polynomials()
    print("")
    test_congruence(6, 4, 3, 64)
