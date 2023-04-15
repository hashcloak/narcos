from sage.all_cmdline import *   # import sage library
from sage.rings.rational_field import Q, QQ

import tnfs
from tnfs.gen.generate_curve_utils import small_factors, str_py_2NAF, str_2NAF, str_binary_form, str_py_binary_form
from tnfs.gen.generate_curve_utils import find_min_max_list_mod_m
from tnfs.gen.generate_curve_utils import find_seed_congruence_for_adicity_QQpoly, find_seed_congruence_for_adicity_ZZpoly
from tnfs.gen.generate_curve_utils import get_next_seed
import tnfs.curve.pairing_friendly_curve

import tnfs.curve.kss
from tnfs.curve.kss import KSS
import tnfs.curve.bls
from tnfs.curve.bls import BLS
import tnfs.curve.bn
from tnfs.curve.bn import BN
import tnfs.curve.aurifeuillean
from tnfs.curve.aurifeuillean import Aurifeuillean
from tnfs.curve.aurifeuillean import allowed_k as allowed_k_auri

allowed_k_kss = [8,16,18,32,36,40,54]
allowed_k_bls = [6,9,12,24,27,48] # TODO 15, 21, (30, 33, 39, 42, 45)
allowed_k_bn = [12]

# read from command-line the embedding degree k, the sizes of p or r and options
# after a major upgrade of sage somewhere in 2020, the older command-line does
# not work anymore, now using sage -python -m option to launch modules properly.
# example:
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 16 -p 320 1024 -s 64
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 16 -r 248 257  --find_all_u --valuation_r 16
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 16 -r 248 257  --find_all_u --negative_u --valuation_r 16
# the pb with KSS18 is partially fixed.
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 18 -r 254 255  --find_all_u --valuation_r 70
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 18 -r 254 255  --find_all_u --negative_u --valuation_r 70
# sage -python -m tnfs.gen.compute_test_vector_curve --bls -k 12 -p 320 1024 -s 64
# sage -python -m tnfs.gen.compute_test_vector_curve --bls -k 12 -r 254 254 --find_all_u --valuation 40 --check_gt_cofactor
# sage -python -m tnfs.gen.compute_test_vector_curve --bn -r 254 255  --find_all_u --valuation 48
# sage -python -m tnfs.gen.compute_test_vector_curve --bn -r 383 383  --find_all_u --valuation 48
# sage -python -m tnfs.gen.compute_test_vector_curve --bn -r 383 383  --find_all_u --valuation_r 64
# sage -python -m tnfs.gen.compute_test_vector_curve --bn -r 383 383  --find_all_u --valuation 48 --negative_u
# sage -python -m tnfs.gen.compute_test_vector_curve --bn -p 256 1024 -s 64
# sage -python -m tnfs.gen.compute_test_vector_curve --auri -k 18 -p 256 1024 -s 64
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 32 -r 384 400 --find_all_u
# sage -python -m tnfs.gen.compute_test_vector_curve --kss -k 36 -r 384 384 --find_all_u
#
# IF YOU DON'T HAVE A TERMINAL: within Sage prompt:
# sage: run tnfs/gen/compute_test_vector_curve --bls -k 48 # and other options

curve = None # KSS or BLS or BN or Aurifeuillean
k = None
D = None
e0 = None
a_auri = None
length_p_min = None
length_p_max = None
length_r_min = None
length_r_max = None
step_length = None
u_mod_m = None
check_pnbits = False
check_rnbits = False
valuation = None
valuation_p = None
valuation_r = None

# some flags and default value
check_all_cofactors = False
check_twist_order = False
check_gt_cofactor = False
allowed_cofactor = 2**16
allow_cofact_r=False ; str_c = "_r_prime"
factor_r = False
negative_u = False;  str_u= "_pos_u_"
find_all_u = False

args=sys.argv
filename_stem = '.'.join((args[0].split('/')[-1]).split('.')[:-1])
print("usage: sage -python -m tnfs.gen.{} [--kss | --bls | --bn | --auri]  [-k <k>] [-D <D>] [-e0 <e0>] [-p <min_pnbits> <max_pnbits> | -r <min_rnbits> <max_rnbits>] [-s step_nbits -u <u_mod_m> --negative_u --find_all_u --allow_cofactor_r --factor_r] [--check_twist_order --check_gt_cofactor | --check_all_cofactors] [--valuation <v> | --valuation_r <vr> | --valuation_p <vp>]".format(filename_stem))
i=1
cmd_options = ["--kss", "--bls", "--bn", "--auri", "-k", "-D", "-e0", "-a", "-p", "-r", "-s", "-u", "--negative_u", "--find_all_u", "--allow_cofactor_r", "--factor_r", "--check_twist_order", "--check_all_cofactors", "--check_gt_cofactor", "--valuation", "--valuation_r", "--valuation_p"]
while i < len(args):
    if args[i] == "--bls":
        curve = BLS; str_curve = "BLS"; str_curve_lowcase = "bls"
        curve_path = tnfs.curve.bls
        allowed_k = allowed_k_bls
    elif args[i] == "--kss":
        curve = KSS; str_curve = "KSS"; str_curve_lowcase = "kss"
        curve_path = tnfs.curve.kss
        allowed_k = allowed_k_kss
    elif args[i] == "--bn":
        curve = BN; str_curve = "BN"; str_curve_lowcase = "bn"
        curve_path = tnfs.curve.bn
        allowed_k = allowed_k_bn
        k = 12
    elif args[i] == "--auri":
        curve = Aurifeuillean; str_curve = "Auri"; str_curve_lowcase = "auri"
        curve_path = tnfs.curve.aurifeuillean
        allowed_k = allowed_k_auri
    elif args[i] == "-k" and (i+1) < len(args) and args[i+1] not in cmd_options:
        k = Integer(args[i+1])
        i += 1
    elif args[i] == "-D" and (i+1) < len(args) and args[i+1] not in cmd_options:
        D = Integer(args[i+1])
        i += 1
    elif args[i] == "-e0" and (i+1) < len(args) and args[i+1] not in cmd_options:
        e0 = Integer(args[i+1])
        i += 1
    elif args[i] == "-a" and (i+1) < len(args) and args[i+1] not in cmd_options:
        a_auri = Integer(args[i+1])
        i += 1
    elif args[i] == "-p" and (i+1) < len(args) and args[i+1] not in cmd_options:
        length_p_min = Integer(args[i+1])
        if (i+2) < len(args) and args[i+2] not in cmd_options:
            length_p_max = Integer(args[i+2])
            i += 1
        else:
            length_p_max = length_p_min
        i += 1
    elif args[i] == "-r" and (i+1) < len(args) and args[i+1] not in cmd_options:
        length_r_min = Integer(args[i+1])
        if (i+2) < len(args) and args[i+2] not in cmd_options:
            length_r_max = Integer(args[i+2])
            i += 1
        else:
            length_r_max = length_r_min
        i += 1
    elif args[i] == "-s" and (i+1) < len(args) and args[i+1] not in cmd_options:
        step_length = Integer(args[i+1])
        i += 1
    elif args[i] == "-u" and (i+1) < len(args) and args[i+1] not in cmd_options:
        u_mod_m = Integer(args[i+1])
        i += 1
    elif args[i] == "--negative_u":
        negative_u = True;  str_u= "_neg_u_"
    elif args[i] == "--find_all_u":
        find_all_u = True
    elif args[i] == "--factor_r":
        factor_r = True
    elif args[i] == "--allow_cofactor_r":
        allow_cofact_r=True ; str_c = "_r_composite"
    elif args[i] == "--check_all_cofactors":
        check_all_cofactors = True
    elif args[i] == "--check_twist_order":
        check_twist_order = True
    elif args[i] == "--check_gt_cofactor":
        check_gt_cofactor = True
    elif args[i] == "--valuation" and (i+1) < len(args) and args[i+1] not in cmd_options:
        valuation = Integer(args[i+1])
        valuation_p = True
        valuation_r = True
        find_all_u = True
        i += 1
    elif args[i] == "--valuation_r" and (i+1) < len(args) and args[i+1] not in cmd_options:
        valuation = Integer(args[i+1])
        valuation_r = True
        valuation_p = False
        find_all_u = True
        i += 1
    elif args[i] == "--valuation_p" and (i+1) < len(args) and args[i+1] not in cmd_options:
        valuation = Integer(args[i+1])
        valuation_p = True
        valuation_r = False
        find_all_u = True
        i += 1
    i += 1

if curve is None:
    raise ValueError("missing curve type: one of --kss or --bls or --bn or --auri")

if k is None:
    raise ValueError("missing embedding degree k, possible values with {}: {}".format(str_curve, allowed_k))
elif k not in allowed_k:
    raise ValueError("embedding degree k={} but possible values with {} are {}".format(k, str_curve, allowed_k))

print("k = {}".format(k))
print("{} curve".format(str_curve))
if length_p_min is not None and length_p_max is not None:
    print("length_p {} - {} bits".format(length_p_min, length_p_max))
    check_pnbits = True
if length_r_min is not None and length_r_max is not None:
    print("length_r {} - {} bits".format(length_r_min, length_r_max))
    check_rnbits = True
if step_length is not None:
    print("step_length = {}".format(step_length))
else:
    step_length = 32
    print("default step_length = {}".format(step_length))
if u_mod_m is not None:
    print("u_mod_m = {}".format(u_mod_m))

if not check_pnbits and not check_rnbits:
    if curve==KSS:
        if k == 8:
            check_pnbits = True ; check_rnbits = False
            length_p_min = 576
            length_p_max = 608
            step_length = 8
        elif k in [16, 18, 32, 36, 40, 54]:
            check_pnbits = False ; check_rnbits = True
            length_r_min = 256
            length_r_max = 512
            step_length = 128
    elif curve == BLS:
        if k == 6:
            check_pnbits = True ; check_rnbits = False
            length_p_min = 800
            length_p_max = 864
            step_length = 8
        elif k == 9:
            check_pnbits = True ; check_rnbits = False
            length_p_min = 576
            length_p_max = 608
            step_length = 8
        elif k == 12:
            check_pnbits = True ; check_rnbits = False
            length_p_min = 384
            length_p_max = 1024
            step_length = 64
        elif k in [24, 27, 48]:
            check_pnbits = False ; check_rnbits = True
            length_r_min = 256
            length_r_max = 512
            step_length = 128
    elif curve == BN and k == 12:
        check_pnbits = True ; check_rnbits = False
        length_p_min = 256
        length_p_max = 1024
        step_length = 64
    elif curve == Aurifeuillean:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 256
        length_r_max = 512
        step_length = 128
    print("{} k = {}, default values".format(str_curve, k))
    if check_pnbits:
        print("length_p {} - {} bits".format(length_p_min, length_p_max))
    if check_rnbits:
        print("length_r {} - {} bits".format(length_r_min, length_r_max))

rnbits = 200
proof.arithmetic(False)

QQx = QQ['x']; (x,) = QQx._first_ngens(1)
if curve == Aurifeuillean:
    px, rx, tx, cx, yx, betax, lambx, D, a_auri, exp_tr, automorphisms, m, u_m, cofactor_r = curve_path.polynomial_params(k, D, a_auri, e0)
else:
    px, rx, tx, cx, yx, betax, lambx, D = curve_path.polynomial_params(k)

assert tx**2 - 4*px == -D*yx**2
assert px == (tx**2 + D*yx**2)/4
assert (px+1-tx) == rx*cx

if curve == Aurifeuillean:
    twx, g2cx, g2twx, polys_cofactor_twist, label_factors, small_cofactors = tnfs.curve.pairing_friendly_curve.poly_cofactor_twist_g1_g2(k, px, rx, tx, cx, yx, D)
    gt_cx = curve_path.poly_cofactor_gt(k, D, a_auri, e0)
else:
    m, u_m, cofactor_r, twx, g2cx, g2twx, polys_cofactor_twist, label_factors, small_cofactors = curve_path.poly_cofactor_twist_g1_g2(k)
    gt_cx = curve_path.poly_cofactor_gt(k)

if not check_all_cofactors:
    print("not checking cofactors")
    if check_twist_order:
        l = lcm([denom(ci) for ci in twx.list()])
        g = gcd([ZZ(ci) for ci in (l*twx).list()])
        if (twx/g).is_irreducible():
            twx = twx/g
            polys_cofactor_twist_ = [twx]
            label_factors_ = ["twx"]
        else:
            factors_tw = twx.factor()
            polys_cofactor_twist_ = [fi for (fi, _) in factors_tw]
            label_factors_ = ["twx"+str(i) for i in range(len(factors_tw))]
    else:
        polys_cofactor_twist_ = []
        label_factors_ = []
else:
    polys_cofactor_twist_ = polys_cofactor_twist
    label_factors_ = label_factors
if check_gt_cofactor:
    print("check GT cofactor")
    polys_cofactor_twist_.append(gt_cx)
    label_factors_.append("gt_cx")
    polys_cofactor_twist.append(gt_cx)
    label_factors.append("gt_cx")

print("polys_cofactor_twist = {}".format(polys_cofactor_twist))
print("label_factors = {}".format(label_factors))
print("polys_cofactor_twist_ = {}".format(polys_cofactor_twist_))
print("label_factors_ = {}".format(label_factors_))

if u_mod_m is not None:
    if u_mod_m in u_m:
        u_m = [u_mod_m]
        print("u_mod_m = {} mod {}".format(u_mod_m, m))
    else:
        if len(u_m) >= 10:
            raise ValueError("u_mod_m = {} given but allowed values mod {} are (first 10) {}".format(u_mod_m, m, u_m[:10]))
        else:
            raise ValueError("u_mod_m = {} given but allowed values mod {} are {}".format(u_mod_m, m, u_m))

if find_all_u:
    print("finding all u")
if allow_cofact_r:
    print("cofactors allowed for r")

filename = "test_vector_{}{}".format(str_curve_lowcase, k)
if curve == Aurifeuillean:
    if a_auri < 0:
        filename += "_D{}_a_{}_e0{}".format(D, -a_auri, exp_tr)
    else:
        filename += "_D{}_a{}_e0{}".format(D, a_auri, exp_tr)        
if check_pnbits:
    filename += "_pnbits_{}_{}".format(length_p_min, length_p_max)
if check_rnbits:
    filename += "_rnbits_{}_{}".format(length_r_min, length_r_max)
if valuation is not None:
    filename += "_val2_{}".format(valuation)
filename += "{}{}_u_".format(str_c, str_u)
for ui_m in u_m[:4]:
    filename += "{}_".format(ui_m)
filename += "mod_{}.py".format(m)
print(filename)
ofile = open(filename, "a+")

if check_pnbits:
    ofile.write("# {}{} curves with seed u = {} mod {} s.t. p has {} to {} bits\n".format(str_curve, k, u_m[:4], m, length_p_min, length_p_max))
if check_rnbits:
    ofile.write("# {}{} curves with seed u = {} mod {} s.t. r has {} to {} bits\n".format(str_curve, k, u_m[:4], m, length_r_min, length_r_max))
ofile.write("test_vector_{}{} = [\n".format(str_curve, k))

if check_pnbits:
    length_min = length_p_min
    length_max = length_p_max
elif check_rnbits:
    length_min = length_r_min
    length_max = length_r_max

q_nbits = length_min # will be the current size of p or r depending on check_pnbits or check_rnbits
if check_pnbits:
    u_min_pos, u_max_pos, u_min_neg, u_max_neg = find_min_max_list_mod_m(px, length_p_min, length_p_max, m, u_m)
    poly_q = px
elif check_rnbits:
    u_min_pos, u_max_pos, u_min_neg, u_max_neg = find_min_max_list_mod_m(rx, length_r_min, length_r_max, m, u_m)
    poly_q = rx
# if polys are even then negative u is useless

if negative_u: # u_min_neg < u_max_neg < 0 but p(u_min_neg) > p(u_max_neg) > 0
    u_max = u_max_neg
    u_min = u_min_neg
    u_step = -m # work on that when several u_mod_m are allowed, compute the deltas between each
    if find_all_u:
        u_start = u_max # so that the shorter p starts, and p-length increases
        u_stop = u_min
else: # 0 < u_min < u_max
    u_max = u_max_pos
    u_min = u_min_pos
    u_step = m
    if find_all_u:
        u_start = u_min
        u_stop = u_max
# note: we always should have u_min < u_max
print("u_max = {} = {:#x}".format(u_max, u_max))
print("u_min = {} = {:#x}".format(u_min, u_min))

if valuation is not None and find_all_u:
    if negative_u:
        rx_ = rx(-x)
        px_ = px(-x)
    else:
        rx_ = rx
        px_ = px
    if m == 1: # polynomials are in ZZ[x]
        if valuation_p and valuation_r:
            dict_mj_uj_p = find_seed_congruence_for_adicity_ZZpoly(px_-1, valuation, 2)
            dict_mj_uj_r = find_seed_congruence_for_adicity_ZZpoly(rx_-1, valuation, 2)
        elif valuation_p:
            dict_mj_uj = find_seed_congruence_for_adicity_ZZpoly(px_-1, valuation, 2)
        elif valuation_r:
            dict_mj_uj = find_seed_congruence_for_adicity_ZZpoly(rx_-1, valuation, 2)
    else:
        # refine u_max and u_min so that they match the constraint =u_m mod m
        # it means, compute the value u_m mod m*2^val
        # there is only u_m = [1] and m=3 for BLS curves
        # more conditions for KSS curves
        print("{} u = {} mod {}".format(curve, u_m, m))
        if negative_u:
            u_m_ = [m-u_m[j] for j in range(len(u_m)-1, -1, -1)]
        else:
            u_m_ = u_m
        if valuation_p and valuation_r:
            dict_mj_uj_p = find_seed_congruence_for_adicity_QQpoly(px_-1, u_m, m, valuation, 2)
            dict_mj_uj_r = find_seed_congruence_for_adicity_QQpoly(rx_-1, u_m, m, valuation, 2)
        elif valuation_p:
            dict_mj_uj = find_seed_congruence_for_adicity_QQpoly(px_-1, u_m, m, valuation, 2)
        elif valuation_r:
            dict_mj_uj = find_seed_congruence_for_adicity_QQpoly(rx_-1, u_m_, m, valuation, 2)
    # ok
    if valuation_p and valuation_r:
        # compute joint congruences: what if mj for p divides mk for r?
        u_m2 = 0
        print("joint valuations")
        dict_mj_uj = {}
        for mj in dict_mj_uj_p:
            list_uj = []
            for mk in dict_mj_uj_r:
                if mj == mk:
                    list_uj += [uj for uj in dict_mj_uj_p[mj] if uj in dict_mj_uj_r[mj]]
                    new_mj = mj
                elif (mj % mk) == 0: # mk divides mj, mj is the largest, check if (uj % mk) in dict_mj_uj_r
                    print("mj % mk == 0 with mj = {} and mk = {}".format(mj.factor(), mk.factor()))
                    list_uj += [uj for uj in dict_mj_uj_p[mj] if (uj % mk) in dict_mj_uj_r[mk]]
                    new_mj = mj
                elif (mk % mj) == 0:
                    print("mk % mj == 0 with mk = {} and mj = {}".format(mk.factor(), mj.factor()))
                    list_uj += [uk for uk in dict_mj_uj_r[mk] if (uk % mj) in dict_mj_uj_p[mj]]
                    new_mk = mk
            list_uj.sort()
            print("joint congruences are {}".format(list_uj))
            list_uj_unique = [list_uj[0]]
            for j in range(1, len(list_uj)):
                if list_uj[j] != list_uj[j-1]:
                    list_uj_unique.append(list_uj[j])
            if len(list_uj_unique) > 0:
                dict_mj_uj[new_mj] = [uj for uj in list_uj_unique]
            print("{} : {}".format(new_mj, list_uj_unique))
    if negative_u:
        for mj in dict_mj_uj:# account for the case uj = 0 -> mj-uj = mj -> reduce mod mj
            dict_mj_uj[mj] = [(mj - dict_mj_uj[mj][j]) % mj for j in range(len(dict_mj_uj[mj])-1, -1, -1)]
    list_m2 = list(dict_mj_uj)
    #print("list_m2 = {}".format(list_m2))
    #print("congruences:")
    #for mj in dict_mj_uj:
    #    print("{} : {}".format(mj, dict_mj_uj[mj]))
    # filter the ones such that the content of rx(mj*x+uj) is one and the content of px(mj*x+uj) is one
    x = px.parent().gen(0)
    new_dict_mj_uj = {}
    for mj in dict_mj_uj:
        new_list_uj = []
        for uj in dict_mj_uj[mj]:
            try:
                cont_pj = gcd([ZZ(pj) for pj in (px(mj*x + uj)).list()])
                cont_rj = gcd([ZZ(rj) for rj in (rx(mj*x + uj)).list()])
            except TypeError as err:
                print("px(mj*x + uj) = {}".format(px(mj*x + uj)))
                print("rx(mj*x + uj) = {}".format(rx(mj*x + uj)))
                cont_pj = 2
                cont_rj = 2
            if cont_pj == 1 and cont_rj == 1:
                new_list_uj.append(uj)
        if len(new_list_uj) > 0:
            new_dict_mj_uj[mj] = [uj for uj in new_list_uj]
    dict_mj_uj = new_dict_mj_uj
    list_m2 = list(dict_mj_uj)
    print("list_m2 = {}".format(list_m2))
    print("congruences:")
    for mj in dict_mj_uj:
        print("{} : {}".format(mj, dict_mj_uj[mj]))

    list_m2 = [m2 for m2 in list_m2 if m2 < max(abs(u_start), abs(u_stop))]
    print("list_m2 = {} smaller than max(abs(u_start), abs(u_stop)) = {}".format(list_m2, max(abs(u_start), abs(u_stop))))
    m2 = min(list_m2)
    print("chosen m2 = {}".format(m2))
    u_m2 = dict_mj_uj[m2]
    # now any u_i = u_m2[i] mod m2 will satisfy 2^val | poly(u_i)
    # so the idea is to look at u = h*m2 + u_m2 for some h
    #h*m2 + u_m2 >= u_start
    #h*m2 + u_m2 <= u_stop
    ui_start = [ui_m2 for ui_m2 in u_m2 if ui_m2 >= (u_start % m2)]
    if len(ui_start) > 0:
        u_start = u_start - (u_start % m2) + min(ui_start)
    else:
        u_start = u_start - (u_start % m2) + m2 + min(u_m2)
    ui_stop = [ui_m2 for ui_m2 in u_m2 if ui_m2 <= (u_stop % m2)]
    if len(ui_stop) > 0:
        u_stop = u_stop - (u_stop % m2) + max(ui_stop)
    else:
        u_stop = u_stop - (u_stop % m2) - m2 + max(u_m2)

print("u_max = {} = {:#x} {}-valuation={}".format(u_max, u_max, 2, valuation))
print("u_min = {} = {:#x} {}-valuation={}".format(u_min, u_min, 2, valuation))
print("u_start = {0} = {0:#x}".format(u_start))
print("u_stop  = {0} = {0:#x}".format(u_stop))


prod_primes = prod(prime_range(10**7))
# save the gcd of parameters to detect possible systematic cofactor
gcd_ci = [0]* len(polys_cofactor_twist)

ctr = 0
low = True
# find lowest seed in the interval, and highest seed.
# One single call to find_min_max is required to get both bounds
while (find_all_u and (u_max > u_min)) or (not find_all_u and (q_nbits <= length_max)):
    print("q_nbits = {}".format(q_nbits))
    if not find_all_u:
        if low:
            # compute start_u, stop_u for given q_nbits
            u_min_pos, u_max_pos, u_min_neg, u_max_neg = find_min_max_list_mod_m(poly_q, q_nbits, q_nbits, m, u_m)
            if negative_u:
                u_start = u_max_neg
                u_stop = u_min_neg
                u_step = -m
            else:
                u_start = u_min_pos
                u_stop = u_max_pos
                u_step = m
            res = get_next_seed(u_start, u_stop, m, u_m, px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_cofactor_twist_,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
        else: # bounds are already computed
            res = get_next_seed(u_stop, u_start, m, u_m, px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_cofactor_twist_,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
        low = not low # flip for next time
    else:
        print("get_next_seed ...")
        if valuation is not None:
            res = get_next_seed(u_start, u_stop, m2, u_m2, px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_cofactor_twist_,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
        else:
            res = get_next_seed(u_start, u_stop, m, u_m, px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_cofactor_twist_,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
        print("got result")
    u, p, r, cofact_tw_orders, tries = res
    if u == None:
        print("nothing found for size {} and u={} mod {}, u_start = {}, u_stop = {}, u_step = {}".format(q_nbits, u_m, m, u_start, u_stop, u_step))
        if not find_all_u:
            low = True
            q_nbits += step_length
        else:
            u_max = u_min-abs(u_step) # this will stop the loop
        continue
    print("u_min   = {}".format(u_min))
    print("u_max   = {}".format(u_max))
    print("u_start = {}".format(u_start))
    print("u_stop  = {}".format(u_stop))
    print("u_step  = {}".format(u_step))
    if allow_cofact_r:
        cofact_r = cofact_tw_orders[0]
        print("u = {:#x} #{} bits\np = {:#x} #{} bits\nr = {:#x} #{} bits\ncofact_r = {} #{} bits".format(u, u.nbits(), p, p.nbits(), r, r.nbits(), cofact_r, cofact_r.nbits()))
    else:
        cofact_r = 1
        print("u = {:#x} #{} bits\np = {:#x} #{} bits\nr = {:#x} #{} bits, prime".format(u, (ZZ(u)).nbits(), p, p.nbits(), r, r.nbits()))
    if valuation is not None:
        if valuation_p:
            print("p-1 = {:#x} 2-valuation {}".format(p-1, (p-1).valuation(2)))
        if valuation_r:
            print("r-1 = {:#x} 2-valuation {}".format(r-1, (r-1).valuation(2)))

    # recompute cofactors
    cofact_tw_orders = [ZZ(c_i(u)) for c_i in polys_cofactor_twist]
    # print cofactors and twist orders
    gcd_c_tw = [0 for _ in range(len(polys_cofactor_twist))]
    j=0
    for ci in cofact_tw_orders:
        small_cofactors_ci, cii = small_factors(ci, prod_primes)
        gcd_ci[j] = gcd(gcd_ci[j], small_cofactors_ci)
        gcd_c_tw[j] = small_cofactors_ci
        print("{} = {}\n  = ({}) * {} {} bits, {} bits".format(label_factors[j], ci, small_cofactors_ci, cii, small_cofactors_ci.nbits(), cii.nbits()))
        print("{}0 = {} = {}".format(label_factors[j], small_cofactors_ci, small_cofactors_ci.factor()))
        print("{}1 is prime: {} {} bits".format(label_factors[j], cii.is_pseudoprime(), cii.nbits()))
        j+=1

    print("u ={: #x}".format(u))
    print("tries = {}".format(tries))

    # find curve parameter a or b
    ok = False
    if (D != 1) and (D != 4) and (D != 3):
        E = curve(k=k, u=u, D=D, e0=exp_tr, aa=a_auri, cofactor_r=cofact_r)
        ok = True
    ab=1
    while not ok:
        try:
            if D == 3:
                if curve == Aurifeuillean:
                    E = curve(k=k, u=u, D=D, e0=exp_tr, aa=a_auri, b=ab, cofactor_r=cofact_r)
                else:
                    E = curve(k=k, u=u, b=ab, cofactor_r=cofact_r)
            elif D == 1 or D == 4:
                if curve == Aurifeuillean:
                    E = curve(k=k, u=u, D=D, e0=exp_tr, aa=a_auri, a=ab, cofactor_r=cofact_r)
                else:
                    E = curve(k=k, u=u, a=ab, cofactor_r=cofact_r)
            ok=True
        except ValueError as err:
            try:
                if D == 3:
                    if curve == Aurifeuillean:
                        E = curve(k=k, u=u, D=D, e0=exp_tr, aa=a_auri, b=-ab, cofactor_r=cofact_r)
                    else:
                        E = curve(k=k, u=u, b=-ab, cofactor_r=cofact_r)
                elif D == 1 or D == 4:
                    if curve == Aurifeuillean:
                        E = curve(k=k, u=u, D=D, e0=exp_tr, aa=a_auri, a=-ab, cofactor_r=cofact_r)
                    else:
                        E = curve(k=k, u=u, a=-ab, cofactor_r=cofact_r)
                ok=True
                ab = -ab
            except ValueError as err:
                ab += 1
        if (ab % 10) == 0:
            print("D={} ab = {}".format(D, ab))
    print(E)
    str_u, Hw_2NAF_u = str_2NAF(u)
    if allow_cofact_r:
        cofact_r_str = "\'cofact_r\':{:7d}, ".format(cofact_r)
    else:
        cofact_r_str = ""
    if valuation_p or valuation_r:
        val_str = "\'val\':{:d}, ".format(2)
        if valuation_p:
            val_str += "\'valp\':{:d}, ".format((p-1).valuation(2))
        if valuation_r:
            val_str += "\'valr\':{:d}, ".format((r-1).valuation(2))
    else:
        val_str = ""
    if find_all_u:
        label_u = "\'label\':\"u={:34s} Hw2NAF={}\"".format(str_u, Hw_2NAF_u)
    else:
        label_u = "\'label\':\"{}{}-{}\"".format(str_curve, k, p.nbits())
    if D==3:
        ofile.write("    {{\'u\':{:#7x}, \'b\':{:2d}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}{}\'deg_h_S\':None,\'cost_S\':None, {}}},#{}\n".format(u, ab, p.nbits(), r.nbits(), cofact_r_str, val_str, label_u, ctr))
    elif D==1 or D==4:
        ofile.write("    {{\'u\':{:#7x}, \'a\':{:2d}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}{}\'deg_h_S\':None,\'cost_S\':None, {}}},#{}\n".format(u, ab, p.nbits(), r.nbits(), cofact_r_str, val_str, label_u, ctr))
    else:
        ofile.write("    {{\'u\':{:#7x}, \'a\':{:2d}, \'b\':{:2d}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}{}\'deg_h_S\':None,\'cost_S\':None, {}}},#{}\n".format(u, E._a, E._b, p.nbits(), r.nbits(), cofact_r_str, val_str, label_u, ctr))
    ofile.flush()
    # find next value of u_max
    if not find_all_u:
        if low:
            q_nbits += step_length
    else:
        #u_start = u + u_step # here, todo: change so that it goes to the next possible value mod m
        if valuation is not None:
            u_mod_m = u % m2
            idx_u_mod_m = u_m2.index(u_mod_m)
            if u > 0:
                if idx_u_mod_m == len(u_m2) - 1:
                    u_start = u - u_mod_m + m2 + u_m2[0]
                else:
                    u_start = u - u_mod_m + u_m2[idx_u_mod_m+1]
            else:
                if idx_u_mod_m == 0:
                    u_start = u - u_mod_m - m2 + u_m2[-1]
                else:
                    u_start = u - u_mod_m + u_m2[idx_u_mod_m-1]
        elif m == 1 or len(u_m) == 1:
            if u > 0:
                u_start = u + m
            else:
                u_start = u - m
        else:
            u_mod_m = u % m
            idx_u_mod_m = u_m.index(u_mod_m)
            if u > 0:
                if idx_u_mod_m == len(u_m) - 1:
                    u_start = u - u_mod_m + m + u_m[0]
                else:
                    u_start = u - u_mod_m + u_m[idx_u_mod_m+1]
            else:
                if idx_u_mod_m == 0:
                    u_start = u - u_mod_m - m + u_m[-1]
                else:
                    u_start = u - u_mod_m + u_m[idx_u_mod_m-1]
        
        if negative_u:
            u_max = u_start
        else:
            u_min = u_start
    ctr+=1
ofile.write("]\n")

ofile.flush()
ofile.close()

print(filename)
