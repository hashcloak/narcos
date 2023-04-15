from sage.all_cmdline import *   # import sage library
from sage.rings.rational_field import QQ

# file enumerate_sparse_T.py from https://gitlab.inria.fr/smasson/cocks-pinch-variant/
import tnfs
from tnfs.simul.enumerate_sparse_T import get_sparse_T_HW_gen, get_sparse_T_uptoHW_gen, get_sparse_T_HW_NAF_gen, get_sparse_T_uptoHW_NAF_gen
from tnfs.gen.generate_curve_utils import small_factors, str_py_2NAF, str_2NAF, str_binary_form, str_py_binary_form
from tnfs.gen.generate_curve_utils import find_u_mod_m, find_min_max_list_mod_m
from tnfs.simul.cpx import guess_cost_STNFS_BLS, guess_cost_ConjTNFS_BLS

import tnfs.curve.kss
from tnfs.curve.kss import KSS
from tnfs.curve.kss import allowed_k as allowed_k_kss
import tnfs.curve.bls
from tnfs.curve.bls import BLS
from tnfs.curve.bls import allowed_k as allowed_k_bls
import tnfs.curve.bn
from tnfs.curve.bn import BN
from tnfs.curve.bn import allowed_k as allowed_k_bn
import tnfs.curve.aurifeuillean
from tnfs.curve.aurifeuillean import Aurifeuillean
from tnfs.curve.aurifeuillean import allowed_k as allowed_k_auri
import tnfs.curve.fotiadis_martindale
from tnfs.curve.fotiadis_martindale import FotiadisMartindale
from tnfs.curve.fotiadis_martindale import allowed_code as allowed_code_fm
from tnfs.curve.fst62 import FST62
from tnfs.curve.fst63 import FST63
from tnfs.curve.fst64 import FST64
from tnfs.curve.fst66 import FST66
from tnfs.curve.fst67 import FST67

import tnfs.curve.pairing_friendly_curve
from tnfs.curve.pairing_friendly_curve import compute_a_b

# read from command-line the embedding degree k, the sizes of p or r and options
# after a major upgrade of sage somewhere in 2020, the older command-line does
# not work anymore, now using sage -python -m option to launch modules properly.
# examples:
# sage -python -m tnfs.gen.generate_sparse_curve --bn -r 250 256 --valuation 32
#
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 48 -r 250 256 --2NAF --find_all_w_up_to -w 5
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 27 -r 384 384 --2NAF --find_all_w_up_to -w 6
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 24 -r 256 256 --2NAF --find_all_w_up_to -w 6
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 21 -r 384 384 --2NAF --find_all_w_up_to -w 4
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 15 -p 928 928 -w 6
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 12 -r 384 384 --2NAF --find_all_w_up_to -w 4
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 9  -p 608 608 --2NAF --find_all_w_up_to -w 4
# sage -python -m tnfs.gen.generate_sparse_curve --bls -k 6  -r 384 384 --2NAF --find_all_w_up_to -w 4
#
# sage -python -m tnfs.gen.generate_sparse_curve --kss -k 16 -r 256 256 -w 6 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --kss -k 16 -p 760 768 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --kss -k 18 -r 256 256 -w 6 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --kss -k 36 -r 384 384 -w 8 --find_all_w_up_to --2NAF
#
# sage -python -m tnfs.gen.generate_sparse_curve --auri -k 18 -r 256 256 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --auri -k 20 -D 1 -a 2 -e0 19 -r 380 384 -w 4 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --auri -k 20 -D 1 -a 2 -e0 9 -r 380 384 -w 4 --find_all_w_up_to --2NAF
#
# sage -python -m tnfs.gen.generate_sparse_curve --fm -code 17 -p 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fm -code 23 -r 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fm -code 25 -r 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fm -code 27 -r 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fm -code 28 -r 384 384 -w 5 --find_all_w_up_to --2NAF
#
# sage -python -m tnfs.gen.generate_sparse_curve --fst62 -k 7 -r 256 256 -w 6 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fst63 -k 14 -r 620 620 -w 4 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fst63 -k 22 -r 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fst64 -k 20 -r 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fst64 -k 28 -r 384 384 -w 4 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fst66 -k 24 -r 384 384 -w 5 --find_all_w_up_to --2NAF
# sage -python -m tnfs.gen.generate_sparse_curve --fst67 -k 24 -r 384 384 -w 5 --find_all_w_up_to --2NAF
#
# IF YOU DON'T HAVE A TERMINAL: within Sage prompt:
# sage: run tnfs/gen/generate_sparse_curve --bls -k 48 # and other options


curve = None # KSS or BLS or BN or Aurifeuillean or FotiadisMartindale or FST62, FST63, FST64, FST66, FST67
k = None
allowed_k = None
code = None # for FotiadisMartindale only
D = None
e0 = None
a_auri = None
length_p_min = None
length_p_max = None
length_r_min = None
length_r_max = None
u_mod_m = None
check_pnbits = False # check that the length of p in bits matches length_p
check_rnbits = False # check that the length of r in bits matches length_r

# for debug purpose
#check_size_only = True
check_size_only = False

# some flags and default value
check_all_cofactors = False
check_twist_order = False
check_gt_cofactor = False
allowed_cofactor = 2**16
find_all_Hw_up_to = False
two_naf = False; str_w = "Hw"
u_hw = 5
valid_curve = False
valuation = None # 2-valuation

args=sys.argv
filename_stem = '.'.join((args[0].split('/')[-1]).split('.')[:-1])
print("usage: sage -python -m tnfs.gen.{} [--kss | --bls | --bn | --auri | --fm | --fst62 | --fst63 | --fst64 | --fst66 | --fst67] [-code <code FM>] [-k <k>] [-D <D>] [-e0 <e0>] [-a <auri_a>] [-p <min_pnbits> <max_pnbits> | -r <min_rnbits> <max_rnbits>] [-u <u_mod_m>] -w <Hamming weight in binary form or 2-NAF> --find_all_w_up_to --2NAF --valuation <minimal 2-valuation of both r-1 and p-1>]".format(filename_stem))
i=1
cmd_options = ["--kss", "--bls", "--bn", "--auri", "--fm", "--fst62", "--fst63", "--fst64", "--fst66", "--fst67", "-code", "-k", "-D", "-e0", "-a", "-p", "-r", "-u", "-w", "--find_all_w_up_to", "--2NAF", "--valuation"]

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
    elif args[i] == "--fm":
        curve = FotiadisMartindale; str_curve = "FM"; str_curve_lowcase = "fm"
        curve_path = tnfs.curve.fotiadis_martindale
        allowed_code = allowed_code_fm
    elif args[i] == "--fst62":
        curve = FST62; str_curve = "FST62_k"; str_curve_lowcase = "fst62_k"
        curve_path = tnfs.curve.fst62
    elif args[i] == "--fst63":
        curve = FST63; str_curve = "FST63_k"; str_curve_lowcase = "fst63_k"
        curve_path = tnfs.curve.fst63
    elif args[i] == "--fst64":
        curve = FST64; str_curve = "FST64_k"; str_curve_lowcase = "fst64_k"
        curve_path = tnfs.curve.fst64
    elif args[i] == "--fst66":
        curve = FST66; str_curve = "FST66_k"; str_curve_lowcase = "fst66_k"
        curve_path = tnfs.curve.fst66
    elif args[i] == "--fst67":
        curve = FST67; str_curve = "FST67_k"; str_curve_lowcase = "fst67_k"
        curve_path = tnfs.curve.fst67
    elif args[i] == "-k" and (i+1) < len(args) and args[i+1] not in cmd_options:
        k = Integer(args[i+1])
        i += 1
    elif args[i] == "-code" and (i+1) < len(args) and args[i+1] not in cmd_options:
        code = Integer(args[i+1])
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
        check_pnbits = True
        if (i+2) < len(args) and args[i+2] not in cmd_options:
            length_p_max = Integer(args[i+2])
            i += 1
        else:
            length_p_max = length_p_min
        i += 1
    elif args[i] == "-r" and (i+1) < len(args) and args[i+1] not in cmd_options:
        length_r_min = Integer(args[i+1])
        check_rnbits = True
        if (i+2) < len(args) and args[i+2] not in cmd_options:
            length_r_max = Integer(args[i+2])
            i += 1
        else:
            length_r_max = length_r_min
        i += 1
    elif args[i] == "-u" and (i+1) < len(args) and args[i+1] not in cmd_options:
        u_mod_m = Integer(args[i+1])
        i += 1
    elif args[i] == "-w" and (i+1) < len(args) and args[i+1] not in cmd_options:
        u_hw = Integer(args[i+1])
        i += 1
    elif args[i] == "--find_all_w_up_to":
        find_all_Hw_up_to = True
    elif args[i] == "--2NAF":
        two_naf = True; str_w = "Hw2naf"
    elif args[i] == "--valuation" and (i+1) < len(args) and args[i+1] not in cmd_options:
        valuation = Integer(args[i+1])
        i += 1
    i += 1

if curve is None:
    raise ValueError("missing curve type: one of --kss --bls --bn --auri --fm --fst62 --fst63 --fst64 --fst66 --fst67")

if curve == FotiadisMartindale:
    if code is None:
        raise ValueError("missing family code, possible values with {}: {}".format(str_curve, allowed_code))
    elif code not in allowed_code:
        raise ValueError("family code ={} but possible values with {} are {}".format(code, str_curve, allowed_code))
elif k is None:
    raise ValueError("missing embedding degree k, possible values with {}: {}".format(str_curve, allowed_k))
elif allowed_k is not None and k not in allowed_k:
    raise ValueError("embedding degree k={} but possible values with {} are {}".format(k, str_curve, allowed_k))

if length_p_min is not None and length_p_max is not None:
    print("length_p {} - {} bits".format(length_p_min, length_p_max))
    check_pnbits = True
if length_r_min is not None and length_r_max is not None:
    print("length_r {} - {} bits".format(length_r_min, length_r_max))
    check_rnbits = True
if u_mod_m is not None:
    print("u_mod_m = {}".format(u_mod_m))

# some default values
if not check_pnbits and not check_rnbits and curve == KSS:
    if k == 8:
        check_pnbits = True ; check_rnbits = False
        length_p_min = 576
        length_p_max = 576
        u_hw = 5
    elif k == 16:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 256
        length_r_max = 256
        u_hw = 5
    elif k == 18:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 256
        length_r_max = 256
        u_hw = 5
    elif k == 32:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 384
        length_r_max = 384
        u_hw = 5
    elif k == 36:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 384
        length_r_max = 384
        u_hw = 5
    elif k == 40:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 384
        length_r_max = 384
        u_hw = 5
    elif k == 54:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 512
        length_r_max = 512
        u_hw = 5
elif not check_pnbits and not check_rnbits and curve == BLS:
    if k == 6:
        check_pnbits = True ; check_rnbits = False
        length_p_min = 832
        length_p_max = 832
        u_hw = 5
    elif k == 9:
        check_pnbits = True ; check_rnbits = False
        length_p_min = 576
        length_p_max = 608
        u_hw = 5
    elif k == 12:
        check_pnbits = True ; check_rnbits = False
        length_p_min = 440
        length_p_max = 448
        u_hw = 5
    elif k == 15:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 384
        length_r_max = 384
        u_hw = 6
    elif k == 21:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 384
        length_r_max = 384
        u_hw = 5
    elif k == 24:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 250
        length_r_max = 256
        u_hw = 5
    elif k == 27:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 384
        length_r_max = 384
        u_hw = 7
    elif k == 48:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 250
        length_r_max = 256
        u_hw = 5
elif not check_pnbits and not check_rnbits and curve == Aurifeuillean:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 256
        length_r_max = 512
        step_length = 128
elif not check_pnbits and not check_rnbits and curve == BN and k == 12:
    check_pnbits = True ; check_rnbits = False
    length_p_min = 440
    length_p_max = 448
    u_hw = 5
elif not check_pnbits and not check_rnbits and curve == FotiadisMartindale:
    if code <= 20:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 250
        length_r_max = 256
        u_hw = 5
    elif code >= 21:
        check_pnbits = False ; check_rnbits = True
        length_r_min = 380
        length_r_max = 384
        u_hw = 5
elif not check_pnbits and not check_rnbits:
    check_pnbits = False ; check_rnbits = True
    length_r_min = 256
    length_r_max = 256
    u_hw = 5

proof.arithmetic(False)

QQx = QQ['x']; (x,) = QQx._first_ngens(1)
if curve == Aurifeuillean:
    px, rx, tx, cx, yx, betax, lambx, D, a_auri, exp_tr, automorphisms, m, u_m, cofactor_r = curve_path.polynomial_params(k, D, a_auri, e0)
elif curve == FotiadisMartindale:
    px, rx, tx, cx, yx, betax, lambx, D, k, m, u_m = curve_path.polynomial_params_from_code(code)
else:
    px, rx, tx, cx, yx, betax, lambx, D = curve_path.polynomial_params(k)

assert tx**2 - 4*px == -D*yx**2
assert px == (tx**2 + D*yx**2)/4
assert (px+1-tx) == rx*cx

if curve == Aurifeuillean:
    twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors = tnfs.curve.pairing_friendly_curve.poly_cofactor_twist_g1_g2(k, px, rx, tx, cx, yx, D)
    gt_cx = curve_path.poly_cofactor_gt(k, D, a_auri, e0)
elif curve == FotiadisMartindale:
    twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors = tnfs.curve.pairing_friendly_curve.poly_cofactor_twist_g1_g2(k, px, rx, tx, cx, yx, D)
elif curve in [FST62, FST63, FST64, FST66, FST67]:
    twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors = tnfs.curve.pairing_friendly_curve.poly_cofactor_twist_g1_g2(k, px, rx, tx, cx, yx, D)
    if curve == FST66 or curve == FST67:
        m, u_m = curve_path.get_m_u_mod_m(k)
    else:
        m, u_m = curve_path.get_m_u_mod_m()
else:
    m, u_m, cofactor_r, twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors = curve_path.poly_cofactor_twist_g1_g2(k)
    gt_cx = curve_path.poly_cofactor_gt(k)

if u_mod_m is not None:
    if u_mod_m in u_m:
        u_m = [u_mod_m]
        print("u_mod_m = {} mod {}".format(u_mod_m, m))
    else:
        raise ValueError("u_mod_m = {} given but allowed values mod {} are {}".format(u_mod_m, m, u_m[:4]))

if two_naf:
    print("Considering 2-NAF")
if find_all_Hw_up_to:
    print("finding all u of Hamming weight up to {}".format(u_hw))
else:
    print("finding all u of Hamming weight exactly {}".format(u_hw))

print("k = {}, D = {}".format(k, D))
print("twx= {}\n   = {}".format(twx, twx.factor()))
print("cx = {}\n   = {}".format(cx, cx.factor()))

label_poly = ["px", "rx", "tx", "cx", "twx", "yx", "betax", "lambx"]
i=0
for poly in [px, rx, tx, cx, twx, yx, betax, lambx]:
    den = lcm([ci.denom() for ci in poly.list()])
    if den == 1:
        print("{} = {}".format(label_poly[i], poly))
    else:
        print("{} = ({})/{}".format(label_poly[i], poly*den, den))
    i += 1

if curve == BLS or curve == BN or curve == KSS:
    print("to get integer px(u): u % m, content(px(u*m+i)), m={}".format(m))
    print(find_u_mod_m(px, m))
    print("to get integer rx(u): u % m, content(rx(u*m+i)), m={}".format(m))
    print(find_u_mod_m(rx, m))
    if cx.degree() > 0:
        print("to get integer cx(u): u % m, content(cx(u*m+i)), m={}".format(m))
        print(find_u_mod_m(cx, m))
    else:
        print("cx = {}".format(cx))
    
    g2cx_den = lcm([fi.denom() for fi in g2cx.list()])
    print("G2_cofactor=({})/({})".format(g2cx*g2cx_den, g2cx_den.factor()))
    print("irreducible: {}".format(g2cx.is_irreducible()))
    print("G2_twist = {}\n    = {}".format(g2twx, g2twx.factor()))
    print("irreducible: {}".format(g2twx.is_irreducible()))

    print("to get integer g2cx(u): u % m, content(g2cx(u*m+i)), m={}".format(m))
    print(find_u_mod_m(g2cx, m))
    print("to get integer g2twx(u): u % m, content(g2twx(u*m+i)), m={}".format(m))
    print(find_u_mod_m(g2twx, m))
else:
    print("G2_cofactor is irreducible: {}".format(g2cx.is_irreducible()))
    print("G2_twist is irreducible: {}".format(g2twx.is_irreducible()))


prod_primes = prod(prime_range(10**7))

if check_pnbits:
    u_min_pos, u_max_pos, u_min_neg, u_max_neg = find_min_max_list_mod_m(px, length_p_min, length_p_max, m, u_m)
    poly_q = px
else:
    u_min_pos, u_max_pos, u_min_neg, u_max_neg = find_min_max_list_mod_m(rx, length_r_min, length_r_max, m, u_m)
    poly_q = rx

# consider the widest possible range of positive and negative u
u_min = min(u_min_pos, -u_max_neg)
u_max = max(-u_min_neg, u_max_pos)

u_nbits_min = u_min.nbits()
u_nbits_max = u_max.nbits()
# if the two most significants bits are 1,1 then in 2-naf, they are 1,0,-1 with one more bit
# examples: 3, 6, 7, 0xc, 0xd, 0xe, 0xf
if two_naf and u_min.digits(2)[-2] == 1:
    u_nbits_min += 1
if two_naf and u_max.digits(2)[-2] == 1:
    u_nbits_max += 1

print("u_min_pos = {:#x} # =-2^{:.2f}".format(u_min_pos, float(RR(log(u_min_pos,2.0)))))
print("u_max_pos = {:#x} # =-2^{:.2f}".format(u_max_pos, float(RR(log(u_max_pos,2.0)))))
print("u_max_pos - u_min_pos = {} = 2^{:.2f}".format(u_max_pos-u_min_pos, float(RR(log(abs(u_max_pos-u_min_pos),2)))))
print([i for i in range(len(u_min_pos.digits(2))-1, -1, -1) if u_min_pos.digits(2)[i] == 1])
print([i for i in range(len(u_max_pos.digits(2))-1, -1, -1) if u_max_pos.digits(2)[i] == 1])

print("u_min_neg = {:#x} # =-2^{:.2f}".format(u_min_neg, float(RR(log(-u_min_neg,2.0)))))
print("u_max_neg = {:#x} # =-2^{:.2f}".format(u_max_neg, float(RR(log(-u_max_neg,2.0)))))
print("u_max_neg - u_min_neg = {} = 2^{:.2f}".format(u_max_neg-u_min_neg, float(RR(log(abs(u_max_neg-u_min_neg),2)))))
print([i for i in range(len((-u_min_neg).digits(2))-1, -1, -1) if (-u_min_neg).digits(2)[i] == 1])
print([i for i in range(len((-u_max_neg).digits(2))-1, -1, -1) if (-u_max_neg).digits(2)[i] == 1])

print("u_min = {:#x} # = 2^{:.2f}".format(u_min, float(RR(log(u_min,2.0)))))
print("u_max = {:#x} # = 2^{:.2f}".format(u_max, float(RR(log(u_max,2.0)))))
print("u_max - u_min = {} = 2^{:.2f}".format(u_max-u_min, float(RR(log(u_max-u_min,2)))))

if valuation is not None:
    # assume that one wants 2^v | p-1 and 2^v | r-1 for BLS and BN curves
    if curve == BN:
        # then it means 2^(v-1) divides u because Gcd(p-1, r-1) = 2 * 3 * u
        u_min = u_min // (2**(valuation-1))
        u_max = u_max // (2**(valuation-1))
    elif curve == BLS:
        # then it means 2^v divides (u-1) because Gcd(p-1, r-1) = (u-1)
        u_min = u_min // (2**(valuation))
        u_max = u_max // (2**(valuation))

print("")
count_u = 0
count_valid_u = 0
count_valid_p_size = 0
count_valid_r_size = 0
print("u_nbits_min = {}, u_nbits_max = {}, u_hw = {}, 2-NAF: {}".format(u_nbits_min, u_nbits_max, u_hw, two_naf))

if not check_size_only:
    filename = "test_vector_sparse_{}{}".format(str_curve_lowcase, k)
    if curve == Aurifeuillean:
        if a_auri < 0:
            filename += "_D{}_a_{}_e0{}".format(D, -a_auri, exp_tr)
        else:
            filename += "_D{}_a{}_e0{}".format(D, a_auri, exp_tr)
    elif curve == FotiadisMartindale:
        filename += "_{}".format(code)
    testvector_name = filename
    if check_pnbits:
        filename += "_pnbits_{}_{}".format(length_p_min, length_p_max)
    if check_rnbits:
        filename += "_rnbits_{}_{}".format(length_r_min, length_r_max)
    filename += "_u"
    for ui_m in u_m[:4]:
        filename += "_{}".format(ui_m)
    filename += "_mod_{}".format(m)
    filename += "_unbits_{}".format(u_nbits_min)
    if u_nbits_max > u_nbits_min:
        filename += "_{}".format(u_nbits_max)
    filename += "_{}_{}.py".format(str_w, u_hw) # where str_w = "Hw" or "Hw2naf"
    print(filename)
    ofile = open(filename, "a+")
    if check_pnbits:
        ofile.write("# {}{} curves with sparse seed u = {} mod {} of {}--{} bits {} {} and s.t. p has {} to {} bits\n".format(str_curve, k, u_m[:4], m, u_nbits_min, u_nbits_max, str_w, u_hw, length_p_min, length_p_max))
    else:
        ofile.write("# {}{} curves with sparse seed u = {} mod {} of {}--{} bits {} {} and s.t. r has {} to {} bits\n".format(str_curve, k, u_m[:4], m, u_nbits_min, u_nbits_max, str_w, u_hw, length_r_min, length_r_max))
    ofile.write(testvector_name + " = [\n")

# save the gcd of parameters to detect possible systematic cofactor
gcd_p = 0
gcd_r = 0
gcd_ci = [0]* len(polys_cofact_twists)

avrg_log2_r = RR(0.0)
avrg_log2_p = RR(0.0)

for u_nbits in range(u_nbits_min, u_nbits_max+1):
    if two_naf and find_all_Hw_up_to:
        iterator_u = get_sparse_T_uptoHW_NAF_gen(u_nbits, u_hw)
        print("for u in get_sparse_T_uptoHW_NAF_gen({}, {})".format(u_nbits, u_hw))
    elif two_naf and not find_all_Hw_up_to:
        iterator_u = get_sparse_T_HW_NAF_gen(u_nbits, u_hw)
        print("for u in get_sparse_T_HW_NAF_gen({}, {})".format(u_nbits, u_hw))
    elif not two_naf and find_all_Hw_up_to:
        iterator_u = get_sparse_T_uptoHW_gen(u_nbits, u_hw)
        print("for u0 in get_sparse_T_uptoHW_gen({}, {})".format(u_nbits, u_hw))
    elif not two_naf and not find_all_Hw_up_to:
        iterator_u = get_sparse_T_HW_gen(u_nbits, u_hw)
        print("for u0 in get_sparse_T_HW_gen({}, {})".format(u_nbits, u_hw))

    #for u0 in get_sparse_T_uptoHW_gen(u_nbits, u_hw):
    #for u0 in get_sparse_T_HW_gen(u_nbits, u_hw):
    #for u0 in get_sparse_T_HW_NAF_gen(u_nbits, u_hw):
    #for u0 in get_sparse_T_uptoHW_NAF_gen(u_nbits, u_hw):
    for u0 in iterator_u:
        count_u += 1
        list_u = [u0, -u0]
        # eliminate seeds that have wrong congruence
        list_u = [ui for ui in list_u if (m == 1) or ((ui % m) in u_m)]
        if check_pnbits or check_rnbits:
            list_u = [ui for ui in list_u if (ui < 0 and ui > u_min_neg and ui < u_max_neg) or (ui > 0 and ui > u_min_pos and ui < u_max_pos)]
        for u in list_u:
            count_valid_u += 1
            p = ZZ(px(ZZ(u)))
            avrg_log2_p += RR(log(RR(p),RR(2.0)))
            
            if (count_valid_u % 10000) == 0:
                print("count_valid_u = {}, count_valid_p_size = {}, count_valid_r_size = {}, latest u = {:x}".format(count_valid_u, count_valid_p_size, count_valid_r_size, u))
            
            if check_pnbits and (p.nbits() < length_p_min or p.nbits() > length_p_max):
                continue
            count_valid_p_size += 1
            
            r = ZZ(rx(ZZ(u)))
            avrg_log2_r += RR(log(RR(r),RR(2.0)))
            if check_rnbits and (r.nbits() < length_r_min or r.nbits() > length_r_max):
                continue
            count_valid_r_size += 1
            
            if check_size_only:
                continue
            
            tw = [ZZ(c_i(ZZ(u))) for c_i in polys_cofact_twists]
            
            i = 0
            
            gcd_pi = gcd(prod_primes, p)
            gcd_ri = gcd(prod_primes, r)
            gcd_p = gcd(gcd_p, gcd_pi)
            gcd_r = gcd(gcd_r, gcd_ri)
            bad_param = gcd_pi > 1 or gcd_ri > 1
            gcd_c_tw = [0 for _ in range(len(polys_cofact_twists))]
            
            j = 0
            for ci in tw:
                small_cofactors_ci, cii = small_factors(ci, prod_primes)
                gcd_ci[j] = gcd(gcd_ci[j], small_cofactors_ci)
                gcd_c_tw[j] = small_cofactors_ci
                if check_all_cofactors:
                    bad_param = bad_param or gcd_c_tw[j] > allowed_cofactor
                j += 1
            bad_param = bad_param or not p.is_pseudoprime() or not r.is_pseudoprime()
            #for ci in tw:
            #    bad_param = bad_param or not ci.is_pseudoprime()
            #bad_param = bad_param or not (tw[0]).is_pseudoprime()

            if not bad_param:
                t = ZZ(tx(ZZ(u)))
                print("\nu = {} mod {}".format(u % m, m))
                print("u = {} # {} bits\np = {} # {} bits {} mod k={}\nr = {} # {} bits".format(u,u.nbits(),p,p.nbits(),p % k, k, r,r.nbits()))
                print("t = {}".format(t))
                if two_naf:
                    str_u, Hw_ = str_2NAF(u)
                    #str_py_u, Hw_ = str_py_2NAF(u)
                else:
                    str_u, Hw_ = str_binary_form(u)
                    #str_py_u, Hw_ = str_py_binary_form(u)
                print("u = {: #x} # {} {}: {}".format(u, str_u, str_w, Hw_))
                print("p = {: #x}".format(p))
                print("r = {: #x}".format(r))
                print("t = {: #x}".format(t))
                
                tw = [ZZ(c_i(ZZ(u))) for c_i in polys_cofact_twists]
                j=0
                for ci in tw:
                    small_cofactors_ci, cii = small_factors(ci, prod_primes)
                    gcd_ci[j] = gcd(gcd_ci[j], small_cofactors_ci)
                    gcd_c_tw[j] = small_cofactors_ci
                    print("{} = {}\n  = ({}) * {}".format(label_factors[j], ci, small_cofactors_ci, cii))
                    print("{}0 = {} = {}".format(label_factors[j], small_cofactors_ci, small_cofactors_ci.factor()))
                    print("{}1 is prime: {}".format(label_factors[j], cii.is_pseudoprime()))
                    j+=1
                if D not in [1, 3, 4]:
                    y = ZZ(yx(ZZ(u)))
                    a, b = compute_a_b(D, p, t, y)
                    E = curve(k=k,u=u,a=a,b=b)
                else:
                    ab=1
                    ok = False
                    while not ok:
                        try:
                            if D == 3:
                                if curve == FotiadisMartindale:
                                    E = curve(code=code,u=u,b=ab)
                                elif curve == Aurifeuillean:
                                    E = curve(k=k, u=u, D=D, e0=e0, aa=a_auri, b=ab)
                                else:
                                    E = curve(k=k,u=u,b=ab)
                            elif D == 1 or D == 4:
                                if curve == FotiadisMartindale:
                                    E = curve(code=code,u=u,a=ab)
                                elif curve == Aurifeuillean:
                                    E = curve(k=k, u=u, D=D, e0=e0, aa=a_auri, a=ab)
                                else:
                                    E = curve(k=k,u=u,a=ab)
                            ok=True
                        except ValueError as err:
                            try:
                                if D == 3:
                                    if curve == FotiadisMartindale:
                                        E = curve(code=code,u=u,b=-ab)
                                    elif curve == Aurifeuillean:
                                        E = curve(k=k, u=u, D=D, e0=e0, aa=a_auri, b=ab)
                                    else:
                                        E = curve(k=k,u=u,b=-ab)
                                elif D == 1 or D == 4:
                                    if curve == FotiadisMartindale:
                                        E = curve(code=code,u=u,a=-ab)
                                    elif curve == Aurifeuillean:
                                        E = curve(k=k, u=u, D=D, e0=e0, aa=a_auri, a=ab)
                                    else:
                                        E = curve(k=k,u=u,a=-ab)
                                ok=True
                                ab = -ab
                            except ValueError as err:
                                ab += 1
                print(E)
                #E.print_parameters_for_RELIC()
                if curve == BLS and k in [12, 24]:
                    lnpk = ceil((k)*RR(log(p,2)))
                    dh_S, stnfs = guess_cost_STNFS_BLS(k, lnpk)
                    dh_C, ctnfs = guess_cost_ConjTNFS_BLS(k, lnpk)
                    stnfs = round(stnfs)
                    ctnfs = round(ctnfs)
                    str_res = "    {{\'u\':{: #x}, \'b\':{:2d}, \'pnbits\':{}, \'rnbits\':{}, \'cost_S\':{}, \'deg_h_S\':{}, \'cost_C\':{}, \'deg_h_C\':{}, \'label\':\"{} {} {}\"}},".format(u,ab,p.nbits(),r.nbits(), stnfs, dh_S, ctnfs, dh_C, str_u, str_w, Hw_)
                else:
                    if D not in [1, 3, 4]:
                        str_res = "    {{\'u\':{: #x}, \'a\':{:2d}, \'b\':{:2d}, \'pnbits\':{}, \'rnbits\':{}, \'cost_S\':None, \'deg_h_S\':None, \'cost_C\':None, \'deg_h_C\':None, \'label\':\"{} {} {}\"}},".format(u, a, b, p.nbits(), r.nbits(), str_u, str_w, Hw_)
                    else:
                        if D == 3:
                            str_ab = "b"
                        elif D==1 or D==4:
                            str_ab = "a"
                        str_res = "    {{\'u\':{: #x}, \'{}\':{:2d}, \'pnbits\':{}, \'rnbits\':{}, \'cost_S\':None, \'deg_h_S\':None, \'cost_C\':None, \'deg_h_C\':None, \'label\':\"{} {} {}\"}},".format(u, str_ab, ab, p.nbits(), r.nbits(), str_u, str_w, Hw_)
                print(str_res)
                ofile.write(str_res+"\n")
                ofile.flush()

    print("count_u = {}\ncount_valid_u = {}".format(count_u,count_valid_u))
    print("count_valid_r_size = {}\n".format(count_valid_r_size))
    print("count_valid_p_size = {}\n".format(count_valid_p_size))

if count_valid_u > 0:
    print("gcd p = {}, gcd r = {}, gcd c = {}, gcd tw = {}, gcd g2c = {}".format(gcd_p, gcd_r, gcd_ci[0], gcd_ci[1], gcd_ci[2]))
    print("avrg log_2 r(u) = {}".format(float(RR(avrg_log2_r/count_valid_u))))
    print("avrg log_2 p(u) = {}".format(float(RR(avrg_log2_p/count_valid_u))))

if not check_size_only:
    ofile.write("]\n")
    ofile.flush()
    ofile.close()
    print(filename)

