import sys
import re
import os
import pprint

from tnfs.simul.cpx import guess_cost_STNFS_BLS, guess_cost_ConjTNFS_BLS
from tnfs.gen.generate_curve_utils import small_factors, str_py_2NAF, str_2NAF, find_min, find_max
from tnfs.gen.generate_curve_utils import get_next_seed
from tnfs.curve.cyclo_kDe import polynomial_params
from tnfs.curve.cyclo_kDe import Cyclo_kDe

#example:
# from top-level of the file tree:
# 1st arg is k, second argument is D, third argument is e0 (in tr = x^e0 + 1 mod r(x))
# 4-th arg is minimum size of p, 5-th (optional) is max size of p
# sage compute_test_vector_Cyclo_kDe.sage 11 1 1 256 512
# sage compute_test_vector_Cyclo_kDe.sage 11 1 6 320 356 # this is equivalent to FST62

# some flags
#CHECK_ALL_COFACTORS = True
CHECK_ALL_COFACTORS = False
CHECK_COFACTORS = False
#CHECK_COFACTORS = True
allowed_cofactor = 2**16

rnbits = 256
#rnbits = 384
#rnbits = 512
step_pnbits = 32
#allow_cofact_r=True ; str_c = "_r_composite"
allow_cofact_r=False ; str_c = "_r_prime"
factor_r = False
#factor_r = True
#negative_u = True;  str_u= "_neg_u_"
negative_u = False;  str_u= "_pos_u_"
FIND_ALL_U = False
#FIND_ALL_U = True
RNBITS = False
#RNBITS = True
#Hw2NAF = 6
Hw2NAF = None

if rnbits >= 256 or rnbits == 1:
    proof.arithmetic(False)

args=sys.argv
print("usage: sage {} <k> <D> <e0> <min_pnbits> <max_pnbits>".format(args[0]))

# read from command-line the embedding degree k
if len(args)>=2:
    k = Integer(args[1])
else:
    raise ValueError("missing embedding degree k")
# read from command-line the discriminant D
if len(args)>=3:
    D = abs(Integer(args[2]))
    if D==4:
        D=1
else:
    raise ValueError("missing discriminant D")
if len(args)>=4:
    e0 = abs(Integer(args[3]))
else:
    raise ValueError("missing exponent 1 <= e0 < k to compute the trace tr = x^e0 +1 mod r(x)")

# read from command-line the size of p
if len(args)>=5:
    start_size=Integer(args[4])
else:
    start_size = 512
if len(args)>=6:
    stop_size=Integer(args[5])
else:
    stop_size=start_size

print("parameters: k={} D={} e0={} min_pnbits={} max_pnbits={}".format(k,D,e0,start_size,stop_size))

# do we need u to be 1 mod 2? YES for at least k=11
px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k,D,e0)
twx = px+1+tx
assert tx^2 - 4*px == -D*yx^2
assert px == (tx^2 + D*yx^2)/4
assert (px+1-tx) == rx*cx

# ideally, factor further cx, and twx. There is at most a denominator 4.
polys_cofact_twists = [cx, twx]
polys_C = []

print("polynomials\npx = {}\nrx = {}\ntx = {}\nyx = {}\ncx = {}".format(px,rx,tx,yx,cx))
# recompute m and u_mod_m
# check if px can take integer values
p_denom = lcm([pi.denom() for pi in px.list()])
if p_denom != 1:
    list_i = []
    list_bad_content = []
    for i in range(0,p_denom):
        pi_denom = lcm([pi.denom() for pi in (px(p_denom*parent(px).0 +i)).list()])
        if pi_denom == 1:
            # then check that the content of px(m*x+i) is 1
            if 1 == gcd([ZZ(pi) for pi in (px(p_denom*parent(px).0+i)).list()]):
                list_i.append(i)
            else:
                list_bad_content.append(i)
    if len(list_i) > 0:
        if p_denom == 4 and list_i == [1,3]:
            m = 2
            u_mod_m = 1
        elif p_denom == 4 and list_i == [0,2]:
            m = 2
            u_mod_m = 0
        elif p_denom == 8 and list_i == [0,4]:
            m = 4
            u_mod_m = 0
        elif p_denom == 8 and list_i == [1,3,5,7]:
            m = 2
            u_mod_m = 1
        else:
            u_mod_m = list_i[0]
            #u_mod_m = list_i[-1]
            m = p_denom
        print("px_denom = {}, u mod {} should be in {}".format(p_denom, p_denom, list_i))
    else:
        print("k={},D={},e0={}, px=({})/{} but no value of u mod {} found".format(k,D,e0,p_denom*px,p_denom,p_denom))
        print("bad content for i in {}".format(list_bad_content))
else:
    m=None
    u_mod_m = None
if k==8 and D==1 and e0==1:
    m=8
    list_i=[0]
    u_mod_m=0
if k==9 and D==3 and e0==5:
    m=6
    list_i=[2]
    u_mod_m=2 # and in this case, then r= r/3 --> modified in Cyclo_kDe, and now b=1
if k==10 and D==5:
    # [0,4,6,10,14,16] mod 20 -> [0,4,6] mod 10
    m = 10
    list_i = [0,4,6]
    u_mod_m = 0
elif k==10 and D==15:
    m=15
    list_i=[1,3,6,13]
    u_mod_m = 1
elif k==10 and D==2 and e0 == 9:
    m=4
    list_i = [0]
    u_mod_m = 0
elif k==12 and D==2:
    m=2
    list_i = [1]
    u_mod_m = 1
elif k==12 and D==3 and e0 == 19:
    m=30
    list_i = [8]
    u_mod_m = 8
elif k==12 and D==3 and e0 == 20:
    m=285
    list_i = [209,266]
    u_mod_m = 209

if Hw2NAF != None:
    strHw = "_Hw{}".format(Hw2NAF)
else:
    strHw=""
if m != None and u_mod_m != None:
    filename = "test_vector_Cyclo_k{}_D{}_e{}_pnbits_{}_{}_rnbits_min_{}{}{}{}_mod_{}{}.py".format(k,D,e0,start_size,stop_size,rnbits,str_c,str_u,u_mod_m,m,strHw)
else:
    filename = "test_vector_Cyclo_k{}_D{}_e{}_pnbits_{}_{}_rnbits_min_{}{}{}_all_ZZ{}.py".format(k,D,e0,start_size,stop_size,rnbits,str_c,str_u,strHw)
    print("output file is " +filename )
    m=1
    u_mod_m = 0
print("output file is " +filename )
ofile = open(filename, "a+")

print("start_size = {} bits, stop_size = {} bits (pnbits)".format(start_size, stop_size))
pnbits=start_size
u_max_all = find_max(px, stop_size, m=m, umodm=u_mod_m)
u_min_all = find_min(px, start_size, m=m, umodm=u_mod_m)
# if polys are even then negative u is useless

# for D=3, polys are not even, searching for negative u is worth doing it
if negative_u: # u_min < u_max < 0
    u_max = max([u0 for u0 in u_min_all if u0 < 0])
    u_min = min([u0 for u0 in u_max_all if u0 < 0])
    print("u_max = {}".format(u_max))
    print("u_min = {}".format(u_min))
    u_step = -m
    if FIND_ALL_U:
        u_start = u_max
        u_stop = u_min
else: # 0 < u_min < u_max
    u_max = max([u0 for u0 in u_max_all if u0 > 0])
    u_min = min([u0 for u0 in u_min_all if u0 > 0])
    print("u_max = {}".format(u_max))
    print("u_min = {}".format(u_min))
    u_step = m
    if FIND_ALL_U:
        u_start = u_min
        u_stop = u_max
# note: we always should have u_min < u_max

# give min and max u for r of given bits
if rnbits >= 256:
    u_min_all_r = find_min(rx, rnbits, m=m, umodm=u_mod_m)
    u_max_all_r = find_max(rx, rnbits, m=m, umodm=u_mod_m)
    u_max_r_neg = max([u0 for u0 in u_min_all_r if u0 < 0])
    u_min_r_neg = min([u0 for u0 in u_max_all_r if u0 < 0])
    u_max_r_pos = max([u0 for u0 in u_max_all_r if u0 > 0])
    u_min_r_pos = min([u0 for u0 in u_min_all_r if u0 > 0])
    
    print("for r being of {} bits".format(rnbits))
    print("u_min = {} pnbits = {}\nu_max = {} pnbits = {}".format(u_min_r_pos, ZZ(px(u_min_r_pos)).nbits(), u_max_r_pos, ZZ(px(u_max_r_pos)).nbits()))
    print("u_min = {} pnbits = {}\nu_max = {} pnbits = {}".format(u_min_r_neg, ZZ(px(u_min_r_neg)).nbits(), u_max_r_neg, ZZ(px(u_max_r_neg)).nbits()))
    if RNBITS:
        if negative_u:
            u_max = u_max_r_neg
            u_min = u_min_r_neg
        else:
            u_max = u_max_r_pos
            u_min = u_min_r_pos
    
ofile.write("test_vector_Cyclo_k{}_D{}_e{}_{}_{} = [\n".format(k,D,e0,start_size,stop_size))
ctr = 0
low = True
while (FIND_ALL_U and (u_max > u_min)) or (not FIND_ALL_U and (pnbits <= stop_size)):
    if not FIND_ALL_U:
        if low:
            # compute start_u, stop_u for given pnbits
            if RNBITS:
                u_max_pnbits = find_max(rx, rnbits, m=m, umodm=u_mod_m)
                u_min_pnbits = find_min(rx, rnbits, m=m, umodm=u_mod_m)
            else:
                u_max_pnbits = find_max(px, pnbits, m=m, umodm=u_mod_m)
                u_min_pnbits = find_min(px, pnbits, m=m, umodm=u_mod_m)
            if negative_u: # u_stop < u_start < 0
                u_start = max([u0 for u0 in u_min_pnbits if u0 < 0])
                u_stop = min([u0 for u0 in u_max_pnbits if u0 < 0])
                u_step = -m
            else: # 0 < u_start < u_stop
                u_start = min([u0 for u0 in u_min_pnbits if u0 > 0])
                u_stop = max([u0 for u0 in u_max_pnbits if u0 > 0])
                u_step = m

            res = get_next_seed(u_start, u_stop, m, [u_mod_m], px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_C,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
        else:
            res = get_next_seed(u_stop, u_start, -m, [u_mod_m], px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_C,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
        low = not low # flip
    else:
        res = get_next_seed(u_start, u_stop, u_step, px, rx, rnbits, allow_cofact_r=allow_cofact_r, polys_cofact_twist=polys_C,verbose=True,allowed_cofactor=allowed_cofactor,factor_r=factor_r)
    u, p, r, cofact_tw_orders, tries = res
    if u == None:
        if not FIND_ALL_U:
            print("nothing found for size {} and u={} mod {}, u_start = {}, u_stop = {}, u_step = {}".format(pnbits,u_mod_m,m, u_start, u_stop, u_step))
            low = True
            pnbits += step_pnbits
        else:
            print("nothing found for size {} and u={} mod {}, u_start = {}, u_stop = {}, u_step = {}".format(pnbits,u_mod_m,m, u_start, u_stop, u_step))
            u_max = u_min-abs(u_step) # this will stop the loop
        continue
    print("u_min   = {}".format(u_min))
    print("u_max   = {}".format(u_max))
    print("u_start = {}".format(u_start))
    print("u_stop  = {}".format(u_stop))
    print("u_step  = {}".format(u_step))
    if allow_cofact_r:
        cofact_r = cofact_tw_orders[0]
    else:
        cofact_r = 1
    print("u = {} #{} bits\np = {} #{} bits\nr = {} #{} bits\ncofact_r = {} #{} bits".format(u,u.nbits(),p,p.nbits(), r, r.nbits(), cofact_r, cofact_r.nbits()))
    # recompute cofactors
    cofact_tw_orders = [ZZ(c_i(u)) for c_i in polys_cofact_twists]
    c,tw = cofact_tw_orders
    print("c   = {}".format(c))
    print("tw  = {}".format(tw))

    print("u ={: #x}".format(u))
    print("tries = {}".format(tries))
    ab=1
    if D==3 or D==1 or D==4:
        ok = False
    else:
        ok = True
        # formula to compute a and b from j is
        # a = -3*x/(x-1728) ; b = 2*x/(x-1728)
        # 1728 * 4*a^3/(4*a^3+27*b^2) == x
        if D == 5 or D==15:
            # j-invariant is a root of HilbertClassPolynomial(-20) = x^2 - 1264000*x - 681472000 of discriminant 2^18 * 5^3 * 13^2 * 17^2 = 5*(2^9*5*13*17)^2
            # it means to compute a root j0, we need to have sqrt(5).
            # modulo p(x), we know sqrt(-5), and moreover, px = 1 mod 4 for u=0,4,6,10,14,16 mod 20 so there is sqrt(5) mod p.
            Fp = GF(p)
            Fpz.<z> = Fp[]
            if D==5: # hilbert_class_polynomial(-D)
                HD = z^2 - 1264000*z - 681472000
            elif D==15:
                HD = z^2 + 191025*z - 121287375
            JJ = HD.roots()
            if len(JJ) == 0:
                raise ValueError("Error no j-invariant found in GF(p)")
            j0 = JJ[0][0]
            ap = -3*j0/(j0-1728) ; bp = 2*j0/(j0-1728)
            print("found ap and bp")
            if (j0/(j0-1728)).is_square():
                omega2 = sqrt(j0/(j0-1728))
                ap = Fp(-3); bp = bp/omega2^3
                a = -3; b = ZZ(bp)
                print("a=-3")
            elif ap.is_square():
                omega2 = sqrt(ap)
                ap = Fp(1); bp = bp/omega2^3
                a = 1; b = ZZ(bp)
                print("a=1")
            else:
                w0 = ZZ(-1)
                w0p = Fp(-1)
                while w0.is_square() or w0p.is_square():
                    w0 += 1
                    w0p += 1
                omega2 = sqrt(ap/w0p) # the product of two non-square is a square mod p
                # a = a/omega2^2 = ap/(ap/w0p) = w0p
                ap = ap/omega2^2
                bp = bp/omega2^3
                a = w0; assert Fp(a) == ap, "Error a mod p != ap"
                b = ZZ(bp)                
            #else:
            #    a = ZZ(ap)
            #    b = ZZ(bp)
        if D==2: # j=8000
            a=-30;b=56
        if D == 11:
            a=-264; b=1694 #(j-invariant is -32768 = -2^15)
        if D==2 or (D==5) or (D==11) or (D==15):
            try:
                print("constructing E, \na={}\nb={}".format(a,b))
                E = Cyclo_kDe(k=k,D=D,e0=e0,u=u,a=a,b=b, cofactor_r=cofact_r)
            except ValueError as err:
                print(err)
                # usually the problem comes from E being actually the quadratic twist of order (p+1+t) instead of the curve itself of order (p+1-t)
                # the quadratic twist has equation (w^4*a, w^6*b) where w is neither a square nor a cube in Fp
                Fp = GF(p)
                w0 = ZZ(-1)
                w0p = Fp(-1)
                while w0.is_square() or w0p.is_square():
                    w0 += 1
                    w0p += 1
                print("trying with a={}^2*a, b={}^3*b".format(w0,w0))
                a = a*w0^2
                b = b*w0^3
                E = Cyclo_kDe(k=k,D=D,e0=e0,u=u,a=a,b=b, cofactor_r=cofact_r)
                
    print("cofact_r={}".format(cofact_r))
    ab=1
    while not ok:
        try:
            if D == 1 or D == 4:
                E = Cyclo_kDe(k=k,D=1,e0=e0,u=u,a=ab, cofactor_r=cofact_r)
            elif D == 3:
                E = Cyclo_kDe(k=k,D=3,e0=e0,u=u,b=ab, cofactor_r=cofact_r)
            ok=True
        except ValueError as err:
            #print(err)
            try:
                ab = -ab
                if D == 1 or D == 4:
                    E = Cyclo_kDe(k=k,D=1,e0=e0,u=u,a=ab, cofactor_r=cofact_r)
                elif D == 3:
                    E = Cyclo_kDe(k=k,D=3,e0=e0,u=u,b=ab, cofactor_r=cofact_r)
                    ok=True
            except ValueError as err:
                ab = -ab + 1
    if D==3 or D==1 or D==4 or D==2 or D==11 or D==5 or D==15:
        print(E)
        if E.p() != p:
            print("strange, p is different\np={}\np={}# from E.p()".format(p,E.p()))
    lnpk = ceil((k)*RR(log(p,2)))
    dh_S, cost_S, dh_C, cost_C = k, 128, k, 200
    if rnbits == 384:
        cost_S = 192
    str_u, Hw_2NAF_u = str_2NAF(u)
    if allow_cofact_r:
        cofact_r_str = "\'cofact_r\':{:7d}, ".format(cofact_r)
    else:
        cofact_r_str = ""
    if FIND_ALL_U:
        label_u = "\'label\':\"u={:34s} Hw2NAF={}\"".format(str_u, Hw_2NAF_u)
    else:
        label_u = "\'label\':\"Cyclo_k{}_D{}_e{}_p{}\"".format(k,D,e0, p.nbits())
    if D==3:
        ofile.write("    {{\'u\':{:#7x}, \'b\':{:2d}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}\'deg_h_S\':{},\'cost_S\':{:3}, \'deg_h_C\':{:},\'cost_C\':{:3}, {}}},#\n".format(u,ab,p.nbits(),r.nbits(), cofact_r_str, dh_S, cost_S, dh_C, cost_C, label_u, ctr))
    elif D==1 or D==4:
        ofile.write("    {{\'u\':{:#7x}, \'a\':{:2d}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}\'deg_h_S\':{},\'cost_S\':{:3}, \'deg_h_C\':{:},\'cost_C\':{:3}, {}}},#\n".format(u,ab,p.nbits(),r.nbits(), cofact_r_str, dh_S, cost_S, dh_C, cost_C, label_u, ctr))
    elif D==2 or D==11 or D==5 or D==15:
        ofile.write("    {{\'u\':{:#7x}, \'a\':{:2d},\'b\':{:2d}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}\'deg_h_S\':{},\'cost_S\':{:3}, \'deg_h_C\':{:},\'cost_C\':{:3}, {}}},#\n".format(u,a,b,p.nbits(),r.nbits(), cofact_r_str, dh_S, cost_S, dh_C, cost_C, label_u, ctr))
    else:
        ofile.write("    {{\'u\':{:#7x}, \'pnbits\':{:3d},\'rnbits\':{:3d}, {}\'deg_h_S\':{},\'cost_S\':{:3}, \'deg_h_C\':{:},\'cost_C\':{:3}, {}}},#\n".format(u,p.nbits(),r.nbits(), cofact_r_str, dh_S, cost_S, dh_C, cost_C, label_u, ctr))
    
    ofile.flush()
    # find next value of u_max
    if not FIND_ALL_U:
        if low:
            pnbits += step_pnbits
    else:
        u_start = u + u_step
        if negative_u:
            u_max = u_start
        else:
            u_min = u_start
    ctr+=1
ofile.write("]\n")
ofile.flush()
ofile.close()
print("output file is " +filename )
