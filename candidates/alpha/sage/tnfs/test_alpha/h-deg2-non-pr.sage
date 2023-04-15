ZZy.<y> = ZZ[]

from tnfs.alpha.alpha_tnfs_2d import *
from tnfs.test_alpha.utils import *

h = y**2+5

deg_h = h.degree()
Kh.<ah> = NumberField(h)
Oh = Kh.maximal_order()
Oh1 = Oh.ideal(1) # ideal containing all, for testing wether two ideals are coprime
#OhX.<X> = PolynomialRing(Oh)
Khx.<x> = PolynomialRing(Kh)
Kh_x1x2.<x1,x2> = PolynomialRing(Kh,2)
print("h = {}".format(h))
monic_f = True

deg_f = 4
f = x**4 - 3*ah*x**3 - (6*ah + 1)*x**2 - (ah + 10)*x - 10*ah
ZZX.<X> = ZZy[]
fyX = X**4 - 3*y*X**3 - (6*y + 1)*X**2 - (y + 10)*X - 10*y

disc_f = f.discriminant()
I_disc_f = Oh.ideal(disc_f)
factors_disc_f = I_disc_f.factor()
ld_f = f.leading_coefficient()
I_ld_f = Oh.ideal(ld_f)
factors_ld_f = I_ld_f.factor()
gcd_disc_f = gcd([ZZ(ei) for ei in disc_f.list()])
bare_disc_f = disc_f/gcd_disc_f
I_bare_disc_f = Oh.ideal(bare_disc_f)
factors_bare_disc_f = I_bare_disc_f.factor()

print("f = {}".format(f))
print("disc(f) = {}".format(disc_f))
print("disc_f = {} * ({})".format(gcd_disc_f, disc_f/gcd_disc_f))
#print("Norm(disc_f) = {}".format(ZZ(disc_f.norm()).factor()))
print("Norm(disc_f/{}) = {}".format(gcd_disc_f, ZZ(bare_disc_f.norm()).factor() ))
#print("factors(disc(f)) = {}".format(factors_disc_f))
print("factors(disc(f)/{}) = {}".format(gcd_disc_f, factors_bare_disc_f))
print("ld(f) = {}".format(ld_f))
print("factors(ld(f)) = {}".format(factors_ld_f))

B = 2000
PrB = prime_range(B)
#list_I = [i[0] for l in PrB for i in Oh.ideal(l).factor()]
# since Magma and SageMath do not output the same list of ideals, here is the list given by Magma, in SageMath format
list_I_raw = [[2,[1,1]],[3,[1,1]],[3,[5,1]],[[0,1]],[7,[3,1]],[7,[4,1]],[11],[13],[17],[19],[23,[8,1]],[23,[15,1]],[[3,-2]],[[3,2]],[31],[37],[[-6,-1]],[[-6,1]],[43,[9,1]],[43,[34,1]],[47,[18,1]],[47,[29,1]],[53],[59],[[-4,3]],[[-4,-3]],[67,[14,1]],[67,[53,1]],[71],[73],[79],[83,[24,1]],[83,[59,1]],[[-3,-4]],[[-3,4]],[97],[[-9,2]],[[-9,-2]],[103,[43,1]],[103,[60,1]],[107,[40,1]],[107,[67,1]],[[-8,-3]],[[-8,3]],[113],[127,[54,1]],[127,[73,1]],[131],[137],[139],[[-12,-1]],[[-12,1]],[151],[157],[163,[22,1]],[163,[141,1]],[167,[50,1]],[167,[117,1]],[173],[179],[[-1,6]],[[-1,-6]],[191],[193],[197],[199],[211],[223,[21,1]],[223,[202,1]],[227,[26,1]],[227,[201,1]],[[-7,6]],[[-7,-6]],[233],[239],[[-14,-3]],[[-14,3]],[251],[257],[263,[28,1]],[263,[235,1]],[[-12,-5]],[[-12,5]],[271],[277],[[-6,-7]],[[-6,7]],[283,[109,1]],[283,[174,1]],[293],[307,[84,1]],[307,[223,1]],[311],[313],[317],[331],[337],[347,[79,1]],[347,[268,1]],[[-13,6]],[[-13,-6]],[353],[359],[367,[27,1]],[367,[340,1]],[373],[379],[383,[83,1]],[383,[300,1]],[[-12,7]],[[-12,-7]],[397],[[-9,8]],[[-9,-8]],[[-2,-9]],[[-2,9]],[419],[[-4,-9]],[[-4,9]],[431],[433],[439],[443,[138,1]],[443,[305,1]],[[-18,5]],[[-18,-5]],[457],[[-21,2]],[[-21,-2]],[463,[213,1]],[463,[250,1]],[467,[205,1]],[467,[262,1]],[479],[487,[143,1]],[487,[344,1]],[491],[499],[503,[178,1]],[503,[325,1]],[[-3,-10]],[[-3,10]],[[-21,4]],[[-21,-4]],[523,[97,1]],[523,[426,1]],[[-19,6]],[[-19,-6]],[547,[33,1]],[547,[514,1]],[557],[563,[75,1]],[563,[488,1]],[[-18,7]],[[-18,-7]],[571],[577],[587,[157,1]],[587,[430,1]],[593],[599],[[-14,9]],[[-14,-9]],[607,[128,1]],[607,[479,1]],[613],[617],[619],[631],[[-6,11]],[[-6,-11]],[643,[150,1]],[643,[493,1]],[647,[44,1]],[647,[603,1]],[653],[659],[[-16,9]],[[-16,-9]],[673],[677],[683,[307,1]],[683,[376,1]],[691],[[-24,-5]],[[-24,5]],[[-23,-6]],[[-23,6]],[719],[727,[214,1]],[727,[513,1]],[733],[739],[743,[348,1]],[743,[395,1]],[751],[757],[[-21,-8]],[[-21,8]],[[-7,-12]],[[-7,12]],[773],[787,[119,1]],[787,[668,1]],[797],[[-27,-4]],[[-27,4]],[811],[[-24,-7]],[[-24,7]],[823,[337,1]],[823,[486,1]],[827,[219,1]],[827,[608,1]],[[-28,3]],[[-28,-3]],[839],[853],[857],[859],[863,[274,1]],[863,[589,1]],[877],[[-6,-13]],[[-6,13]],[883,[384,1]],[883,[499,1]],[887,[193,1]],[887,[694,1]],[907,[334,1]],[907,[573,1]],[911],[919],[[-18,-11]],[[-18,11]],[937],[[-21,10]],[[-21,-10]],[947,[330,1]],[947,[617,1]],[953],[967,[295,1]],[967,[672,1]],[971],[977],[983,[133,1]],[983,[850,1]],[991],[997],[[-17,12]],[[-17,-12]],[1013],[1019],[[-29,-6]],[[-29,6]],[1031],[1033],[1039],[[-27,8]],[[-27,-8]],[1051],[[-9,-14]],[[-9,14]],[1063,[383,1]],[1063,[680,1]],[[-32,-3]],[[-32,3]],[1087,[446,1]],[1087,[641,1]],[1091],[1093],[1097],[1103,[105,1]],[1103,[998,1]],[[-33,2]],[[-33,-2]],[1117],[1123,[58,1]],[1123,[1065,1]],[[-2,-15]],[[-2,15]],[1151],[1153],[1163,[221,1]],[1163,[942,1]],[1171],[[-24,-11]],[[-24,11]],[1187,[282,1]],[1187,[905,1]],[1193],[[-34,3]],[[-34,-3]],[1213],[1217],[1223,[424,1]],[1223,[799,1]],[[-27,10]],[[-27,-10]],[1231],[1237],[[-23,-12]],[[-23,12]],[1259],[1277],[1279],[1283,[62,1]],[1283,[1221,1]],[[-3,-16]],[[-3,16]],[1291],[1297],[[-36,-1]],[[-36,1]],[1303,[51,1]],[1303,[1252,1]],[1307,[140,1]],[1307,[1167,1]],[1319],[[-14,-15]],[[-14,15]],[1327,[404,1]],[1327,[923,1]],[[-9,-16]],[[-9,16]],[1367,[64,1]],[1367,[1303,1]],[1373],[[-16,15]],[[-16,-15]],[1399],[[-33,8]],[[-33,-8]],[1423,[196,1]],[1423,[1227,1]],[1427,[458,1]],[1427,[969,1]],[[-32,9]],[[-32,-9]],[1433],[1439],[1447,[611,1]],[1447,[836,1]],[1451],[1453],[1459],[1471],[[-6,17]],[[-6,-17]],[1483,[264,1]],[1483,[1219,1]],[1487,[102,1]],[1487,[1385,1]],[[-38,-3]],[[-38,3]],[1493],[1499],[1511],[1523,[364,1]],[1523,[1159,1]],[1531],[1543,[68,1]],[1543,[1475,1]],[[-37,6]],[[-37,-6]],[1553],[1559],[1567,[564,1]],[1567,[1003,1]],[1571],[1579],[1583,[303,1]],[1583,[1280,1]],[1597],[[-39,-4]],[[-39,4]],[1607,[363,1]],[1607,[1244,1]],[[-22,-15]],[[-22,15]],[1613],[1619],[[-1,18]],[[-1,-18]],[1627,[57,1]],[1627,[1570,1]],[1637],[1657],[1663,[173,1]],[1663,[1490,1]],[1667,[108,1]],[1667,[1559,1]],[[-7,-18]],[[-7,18]],[1693],[1697],[1699],[[-27,-14]],[[-27,14]],[[-21,-16]],[[-21,16]],[1723,[269,1]],[1723,[1454,1]],[1733],[[-11,18]],[[-11,-18]],[1747,[491,1]],[1747,[1256,1]],[1753],[1759],[1777],[1783,[855,1]],[1783,[928,1]],[1787,[679,1]],[1787,[1108,1]],[[-13,18]],[[-13,-18]],[[-26,-15]],[[-26,15]],[1811],[1823,[135,1]],[1823,[1688,1]],[1831],[1847,[547,1]],[1847,[1300,1]],[[-41,-6]],[[-41,6]],[1867,[732,1]],[1867,[1135,1]],[1871],[1873],[1877],[1879],[[-42,-5]],[[-42,5]],[[-36,11]],[[-36,-11]],[1907,[283,1]],[1907,[1624,1]],[1913],[1931],[1933],[[-12,19]],[[-12,-19]],[1951],[1973],[1979],[1987,[63,1]],[1987,[1924,1]],[1993],[1997],[1999] ]

list_I = [Oh.ideal([Kh(Igen) for Igen in Ia]) for Ia in list_I_raw]

print("there are {} prime ideals above the {} prime numbers up to B={}".format(len(list_I), len(PrB), B))

#N = 10**5 # 10^5 = 1 min on muscox
N = 10**8 #N = 10**8 div 32; # catrel
cores = 8 #cores = 128
Ni = N // cores
A = 10**6
print("computing experimentally the valuation for N={} random coprime pairs".format(N))
#print("id = {}, {} samples".format(id, Ni))
#sum_val_I, ratio_coprime_pairs_Oh = MonteCarlo_count_val(Ni,A,Oh,Oh1,deg_h,f,list_I)

ratio_coprime_pairs_Oh = 0.538959693680161748182422735181
sum_val_I = [ 133337292, 0, 24996308, 83336012, 29163957, 14587525, 0, 0, 0, 276846, 4356821, 0, 3334049, 6901968, 0, 0, 9762294, 9760688, 0, 2327637, 0, 4260223, 0, 28847, 3276838, 3279910, 1491800, 1493539, 0, 0, 15994, 2409082, 2409891, 0, 1123108, 21155, 1981221, 3959372, 0, 0, 935132, 0, 919698, 1835347, 15635, 0, 1574813, 0, 10765, 5236, 669525, 0, 4385, 0, 0, 613395, 0, 0, 6578, 6268, 1106025, 1104418, 10772, 2649, 0, 2474, 0, 0, 449068, 0, 0, 873031, 0, 0, 1709, 0, 414985, 0, 0, 1137067, 0, 0, 372079, 1350, 1329, 355967, 711516, 352667, 353894, 2343, 326196, 325580, 0, 0, 0, 1744, 905, 287860, 1153834, 0, 285865, 1616, 0, 271100, 544887, 0, 703, 261418, 261367, 256637, 0, 0, 249921, 498994, 0, 489087, 0, 237410, 0, 1042, 0, 0, 452164, 452181, 0, 445696, 0, 0, 434346, 216041, 0, 0, 0, 456, 0, 204960, 879, 0, 198424, 397880, 0, 393142, 382616, 383568, 0, 191011, 0, 369578, 731049, 0, 649, 356035, 354794, 176350, 351564, 0, 607, 170764, 169605, 258, 267, 167130, 333050, 330303, 329742, 0, 1047, 270, 276, 156253, 311157, 310424, 0, 309149, 309028, 206, 235, 302449, 151328, 262, 0, 292440, 0, 185, 0, 0, 0, 0, 376, 275602, 0, 0, 170, 134813, 134308, 186, 0, 131418, 0, 129358, 0, 0, 254853, 0, 326, 123747, 0, 0, 121280, 244051, 242846, 0, 121036, 0, 241926, 120366, 0, 248, 298, 0, 0, 232335, 0, 0, 0, 225768, 113429, 224969, 112640, 110290, 0, 0, 122, 108010, 107510, 124, 0, 425479, 210868, 211287, 203, 0, 0, 0, 206, 0, 101842, 198, 0, 198380, 98919, 0, 0, 196210, 0, 0, 0, 191, 0, 95029, 183, 0, 94800, 188306, 188194, 0, 93621, 0, 367313, 0, 80, 154, 180337, 90198, 359921, 0, 0, 0, 0, 176702, 88711, 0, 88, 0, 0, 128, 84760, 84843, 83954, 84340, 71, 0, 0, 130, 0, 163868, 0, 0, 162548, 67, 0, 319841, 0, 57, 130, 0, 155672, 77892, 77318, 0, 64, 113, 153113, 0, 76917, 76753, 153143, 306809, 0, 0, 76058, 75380, 75209, 0, 0, 73157, 0, 44, 72611, 0, 125, 142084, 71205, 280843, 70035, 0, 0, 0, 70174, 53, 0, 69057, 69228, 96, 0, 64, 183, 134463, 67366, 0, 67555, 67186, 0, 134547, 0, 183, 48, 73, 0, 65233, 0, 259386, 0, 0, 64661, 81, 0, 127805, 0, 40, 87, 63286, 0, 46, 62337, 0, 61964, 0, 0, 0, 0, 0, 61403, 62022, 0, 122479, 78, 30, 119867, 60119, 0, 60001, 59784, 59992, 81, 67, 0, 0, 58409, 0, 116572, 0, 0, 0, 57613, 0, 57474, 0, 32, 0, 28, 112274, 223383, 112200, 56191, 111582, 0, 55751, 0, 43, 54841, 110023, 57, 108071, 54029, 107466, 53585, 0, 214165, 52, 80, 23, 0, 52653, 0, 52427, 52583, 104555, 104598, 33, 22, 58, 51311, 51362, 17, 54, 0, 0, 100595, 22, 20, 0 ]

print("\nratio_coprime_pairs = {}".format(ratio_coprime_pairs_Oh))
print("\nsum_val_I = {}".format(sum_val_I))
#print("\nratio_coprime_pairs_{} = {}".format(id, ratio_coprime_pairs_Oh))
#print("\nsum_val_I_{} = {}".format(id, sum_val_I))

#quit

thr_val, thr_val_non_Pr, idx_bad_I, idx_proj_I = compute_expected_valuation(list_I,f,I_disc_f,Kh,Oh,Oh1,Kh_x1x2)

av_val = [float(sum_val_I[j]/N) for j in range(len(list_I))]
ratio = ratio_val(av_val, thr_val)
ratio_non_Pr = ratio_val(av_val, thr_val_non_Pr)

wrong_values = print_ratio(ratio, thr_val, av_val, 0.2, N)
print("possibly wrong values at indices {}".format(wrong_values))
print("ideals concerned:")
for idx_wrong in wrong_values:
    I = list_I[idx_wrong]
    Norm_I = abs(ZZ(I.norm()))
    I_pr = I.is_principal()
    if I_pr:
        gen = I.gens_reduced()[0]
	gen_s = "{:11s}".format(str(Kh(gen)))
    else:
        gen1, gen2 = I.gens_two()
	gen_s = "{},{}".format(gen1,gen2)
    print("I{:2d}, Norm {:6d}, <{}>".format(idx_wrong, Norm_I, gen_s))

print("the bad ideals are at index: {}".format(idx_bad_I))

s = ""
for i in idx_bad_I:
    s += "{:.5f}, ".format(float(ratio[i]))
print("the ratio between experimental and theoretical expected valuation for bad ideals is:  ["+s+"]")
print("these bad ideals are principal: {}".format([list_I[i].is_principal() for i in idx_bad_I]))

print("projective roots:")
print("the ideals dividing ld(f) are at index: {}".format(idx_proj_I))
s = ""
for i in idx_proj_I:
    s += "{:.5f}, ".format(float(ratio[i]))
print("the ratio between experimental and theoretical expected valuation for these ideals is:  ["+s+"]")
print("these ideals are principal: {}".format([list_I[i].is_principal() for i in idx_proj_I]))


print("\nwith I.gens_two():")
wrong_values, idx_wrong = print_err_val_non_principal_ideal(thr_val, thr_val_non_Pr)

print("possibly wrong values at indices {}".format(idx_wrong))
print("ideals concerned:")
for idx in idx_wrong:
    I = list_I[idx]
    Norm_I = abs(ZZ(I.norm()))
    I_pr = I.is_principal()
    gen1, gen2 = I.gens_two()
    gen_s_two = "{},{}".format(gen1,gen2)
    if I_pr:
        gen = I.gens_reduced()[0]
	gen_s = "{:11s}".format(str(Kh(gen)))
        print("I{:2d}, Norm {:6d}, <{}> (gens_two <{}>)".format(idx, Norm_I, gen_s, gen_s_two))
    else:
        print("I{:2d}, Norm {:6d}, <{}>".format(idx, Norm_I, gen_s_two))

# there is a problem at index 12
# theoretical valuation computed with "is_principal" is 1/30 = 1/(29+1)
# theoretical valuation computed without principal is 29/840 = 29/(30*28) = (28+1)/(N^2-1)
# = 1/30 + 1/840 = 1/30*(1+1/28)
# so apparently the algorithm counts one more recursive step without principal ideal, why?

l = 29
I = list_I[12]
Norm_I = I.norm()
average_valuation_homogeneous_coprime(f,I_disc_f,I,Norm_I,Kh,Oh,Oh1,Kh_x1x2,test_principal=True)
average_valuation_homogeneous_coprime(f,I_disc_f,I,Norm_I,Kh,Oh,Oh1,Kh_x1x2,test_principal=False)

test_principal=False
Fq = I.residue_field('z')
FqZ = Fq['Z']; (Z,) = FqZ._first_ngens(1)
delta, gam = I.gens_two() # delta=l in ZZ, gam with 'w'

if gam == 0 or (gam.norm() == Norm_I) or (test_principal and I.is_principal()):
    print("I = {} is principal".format(I))
    print("r=13 mod 29, f1 = f(13+gen*x), val_I(f1) = 1, fv1 = f1//gen, FqZ(fv1) = 15 (constant poly)")
    gen = I.gens_reduced()[0]
    fv = f(13+gen*x)
    # fv (264*ah - 599)*x^4 + (-269*ah - 8166)*x^3 + (-10803*ah - 18523)*x^2 + (-22538*ah + 9476)*x - 7628*ah + 28262
    fv.content_ideal().norm().factor()
    fv.content_ideal().valuation(I)
    fv1 = fv // gen
    # (-14*ah - 153)*x^4 + (-591*ah - 752)*x^3 + (-2395*ah + 1809)*x^2 + (-1678*ah + 8752)*x + 1160*ah + 5554
    FqZ(fv1) # this is a constant: 15
    # fv0_FqZ = FqZ([Fq(coeff_i) for coeff_i in fv0.list()])
    # 15
    #end.
else:
    print("I = {} not principal".format(I))
    print("average_valuation_affine_TNFS_I_not_principal(...)")
    Fq_z1z2 = Fq['z1, z2']; (z1, z2,) = Fq_z1z2._first_ngens(2)
    r = average_valuation_affine_TNFS_I_not_principal(f, I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1) * Norm_I
    print("obtained r = {}".format(r))
    r = r/(Norm_I+1)

    # average_valuation_affine_TNFS_I_not_principal(f, I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1)
    L = [(f,0)]
    v = QQ(0)
    while len(L) > 0:
        print("depth = {}".format(L[0][1]))
        v += average_valuation_affine_TNFS_I_not_principal_rec(I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1, L)
        print("v = {}".format(v))
    print("final value v = {}".format(v))

    fv = f(13 + delta*x1 + gam*x2)
    fv.content_ideal().valuation(I) # this is 1
    val_i = min([valuation(Oh.ideal(fi), I) for fi in fv.coefficients()]) # this is also 1
    fv // gam
    fv1 = fv // gam
    lcm_coeffs = lcm([ZZ(fv_ij.denominator()) for fv_ij in flatten([fv_i.list() for fv_i in fv1.coefficients()])]) # this is 6
    co_mult = gcd(ZZ(lcm_coeffs),ZZ(abs(ZZ(Oh.ideal(gam).norm())/Norm_I))**val_i) # this is 6
    fv1 = co_mult*fv1
    dfv1 = fv1.derivative(fv1.variables()[0])
    dfv2 = fv1.derivative(fv1.variables()[1])
    fv_Fq_z1z2 = map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, fv1)
    dfv1_Fq_z1z2 = map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, dfv1)
    dfv2_Fq_z1z2 = map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, dfv2)
    # and here, actually fv_Fq_z1z2 is univariate in z2, and its derivative in z1 equals 0
    roots_bivariate_enum(fv_Fq_z1z2, Fq, FqZ, Oh, Oh1, delta, gam)
    
#alpha_f = alpha_TNFS_2d(fyX,h,B,test_principal=True)
#print("now computes alpha: {}".format(alpha_f))
#print("in basis 2: {:.3f}/log(2) = {:.3f}\n".format(float(alpha_f), float(alpha_f/log(2))))
#print("computes alpha without testing for principal ideals:")
#alpha_f_ = alpha_TNFS_2d(fyX,h,B,test_principal=False)
#print("             alpha: {}".format(alpha_f_))

"""

# nice latex printing
RR = RealField()

print("bad ideals")
for i in idx_bad_I:
    I = list_I[i]
    #d_I = I.degree()
    d_I = I.residue_class_degree()
    gens = I.gens_two()
    l = gens[0]
    l = ZZ(l)
    if len(gens) > 1:
	fi = gens[1]
	fi = ZZy([ZZ(ai) for ai in fi.list()])
	gen_s = "{:3d},{:7s}".format(l,str(fi))
    else:
	gen_s = "{:3d}        ".format(l)
    if d_I > 1:
	print("${:3d}^{}$ & $\\langle {:12s} \\rangle$ & {} & {:18s} & {:.7f} & {:.7f} & {:.7f} \\\\".format(l, d_I, gen_s, I_disc_f.valuation(I), thr_val[i], float(thr_val[i]), float(av_val[i]), float(ratio[i])))
    else:
	print("${:3d}  $ & $\\langle {:12s} \\rangle$ & {} & {:18s} & {:.7f} & {:.7f} & {:.7f} \\\\".format(l, gen_s, I_disc_f.valuation(I), thr_val[i], float(thr_val[i]), float(av_val[i]), float(ratio[i])))

print("regular ideals")
for i in range(100):
    if thr_val[i] == 0 or i in idx_proj_I or i in idx_bad_I:
	continue
    I = list_I[i]
    #d_I = I.degree()
    d_I = I.residue_class_degree()
    #gens = I.generators()
    gens = I.gens_two()
    l = gens[0]
    l = ZZ(l)
    if len(gens) > 1:
	fi = gens[1]
	fi = ZZy([ZZ(ai) for ai in fi.list()])
	gen_s = "{:3d},{:7s}".format(l,str(fi))
    else:
	gen_s = "{:3d}        ".format(l)
    if d_I > 1:
	print("${:3d}^{}$ & $\\langle {:12s} \\rangle$ & {} & {:18s} & {:.7f} & {:.7f} & {:.7f} \\\\".format(l, d_I, gen_s, I_disc_f.valuation(I), thr_val[i], float(thr_val[i]), float(av_val[i]), float(ratio[i])))
    else:
	print("${:3d}  $ & $\\langle {:12s} \\rangle$ & {} & {:18s} & {:.7f} & {:.7f} & {:.7f} \\\\".format(l, gen_s, I_disc_f.valuation(I), thr_val[i], float(thr_val[i]), float(av_val[i]), float(ratio[i])))
"""
