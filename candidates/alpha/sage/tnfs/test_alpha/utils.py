
# compute the valuation
# as in Tab. 3.1 Shi Bai thesis

# compute the average valuation of each prime ideal only by MonteCarlo technique
# N: number of samples
# A: bound, a_ij in [-A,A]

from sage.rings.real_mpfr import RR
from sage.rings.rational_field import Q, QQ
from sage.rings.integer_ring import Z, ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fast_arith import prime_range
from sage.misc.prandom import randint
from sage.functions.log import log
from sage.arith.functions import lcm
from sage.arith.misc import GCD, gcd
from sage.arith.misc import valuation
from sage.misc.flatten import flatten

from tnfs.alpha.alpha_tnfs_2d import *

def MonteCarlo_count_val(N,A,Oh,Oh1,deg_h,f,list_I):
    random_pairs = 0
    coprime_random_pairs = 0
    sum_val_I = [ ZZ(0) for I in list_I ]
    for i in range(1,N+1):
	if (i % 10**5) == 0 :
	    print("{} ".format(i))
        coprime_I = False
        while not coprime_I:
            gcd_1 = False
            while not gcd_1:
	        # generate two coprime ideals in Oh (of any degree)
		aa_coeffs = [randint(-A,A) for j in range(deg_h)]
		bb_coeffs = [randint(-A,A) for j in range(deg_h-1)] + [randint(0,A)]
		random_pairs += 1
	        gcd_1 = (gcd(aa_coeffs + bb_coeffs) == 1)
	    aa = Oh.ideal(Oh(aa_coeffs));  bb = Oh.ideal(Oh(bb_coeffs))
	    coprime_I = ((aa + bb) == Oh1) # ideals are coprime
	coprime_random_pairs += 1
	RI = Oh.ideal(f.resultant(aa+bb*f.variables()[0]))
	for j in range(len(list_I)):
	    sum_val_I[j] += RI.valuation(list_I[j])
    return sum_val_I, 1.0*coprime_random_pairs/random_pairs

def compute_expected_valuation(list_I,f,I_disc_f,Kh,Oh,Oh1,Kh_x1x2, print_first_N_ideals=0):
    idx_bad_I = []
    idx_proj_I = []
    val_f_I = [QQ(0) for i in range(len(list_I))]
    val_f_I_non_Pr = [QQ(0) for i in range(len(list_I))]
    bad_ideal = False
    proj_ideal = False

    for i in range(len(list_I)):
	I = list_I[i]
	Norm_I = abs(ZZ(I.norm()))
	#d_I = I.degree()
        d_I = I.residue_class_degree()
        I_Pr = I.is_principal()
        gen1, gen2 = I.gens_two()
	gens_two_str = "{},{}".format(gen1,gen2)
	if I_Pr :
            genI = I.gens_reduced()[0]
	    gen_s = "{:11s}".format(str(Kh(genI)))
        else:
            gen_s = gens_two_str
	    
	val_f_I[i], bad_ideal, proj_ideal, time_bad = average_valuation_homogeneous_coprime(f,I_disc_f,I,Norm_I,Kh,Oh,Oh1,Kh_x1x2,test_principal=True)
	val_f_I_non_Pr[i], bad_ideal_, proj_ideal_, time_bad_ = average_valuation_homogeneous_coprime(f,I_disc_f,I,Norm_I,Kh,Oh,Oh1,Kh_x1x2,test_principal=False)
	if bad_ideal :
	    idx_bad_I.append(i)
	    print("I{:3d}, Norm {:6d}, <{:16s}>: bad,  val = {}, {} (gens_two = <{}>)".format(i, Norm_I, gen_s, val_f_I[i], val_f_I_non_Pr[i], gens_two_str))
	elif i <= print_first_N_ideals :
	    print("I{:3d}, Norm {:6d}, <{:16s}>: good, val = {}, {} (gens_two = <{}>)".format(i, Norm_I, gen_s, val_f_I[i], val_f_I_non_Pr[i], gens_two_str))
	if proj_ideal :
	    idx_proj_I.append(i)
    return val_f_I, val_f_I_non_Pr, idx_bad_I, idx_proj_I


# exp_val is the average valuation obtained for N>= 10^6 random coprime pairs (a,b)
# if the theoretical valuation is lower than 5/N then it could be that there was no hit while thr_val[i] > 0
def ratio_val(av_val, thr_val):
    size = len(av_val)
    ratio = [float(1.0) for i in range(size)]
    for i in range(size):
	if thr_val[i] != 0 :
	    ratio[i] = av_val[i]/thr_val[i]
	elif av_val[i] == 0 :
	    ratio[i] = 1
	else:
	    ratio[i] = 0
    return ratio

def print_ratio(ratio, thr_val, av_val, approx, N, min_hit_I=50):
    #approx = 0.1
    idx_bad_ratio = []
    for i in range(len(ratio)):
	if (ratio[i] < (1-approx) or ratio[i] > (1+approx)) and (thr_val[i] > min_hit_I/N) :
	    print("idx {:3d} error: ratio {:.5f}, theoretical val {}, average val {}".format(i, float(ratio[i]), float(thr_val[i]), av_val[i]))
	    idx_bad_ratio.append(i)
    return idx_bad_ratio

def print_err_val_non_principal_ideal(I_pr_val, I_non_pr_val):
    idx_err = []
    ratio_err = []
    for i in range(len(I_pr_val)):
        if I_pr_val[i] != I_non_pr_val[i]:
            ratio_err.append(I_pr_val[i] / I_non_pr_val[i])
            idx_err.append(i)
	    print("idx {:3d} error: Ipr val {}, Inpr val {}, ratio pr/npr {}".format(i, I_pr_val[i], I_non_pr_val[i], I_pr_val[i] / I_non_pr_val[i]))
    return ratio_err, idx_err
    
def bad_ideals(list_I, disc_f, Oh1):
    idx_bad = []
    i = 1
    for I in list_I:
	if is_bad_ideal(I, disc_f, Oh1) :
	    idx_bad.append(i)
	i += 1
    return idx_bad

def GCD_ideals_from_coeffs(f, Oh):
    return sum([Oh.ideal(coeff_i) for coeff_i in f.coefficients()])


def get_random_f(deg_f, monic_f, Oh1, Khx, fi_max=12):
    Oh = Order(Oh1)
    Kh = NumberField(Oh)
    deg_h = Degree(Kh)
    ok = False
    if True:
    #while not ok:
        f_irr = False
        while not f_irr:
	    if monic_f :
		f=Khx([Kh([randint(-fi_max,fi_max) for i in range(deg_h)]) for j in range(deg_f)]+[Kh(1)])# monic f
	    else:
                content_1 = False
                while not content_1:
		    f=Khx([Kh([randint(-fi_max,fi_max)for i in range(deg_h)])for j in range(deg_f+1)])
                    # non-monic f, test that Gcd(coeffs) = 1 and ideal<coeff_i> are all coprime
	            gcd_ZZ = gcd([ZZ(fij) for fij in flatten([fi.coefficients() for fi in f.coefficients()])])
		    if gcd_ZZ > 1 :
			f = f // gcd_ZZ
		    content_1 = (sum([Oh.ideal(coeff_i) for coeff_i in f.coefficients()]) == Oh1)
	    f_irr = f.is_irrreducible()
	disc_f = f.discriminant()
	I_disc_f = Oh.ideal(disc_f)
	#factors_disc_f = I_disc_f.factor()
	ld_f = f.leading_coefficient()
	I_ld_f = Oh.ideal(ld_f)
	#factors_ld_f = I_ld_f.factor()
    #ok = (len(factors_disc_f) >= 6) and (len(factors_ld_f) >= 4)
    return f, disc_f, I_disc_f, ld_f, I_ld_f

