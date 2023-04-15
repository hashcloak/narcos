# Polynomial Selection, computation of alpha, and Murphy _E_ value for STNFS

How to cite this work:

   **On the alpha value of polynomials in the tower number field sieve algorithm**,
    _Aurore Guillevic and Shashank Singh_,
	Mathematical Cryptology, Vol 1 No 1 (Feb 2021),
	[open access](https://journals.flvc.org/mathcryptology/issue/view/5910),
    eprint [2019/885](https://eprint.iacr.org/2019/885), 2019.

A complete example is available in `example_bn_254.py`, in SageMath, run
```python
sage: load("example_bn_254.py")
```
or
```python
sage: run example_bn_254
```
or in command-line:
```shell
sage example_bn_254.py
```

Related papers:
 -  **A short-list of pairing-friendly curves resistant to Special TNFS at the 128-bit security level**,
    _Aurore Guillevic_,
    in A. Kiayias and M. Kohlweiss and P. Wallden and V. Zikas, eds,
    [PKC 2020](https://pkc.iacr.org/2020/program.php#session-189),
    LNCS 12111, pp. 535-564, Springer.
    DOI [10.1007/978-3-030-45388-6_19](https://doi.org/10.1007/978-3-030-45388-6_19),
    eprint [2019/1371](https://eprint.iacr.org/2019/1371).  
    The curves in this paper are given in the file
    `example_curves_short_list.sage` which can be loaded from sage:
```python
load("example_curves_short_list.sage")
```
-   **Optimized and secure pairing-friendly elliptic curves suitable for one layer proof composition**,
    _Youssef El Housni and Aurore Guillevic_,
    in Stephan Krenn and Haya Shulman and Serge Vaudenay, eds,
    [CANS 2020](https://cans2020.at/program-day-2/),
    LNCS 12579, pp. 259-279, Springer.
    DOI [10.1007/978-3-030-65411-5_13](https://doi.org/10.1007/978-3-030-65411-5_13),
    eprint [2020/351](https://eprint.iacr.org/2020/351).  
    This paper uses this code to estimate the cost of a DL computation in the
    target field of a curve **BW6_761**.
    An example with this curve is available in
    `example_bw6_761.py`, run one of the following, from shell or within SageMath:
```shell
sage example_bw6_761.py
```
```python
sage: load("example_bw6_761.py")
```
 -  **Curves with fast computations in the first pairing group**
    _Rémi Clarisse and Sylvain Duquesne and Olivier Sanders_,
    in Stephan Krenn and Haya Shulman and Serge Vaudenay, eds,
    [CANS 2020](https://cans2020.at/program-day-2/),
    LNCS 12579, pp. 280--298, Springer.
    DOI [10.1007/978-3-030-65411-5_14](https://doi.org/10.1007/978-3-030-65411-5_14),
    eprint [2020/760](https://eprint.iacr.org/2020/760)   
    The authors focus on pairing-friendly curves of prime embedding degrees
    11, 13, 17, 19. They obtain the smallest possible **G**<sub>1</sub> at the
    128-bit security level (286 bits).
    They used this code to estimate the DL security of their curves.

Prequel paper:
 -  **Cocks-Pinch curves of embedding degrees five to eight and optimal ate pairing computation**,
    _Aurore Guillevic and Simon Masson and Emmanuel Thomé_,
    Designs, Codes and Cryptography, **88**(6), pp. 1047-1081, 2020.
    DOI [10.1007/s10623-020-00727-w](https://doi.org/10.1007/s10623-020-00727-w)
    eprint [2019/431](https://eprint.iacr.org/2019/431)

Four Cocks-Pinch curves are presented in this paper. The data to simulate NFS
and TNFS and obtain an estimate of the security level is given in the files
[`example_cocks_pinch_nfs.py`](./example_cocks_pinch_nfs.py) and
[`example_cocks_pinch_tnfs.py`](./example_cocks_pinch_tnfs.py).
The procedure is to run the simulation for each possible polynomial selection
and let the code optimise the parameters untils it obtains a balanced cost
between relation collection and linear algebra.
For NFS, there are the following polynomial selections:
-   Conjugation [BGGM15 at Eurocrypt'15](https://eprint.iacr.org/2016/605)
-   Sarkar--Singh [Algorithm A](https://eprint.iacr.org/2015/944)
-   Joux--Lercier [Math Comp'2003 vol.**72**(242)](https://www.ams.org/journals/mcom/2003-72-242/S0025-5718-02-01482-5/S0025-5718-02-01482-5.pdf)
-   generalized Joux--Lercier [BGGM15](https://eprint.iacr.org/2016/605)
-   JLSV1 [Joux-Lercier-Smart-Vercauteren CRYPTO'06](https://www.iacr.org/archive/crypto2006/41170323/41170323.pdf)
-   Special, when the characteristic _p_ has a special form
-   Base-_m_ (folklore, obsolete)

For TNFS, not all possibilities are implemented. What is not implemented is when
gcd(deg(_h_), k/deg(_h_)) > 1, i.e. the cases addressed in
[Kim-Jeong PKC'2017](https://eprint.iacr.org/2016/526).
The available ones are:
-   Conjugation, all cases
-   generalized Joux--Lercier, when GCD is 1
-   Special, when the characteristic _p_ has a special form
-   Sarkar--Singh, when GCD is 1


# Organisation of code (from January 6th, 2020)
The Sage code is now reorganised, a package `tnfs` contains the Python code,
with subdirectories
 -   `alpha/`: the implementation of alpha for TNFS;
 -   `test_alpha/`: tests of alpha and precomputed valuations of ideals;
 -   `poly_h/`: all the precomputed irreducible momic polynomials _h_;
 -   `curve/`: classes defining families of curves, such as BN, BLS12, KSS16. To
     init a curve, one need to instantiate a class with a seed. The class
     `Cyclo_kDe` is a more generic class that requires the embedding degree _k_,
     the CM discriminant _D_ and the exponent when computing the trace with the
     Brezing--Weng method.
 -   `param/`: precomputed test vectors of seeds of curves.
 -   `gen/`: code for generating sparse seeds for BLS and KSS curves (more to
     come).
 -   `simul/`: code for polynomial selection and simulation of NFS and TNFS.


# Compatibility with Python 2.7 (until SageMath 8.) and Python 3.7 (from SageMath 9.)
The code was first developped and tested with Python 2.7. Since July 2020, it
has been changed to be compatible with both Python 2.7 and Python 3.7. The
adaptations are the following:
-   change tabulations to spaces in indentation
-   for dictionary structures, change `dic.has_key(k)` to `k in dic`
-   for print function, always use parenthesis `print("{}".format(k))`
-   replace `<>` with `!=`
-   absolute paths with import statements (e.g. `from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d`)
-   about `import ValueError`: if Python version is 3 or higher, it is
    automatically imported. Otherwise, import it.
```python
from sys import version_info
if sys.version_info[0] < 3:
    from exceptions import ValueError
```
-   `None` no longer is mapped to `0` by default, a comparison `if x != 1` when
    `x` is `None` will fail with Python 3, change to `if x is not None and x != 1`
    when appropriate.

# Polynomial selection

The files `tnfs/simul/polyselect.py`, `tnfs/simul/polyselect_utils.py` and
`tnfs/simul/polyselect_h.py` contain code to compute polynomials for TNFS. Not
all combinations are available yet.

## Polynomial selection for _h_ (base of the tower)
The polynomial _h_ is univariate, monic and irreducible modulo _p_, where _p_ is
the characteristic of the finite field.
Tables of polynomials _h_ are precomputed for degrees 2, 3, 4, 5, 6, 7, 8, 9,
10, 11, 12, 13, 14, 16. For other degrees, the user needs to generate a table of
polynomials _h_.
The tables for degrees 2 to 9 were generated with the command
```python
from tnfs.simul.polyselect_utils import get_list_irr_poly
for d in range(2,10):
    for c in range(1,max(1,7-d)+1):
        get_list_irr_poly(d, max_coeff=c, compute_zeta=True, output_file="tnfs/poly_h/tab_h_{}_{}.py".format(d,c))
```

Up to degree, say, 15, it is feasible to enumerate all polynomials _h_ of given
degree and coefficients in {-1,0,1}, but there are too many such polynomials
already for degree 10.
Run the function
```python
tab_h_8 = get_list_irr_poly(deg=8)
```
To save the table directly to a file, without printing it and without storing
it (for large degrees), run
```python
get_list_irr_poly(deg=10, output_file="tab_h_10.py", onthefly=True, compute_zeta=True)
```
For large degrees (from 11), there are many irreducible polynomials, and this is
enough to restrict the search to sparse polynomials. For this, run
```python
from tnfs.simul.polyselect_h import get_list_irr_polys_nonzero_coeffs
for s in range(2,7):
    get_list_irr_polys_nonzero_coeffs(12, s, compute_zeta=True, output_file="tnfs/poly_h/tab_h_12_1_sparse_{}.py".format(s))
```

For degrees larger than 16, the exact computation of zeta is not satisfying.
For degrees 18 and 24, the following code was run.    
To generate only the polynomials, without computing zeta:
```python
from tnfs.simul.polyselect_h import get_list_irr_polys_nonzero_coeffs

get_list_irr_polys_nonzero_coeffs(deg=24, nonzero_coeffs=2, compute_zeta=False, output_file="tnfs/poly_h/tab_h_24_1_sparse_2_no_zeta.py", onthefly=True, verbose=False)
# results written in file tnfs/poly_h/tab_h_24_1_sparse_2_no_zeta.py
# There were 0 irreducible polynomials

get_list_irr_polys_nonzero_coeffs(deg=24, nonzero_coeffs=3, compute_zeta=False, output_file="tnfs/poly_h/tab_h_24_1_sparse_3_no_zeta.py", onthefly=True, verbose=False)
# results written in file tnfs/poly_h/tab_h_24_1_sparse_3_no_zeta.py
# There were 33 irreducible polynomials

get_list_irr_polys_nonzero_coeffs(deg=24, nonzero_coeffs=4, compute_zeta=False, output_file="tnfs/poly_h/tab_h_24_1_sparse_4_no_zeta.py", onthefly=True, verbose=False)
# results written in file tnfs/poly_h/tab_h_24_1_sparse_4_no_zeta.py
# There were 168 irreducible polynomials

get_list_irr_polys_nonzero_coeffs(deg=24, nonzero_coeffs=5, compute_zeta=False, output_file="tnfs/poly_h/tab_h_24_1_sparse_5_no_zeta.py", onthefly=True, verbose=False)
# results written in file tnfs/poly_h/tab_h_24_1_sparse_5_no_zeta.py
# There were 6878 irreducible polynomials

get_list_irr_polys_nonzero_coeffs(deg=24, nonzero_coeffs=6, compute_zeta=False, output_file="tnfs/poly_h/tab_h_24_1_sparse_6_no_zeta.py", onthefly=True, verbose=False)
# results written in file tnfs/poly_h/tab_h_24_1_sparse_6_no_zeta.py
# There were 29538 irreducible polynomials
```
To obtain about 100 irreducible degree _d_ polynomials _h_ modulo a given prime _p_, one
needs roughly 100_d_ irreducible polynomials over **Q** (an irreducible
polynomial of degree _d_ has roughly a probability 1/_d_ to be irreducible
modulo a given prime _p_).

Polynomials of degree 24 and 3 non-zero coefficients were generated with
the following:
```python
from tnfs.simul.polyselect_h import get_list_irr_polys_nonzero_coeffs
get_list_irr_polys_nonzero_coeffs(deg=24, nonzero_coeffs=3, compute_zeta=True, inv_zeta_Kh_min=0.0, simul_zeta=True, N_simul_zeta=10**4, output_file="tnfs/poly_h/tab_h_24_1_sparse_3_prec_10e4_not_sorted.py", onthefly=True, verbose=False)
```
then sorted with `sort -r` t finally obtain `tnfs/poly_h/tab_h_24_1_sparse_3_prec_10e4.py`.
There is no irreducible polynomial of degree 24 and at most two non-zero
coefficients. All the irreducible polynomials of degree 24 and 4 non-zero
coefficients have value of inverse of zeta smaller than 0.75.


## Computation of zeta
The value of 1/zeta(_Kh_, 2) is required for each polynomial _h_, where _Kh_ is
the number field over **Q** defined by _h_. It can be computed with
[PARI](http://pari.math.u-bordeaux.fr/), or directly from SageMath.
```python
#http://pari.math.u-bordeaux.fr/dochtml/html-stable/_L_minusfunctions.html
#pari.allocatemem(12702240768)
pari.default("realbitprecision")
# usually this is 53 bits of precision by default. To change it:
pari.default("realbitprecision",80)
```

## Polynomial selection for (_f_,_g_) (top of the tower)

The class `Polyselect` defined in `tnfs/simul/polyselect.py` contains methods to
compute pairs of polynomials for
Conjugation, Joux-Lercier, generalized Joux-Lercier, Sarkar-Singh, Special
(auxiliary polynomials should be provided).
An initialized instance of a class BN, BLS12, KSS16, KSS18 or BLS24 (or any
other user-defined class of the same format) is required.
Example:
```python
import tnfs
from tnfs.curve.bn import BN
from tnfs.simul.polyselect import Polyselect
# get the polynomial p(x) defining the characteristic p for BN curves
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bn.polynomial_params()

# create curve
seed = -(2**62+2**55+1) # see testvector_sparseseed.py for other values
E = BN(u=seed, b=2)
print(E)
E.print_parameters()
# create an instance of Polyslect class
poly_init = Polyselect(E, deg_h=6)
poly_init.compute_h(deg_h=6)
list_h = poly_init.get_h()[6]
for hi in list_h:
    print hi
inv_zeta_Kh, w, hc = list_h[0]
ZZy = ZZ['y']; (y,) = ZZy._first_ngens(1)
h = ZZy(hc)
# Special Joux-Pierrot construction
# since deg_h = 6, then deg_g = 12/6 = 2
# since gcd(deg_h, deg_g) = 2 > 1, g should have algebraic coefficients: with_y=True
f, g, max_fij, max_gij, aut_fg = poly_init.TNFS_Special(deg_g=2, h=h, poly_p=px, u=seed, with_y=True)
print(f)
print(g)
print(h)
```

## Other polynomial selections from cado-nfs
To generate two quadratic polynomials with Montgomery's method, one can use the
following binary code from cado-nfs (first run `make twoquadratics` as it is not
included in the default building set, then `cd` to the build directory, by
default it is `cado-nfs/build/${HOSTNAME}/` where `$HOSTNAME` is the machine
hostname).
```shell
./polyselect/twoquadratics -N <p> -maxP 1000000 -skewness 1.0 -q -Bf 2.16e22 -Bg 2.16e22 -area 3.40e38
```
where `<p>` should be replaced by an integer (a prime number _p_ for discrete
logarithm).    

A generalization of Montgomery's method produces two cubic polynomials, but
their coefficient size is no longer optimal. To run this method:

```shell
cd cado-nfs
make twocubics
cd build/`hostname`
./polyselect/twocubics -N <p> -P 1000000 -skewness 1.0 -q -Bf 2.16e22 -Bg 2.16e22 -area 3.40e38
```

To generate polynomials with Joux-Lercier polynomial selection, one can run
```script
./polyselect/dlpolyselect -N <p> -df 3 -dg 2 -bound 4 -t 2 -Bf 2.16e22 -Bg 2.16e22 -area 3.40e38
```

# Simulation

## Simulation of STNFS: example for a BN curve

With the above piece of code already executed, now use
```python
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d
from tnfs.simul.polyselect_utils import automorphism_factor

# parameters
p = E.p()
ell = E.r()
Fp = E.Fp()
Fpz,z = E.Fpz()
Rxy = f.parent()

cost = 102 # see testvector_sparseseed.py

deg_h = 6
deg_g = 2
if (gcd(deg_g,deg_h) == 1):
    aut_h = automorphism_factor(h_coeffs)
else:
    aut_h = 1

aut = aut_h*aut_fg

# computing alpha, this takes at least a few seconds

alpha_f = float(alpha_TNFS_2d(f,h,1000))
alpha_g = float(alpha_TNFS_2d(g,h,1000))
sum_alpha = alpha_f+alpha_g
print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f,alpha_g,sum_alpha))

#initialisation of data
simul = Simulation_TNFS(p,ell,Fp,Fpz,h,f,g,Rxy,cost,aut,inv_zeta_Kh,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

simul.print_params()
simul.simulation(samples=100000) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
simul.print_results()

```

# Generation of parameters (seeds) for pairing-friendly curves
The file `compute_test_vector_cyclo_kDe.sage` can be run from a shell to compute
an array of seeds so that _p_ has some given size. It contains parameters, such
as checking for co-factors, having the need negative or positive. Example:
```shell
sage compute_test_vector_cyclo_kDe.sage 12 3 1 384 512
```
will generate BLS12 parameters where _p_ is 384 to 512 bit long by steps of 32
bits, and store the data in
`test_vector_Cyclo_k12_D3_e1_pnbits_384_512_rnbits_min_256_r_prime_pos_u_1_mod_3.py`.
In addition, the seed _u_ will be positive (switch the flag `negative_u` to True to
change it), and _r_ will be prime.

# Generation of sparse seeds
The files in `tnfs/gen` search for sparse seeds for BLS and KSS curves and
compute the cofactor and the order of the quadratic twist.

To generate sparse seeds for specific curves: see the examples in
[`tnfs/gen/generate_sparse_curve.py`](./tnfs/gen/generate_sparse_curve.py), e.g.
```shell
sage -python -m tnfs.gen.generate_sparse_curve --bls -k 24 -r 256 256 --2NAF -w 6 --find_all_w_up_to
```

# Instantiate a Fotidadis--Martindale curve

E.g. in ePrint [2019/555](eprint.iacr.org/2019/555), there are seeds:
Family 17, seed `ua=-2^64-2^63-2^11-2^10` and `ub=-2^72-2^71-2^36`.
Family 25, `u = -2^64-2^35+2^11-1`, Family 23, `u = 2^48+2^28+2^26`.
```python
from tnfs.curve.fotiadis_martindale import FotiadisMartindale
ua = -2^64-2^63-2^11-2^10
ub = -2^72-2^71-2^36
uc = 2^48+2^28+2^26
ud = -2^64-2^35+2^11-1
Ea = FotiadisMartindale(ua, 17)
Eb = FotiadisMartindale(ub, 17)
Ec = FotiadisMartindale(uc, 23)
Ed = FotiadisMartindale(ud, 25)
for Ei in Ea, Eb, Ec, Ed:
	Ei.print_parameters()
```
