"""
Curves of the paper
  A short-list of pairing-friendly curves resistant to Special TNFS at the 128-bit security level
eprint 2019/1371
Author: Aurore Guillevic

"""
from sage.rings.rational_field import QQ

import tnfs
from tnfs.curve.cyclo_kDe import Cyclo_kDe
from tnfs.curve.cyclo_kDe import polynomial_params
from tnfs.curve.bn import BN
from tnfs.curve.bls12 import BLS12
from tnfs.curve.kss16 import KSS16
from tnfs.curve.kss18 import KSS18
from tnfs.curve.bls24 import BLS24
from tnfs.simul.enumerate_sparse_T import bit_positions_2naf, bits_2naf, bit_positions

proof.arithmetic(False)
QQx = QQ['x']; (x,) = QQx._first_ngens(1)

def print_params(px,rx,betax,lambx):
    px_den = lcm([pi.denom() for pi in px.list()])
    print("px = ({})/{}".format(px*px_den, px_den))
    print("rx = {}".format(rx))
    betax_den = lcm([bi.denom() for bi in betax.list()])
    if betax_den != 1:
        print("betax = ({})/{}".format(betax_den*betax, betax_den))
    else:
        print("betax = {}".format(betax))
    print("lambx = {}".format(lambx))

# Curve FM15
print("\nk=10, D=15, e=1")
k=10
D=15
e0=1
u1=2^32-2^26-2^17+2^10-1
a1=-3
b1=28571049135339031392428419385121263846381182765771408289731088478621591331823569822443993165523925569565418970657909460228657624495497
E1 = Cyclo_kDe(k,D,e0,u1,a=a1,b=b1)
print(E1)
E1.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.cyclo_kDe.polynomial_params(k=k,D=D,e0=e0)
assert (betax^2+15) % px == 0
assert (lambx^2+15) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (2*X^4+X^3-4*X^2+X+2)/(X^3-X) == lambx(X)

print("compute the endomorphism: D=-15 -> CM by (-1+/-sqrt(-15))/2 (x^2+x+4) (or (1+/-sqrt(-15))/2 (x^2-x+4) -> 4-isogeny")
p = E1.p()
r = E1.r()
c = E1.c()
tr = E1.tr()
lambx = (lambx-1)/2
lambx_ = -lambx-1
assert ((lambx**2 + lambx + 4) % rx) == 0
assert ((lambx_**2 + lambx_ + 4) % rx) == 0
lambr = ZZ(lambx(u1))
assert ((lambr**2 + lambr + 4) % r) == 0
lambr_ = -lambr-1
assert ((lambr_**2 + lambr_ + 4) % r) == 0
betax = (betax-1)/2
assert ((betax**2 + betax + 4) % px) == 0
betap = betax(u1)
if not betap in ZZ:
    if betap+p/3 in ZZ:
        betap = betap + p/3
    else:
        betap = betap - p/3
betap = ZZ(betap)
assert ((betap**2 + betap + 4) % p) == 0

for x0,_ in E1.torsion_polynomial(4).roots():
    y02 = x0**3 + E1.a4()*x0 + E1.a6()
    if y02.is_square():
        y0 = y02.sqrt()
        P4 = E1((x0,y0))
        i4 = E1.isogeny(P4)
        j4 = i4.codomain().j_invariant()
        if j4 == E1.j_invariant():
            iso4 = i4
            PP = P4
            print(i4)
            print(i4.rational_maps())
            print(P4)
            print(2*P4)
            continue
P = c*E1.random_element()
while P == E1(0):
    P = c*E1.random_element()
assert P != E1(0) and r*P == E1(0)
E4i = iso4.codomain()
phi4 = E4i.isomorphism_to(E1)
Q = phi4(iso4(P))
R = phi4(iso4(Q))
if R + Q + 4*P != E1(0):
    if (R - Q + 4*P) == E1(0):
        lambr = lambr + 1
        lambr_ = lambr_ + 1
        assert ((lambr**2 - lambr + 4) % r) == 0
        assert ((lambr_**2 - lambr_ + 4) % r) == 0
        print("endomorphism has characteristic polynomial X^2-X+4 on r-torsion points")
else:
    print("endomorphism has characteristic polynomial X^2+X+4 on r-torsion points")

if (lambr*P == Q):
    print("eigenvalue lambr = {}".format(lambr))
elif (lambr_*P == Q):
    print("eigenvalue lambr_ = {}".format(lambr_))
print(iso4.rational_maps())
print(phi4)
"""
QQx<x> := PolynomialRing(QQ);
rx := CyclotomicPolynomial(30);
K<w> := NumberField(rx);
rr := Roots(x^2+x+4, K);
l1 := QQx ! Eltseq(rr[1][1]);
l2 := QQx ! Eltseq(rr[2][1]);
l1; // -x^7 + x^5 + 2*x^4 + x^3 + x^2 - 2*x - 2
l2; // x^7 - x^5 - 2*x^4 - x^3 - x^2 + 2*x + 1
if LeadingCoefficient(l1) lt 0 then
    l_ := l1;
    l1 := l2;
    l2 := l_;
end if;
M := Matrix(QQx, 2, 2, [rx, 0, -l1, 1]);
R := LLL(M);
R;
[          4*x^3 - 4*x x^4 + x^3 - 2*x^2 + 1]
[  x^4 - 2*x^2 + x + 1              -x^3 + x]
a0 := R[1][1]; a1 := R[1][2];
a2 := R[2][1]; a3 := R[2][2];
assert ((a0 + a1*l1) mod rx) eq 0;
(a0 + a1*l1) div rx; // x^3 - 3*x + 1
assert ((a2 + a3*l1) mod rx) eq 0;
(a2 + a3*l1) div rx; // -x^2 + x + 1
M := Matrix(QQx, 2, 2, [rx, 0, -l2, 1]);
R := LLL(M);
R;
[         -4*x^3 + 4*x   x^4 - 2*x^2 + x + 1]
[x^4 + x^3 - 2*x^2 + 1               x^3 - x]
b0 := R[1][1]; b1 := R[1][2];
b2 := R[2][1]; b3 := R[2][2];
assert ((b0 + b1*l2) mod rx) eq 0;
(b0 + b1*l2) div rx; // -x^3 + x^2 + 2*x - 2
assert ((b2 + b3*l2) mod rx) eq 0;
(b2 + b3*l2) div rx; // -x^2 + x + 1
"""
a0 = 4*x**3 - 4*x
a1 = x**4 + x**3 - 2*x**2 + 1
assert a1 == (x+1)*a0/4 -x**2+x+1
a2 = x**4 - 2*x**2 + x + 1
a3 = -x**3 + x
assert a2 == -x*a3 -x**2+x+1

b0 = -4*x**3 + 4*x
b1 = x**4 - 2*x**2 + x + 1
assert b1 == -x*b0/4 -x**2+x+1
b2 = x**4 + x**3 - 2*x**2 + 1
b3 = x**3 - x
assert b2 == (x+1)*b3 -x**2+x+1

assert ((a0 + a1*lambx) % rx) == 0
assert ((a2 + a3*lambx) % rx) == 0
assert ((b0 + b1*lambx_) % rx) == 0
assert ((b2 + b3*lambx_) % rx) == 0

# Curve k11
print("\nk=11, D=3, e=8")
k=11
D=3
e0=8
#u=-2^13+2^10-2^8-2^5-2^3-2
u2=-0x1d2a
a2=0
b2=13
E2 = Cyclo_kDe(k,D,e0,u2,a=a2,b=b2)
print(E2)
E2.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=11,D=3,e0=8)
assert (betax^2+betax +1) % px == 0
assert (lambx^2+lambx +1) % rx == 0
lambx_ = -lambx-1
assert (lambx_^2+lambx_ +1) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (X^10-X^9+X^7-X^6+X^4-X^3+X-1)/(X^9-X^8+X^6-X^5+X^3-X^2+1) == lambx(X)
"""
QQx<x> := PolynomialRing(QQ);
rx := CyclotomicPolynomial(33);
K<w> := NumberField(rx);
rr := Roots(x^2+x+1, K);
l1 := QQx ! Eltseq(rr[1][1]);
l2 := QQx ! Eltseq(rr[2][1]);
l1; // -x^11 - 1
l2; // x^11
if LeadingCoefficient(l1) lt 0 then
    l_ := l1;
    l1 := l2;
    l2 := l_;
end if;
M := Matrix(QQx, 2, 2, [rx, 0, -l1, 1]);
R := LLL(M);
R;
[     x^9 - x^8 + x^6 - x^5 + x^3 - x^2 + 1     x^10 - x^8 + x^7 - x^5 + x^4 - x^2 + x]
[x^10 - x^9 + x^7 - x^6 + x^4 - x^3 + x - 1     -x^9 + x^8 - x^6 + x^5 - x^3 + x^2 - 1]
a0 := R[1][1]; a1 := R[1][2];
a2 := R[2][1]; a3 := R[2][2];
assert ((a0 + a1*l1) mod rx) eq 0;
(a0 + a1*l1) div rx; // x + 1
assert ((a2 + a3*l1) mod rx) eq 0;
(a2 + a3*l1) div rx; // -1
M := Matrix(QQx, 2, 2, [rx, 0, -l2, 1]);
R := LLL(M);
R;
[    -x^9 + x^8 - x^6 + x^5 - x^3 + x^2 - 1 x^10 - x^9 + x^7 - x^6 + x^4 - x^3 + x - 1]
[    x^10 - x^8 + x^7 - x^5 + x^4 - x^2 + x      x^9 - x^8 + x^6 - x^5 + x^3 - x^2 + 1]
b0 := R[1][1]; b1 := R[1][2];
b2 := R[2][1]; b3 := R[2][2];
assert ((b0 + b1*l2) mod rx) eq 0;
(b0 + b1*l2) div rx; // -x
assert ((b2 + b3*l2) mod rx) eq 0;
(b2 + b3*l2) div rx; // -1
"""
a0 = x^9 - x^8 + x^6 - x^5 + x^3 - x^2 + 1
a1 = x^10 - x^8 + x^7 - x^5 + x^4 - x^2 + x
assert a1 == (x+1)*a0 - 1
a2 = x^10 - x^9 + x^7 - x^6 + x^4 - x^3 + x - 1
a3 = -x^9 + x^8 - x^6 + x^5 - x^3 + x^2 - 1
assert a2 == -x*a3 - 1

b0 = -x^9 + x^8 - x^6 + x^5 - x^3 + x^2 - 1
b1 = x^10 - x^9 + x^7 - x^6 + x^4 - x^3 + x - 1
assert b1 == -x*b0 - 1
b2 = x^10 - x^8 + x^7 - x^5 + x^4 - x^2 + x
b3 = x^9 - x^8 + x^6 - x^5 + x^3 - x^2 + 1
assert b2 == (x+1)*b3 - 1

assert ((a0 + a1*lambx) % rx) == 0
assert ((a2 + a3*lambx) % rx) == 0
assert ((b0 + b1*lambx_) % rx) == 0
assert ((b2 + b3*lambx_) % rx) == 0

#second curve k11
print("\nk=11, D=11, e0=4")
k=11
D=11
e0=4
u3=-2^26+2^21+2^19-2^11-2^9-1;
a3=-1056;b3=13552
E3 = Cyclo_kDe(k=11, D=11, e0=4, u=u3,a=a3,b=b3)
E3.j_invariant() == -32768
print(E3)
E3.print_parameters()
Fp = E3._Fp
w = Fp(a3/(2)).sqrt()
a3_ = 2
b3_ = b3/w^3
j3 = 1728*(4*a3^3)/(4*a3^3+27*b3^2)
Fp(j3) == 1728*(4*a3_^3)/(4*a3_^3+27*b3_^2)
E3_ = Cyclo_kDe(k=11, D=11, e0=4, u=u3,a=a3_,b=b3_)
print(E3_)
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=11,D=11,e0=4)
assert (betax^2+11) % px == 0
assert (lambx^2+11) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (2*X^5+X^4-2*X^3+2*X^2-X-2)/(X^4+X) == lambx(X)
# let's compute the endomorphism of characteristic polynomial x^2 +/-x + 3 from a 3-isogeny
# a 3-torsion point is (0, sqrt(E.a6())
p = E3_.p()
r = E3_.r()
c = E3_.c()
tr = E3_.tr()
lambx = (lambx-1)/2
lambx_ = -lambx-1
assert ((lambx**2 + lambx + 3) % rx) == 0
assert ((lambx_**2 + lambx_ + 3) % rx) == 0
lambr = ZZ(lambx(u3))
assert ((lambr**2 + lambr + 3) % r) == 0
lambr_ = -lambr-1
assert ((lambr_**2 + lambr_ + 3) % r) == 0
betax = (betax-1)/2
assert ((betax**2 + betax + 3) % px) == 0
betap = betax(u3)
if betap not in ZZ:
    d = betap.denom()
    dd = (betap.numer() % d)
    dp = (p % d)
    g, ud, vd = xgcd(dp, d)
    #ud * dp + vd * d == 0
    betap = betap + (d - dd)*ud*p/d
betap = ZZ(betap)
assert ((betap**2 + betap + 3) % p) == 0

for x0,_ in E3_.torsion_polynomial(3).roots():
    y02 = x0**3 + E3_.a4()*x0 + E3_.a6()
    if y02.is_square():
        y0 = y02.sqrt()
        P3 = E3_((x0,y0))
        i3 = E3_.isogeny(P3)
        j4 = i3.codomain().j_invariant()
        if j4 == E3_.j_invariant():
            iso3 = i3
            PP = P3
            print(i3)
            print(i3.rational_maps())
            print(P3)
            continue
P = c*E3_.random_element()
while P == E3_(0):
    P = c*E3_.random_element()
assert P != E3_(0) and r*P == E3_(0)
E3i = iso3.codomain()
phi3 = E3i.isomorphism_to(E3_)
Q = phi3(iso3(P))
R = phi3(iso3(Q))
if R + Q + 3*P != E3_(0):
    if (R - Q + 3*P) == E3_(0):
        lambr = lambr + 1
        lambr_ = lambr_ + 1
        assert ((lambr**2 - lambr + 3) % r) == 0
        assert ((lambr_**2 - lambr_ + 3) % r) == 0
        print("endomorphism has characteristic polynomial X^2-X+3 on r-torsion points")
else:
    print("endomorphism has characteristic polynomial X^2+X+3 on r-torsion points")

if (lambr*P == Q):
    print("eigenvalue lambr = {}".format(lambr))
elif (lambr_*P == Q):
    print("eigenvalue lambr_ = {}".format(lambr_))
print(iso3.rational_maps())
print(phi3)
"""
QQx<x> := PolynomialRing(QQ);
rx := CyclotomicPolynomial(11);
K<w> := NumberField(rx);
rr := Roots(x^2+x+3, K);
l1 := QQx ! Eltseq(rr[1][1]);
l2 := QQx ! Eltseq(rr[2][1]);
l1; // -x^9 - x^5 - x^4 - x^3 - x - 1
l2; // x^9 + x^5 + x^4 + x^3 + x
if LeadingCoefficient(l1) lt 0 then
    l_ := l1;
    l1 := l2;
    l2 := l_;
end if;
M := Matrix(QQx, 2, 2, [rx, 0, -l1, 1]);
R := LLL(M);
R;
[              3*x^4 + 3*x x^5 + x^4 - x^3 + x^2 - 1]
[  x^5 - x^3 + x^2 - x - 1                  -x^4 - x]
a0 := R[1][1]; a1 := R[1][2];
a2 := R[2][1]; a3 := R[2][2];
assert ((a0 + a1*l1) mod rx) eq 0;
(a0 + a1*l1) div rx; // x^4 - 2*x^2 + 2*x
assert ((a2 + a3*l1) mod rx) eq 0;
(a2 + a3*l1) div rx; // -x^3 + x^2 - 1
M := Matrix(QQx, 2, 2, [rx, 0, -l2, 1]);
R := LLL(M);
R;
[             -3*x^4 - 3*x   x^5 - x^3 + x^2 - x - 1]
[x^5 + x^4 - x^3 + x^2 - 1                   x^4 + x]
b0 := R[1][1]; b1 := R[1][2];
b2 := R[2][1]; b3 := R[2][2];
assert ((b0 + b1*l2) mod rx) eq 0;
(b0 + b1*l2) div rx; // -x^4 + x^3 + x^2 - 2*x + 1
assert ((b2 + b3*l2) mod rx) eq 0;
(b2 + b3*l2) div rx; // -x^3 + x^2 - 1
"""
a0 = 3*x**4 + 3*x
a1 = x**5 + x**4 - x**3 + x**2 - 1
assert a1 == (x+1)/3*a0 - x**3 - x - 1
a2 = x**5 - x**3 + x**2 - x - 1
a3 = -x**4 - x
assert a2 == -x*a3 -x**3-x-1
b0 = -3*x**4 - 3*x
b1 = x**5 - x**3 + x**2 - x - 1
assert b1 == -x/3*b0 -x**3-x-1
b2 = x**5 + x**4 - x**3 + x**2 - 1
b3 = x**4 + x
assert b2 == (x+1)*b3 -x**3-x-1

assert ((a0 + a1*lambx) % rx) == 0
assert ((a2 + a3*lambx) % rx) == 0
assert ((b0 + b1*lambx_) % rx) == 0
assert ((b2 + b3*lambx_) % rx) == 0

# BN curve of 446 bits from Pereira et al.
print("\nk=12 D=3 BN curve")
k=12
D=3
u4=2^110+2^36+1
a4=0
b4=2^8+1
E4 = BN(u=u4,b=b4)
print(E4)
E4.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bn.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

# BLS12 curve
print("\nk=12 D=3 BLS12 curve")
k=12
D=3
u5=-(2^74+2^73+2^63+2^57+2^50+2^17+1)
a5=0
b5=1
E5 = BLS12(u=u5,b=b5)
print(E5)
E5.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bls12.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

print("\nk=12, D=3, e0=17 Fotiadis-Martindale #17 curve")
k=12
D=3
e0=17
u6=-2^72-2^71-2^36
a6=0
b6=-2
E6 = Cyclo_kDe(k=12, D=3, e0=17, u=u6, b=-2)
print(E6)
E6.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=12,D=3,e0=17)
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

print("\nk=13, D=3, e0=9")
u7 = 0x8b0
a7=0
b7=-17
E7 = Cyclo_kDe(k=13, D=3, e0=9, u=u7, b=-17)
print(E7)
E7.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=13,D=3,e0=9)
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (X^11 - X^10 + X^8 - X^7 + X^5 - X^4 + X^2 - X)/(X^12 - X^11 + X^9 - X^8 + X^6 - X^5 + X^3 - X^2 + 1) == lambx(X)

print("\nk=14, D=3, e0=5")
#u8 = 2^21+2^19+2^10-2^6
u8 = 0x2803c0
a8=0
b8=-4
E8 = Cyclo_kDe(k=14, D=3, e0=5, u=u8, b=-4)
print(E8)
E8.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=14,D=3,e0=5)
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (X^5 + X^4 - X^2 - X)/(X^6 - X^4 - X^3 + X + 1) == lambx(X)

print("\nKSS16 curve k=16 D=1")
u9=-2^34+2^27-2^23+2^20-2^11+1
a9=1
b9=0
E9 = KSS16(u=u9,a=a9)
print(E9)
E9.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.kss16.polynomial_params()
assert (betax^2+1) % px == 0
assert (lambx^2+1) % rx == 0
print_params(px,rx,betax,lambx)

print("\nKSS18 curve k=18 D=3")
u10=2^44+2^22-2^9+2
a10=0
b10=3
E10= KSS18(u=u10,b=b10)
print(E10)
E10.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.kss18.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

# BLS24 curve
print("\nk=24 D=3 BLS24 curve")
k=24
D=3
u11=-2^32+2^28-2^23+2^21+2^18+2^12-1
a11=0
b11=1
E11 = BLS24(u=u11,b=b11)
print(E11)
E11.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bls24.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

print("\n192-bit security level\n")
# more to come
print("\nCyclo (6.4) k=20 D=1 trace = u+1")
E20a = Cyclo_kDe(k=20, D=1, e0=1, u=0xf1a1bf38809c5d, a=-1)
print(E20a)
E20a.print_parameters()
E20b = Cyclo_kDe(k=20, D=1, e0=1, u=0xffffffffffd1ed, a=-1)
print(E20b)
E20b.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=20,D=1,e0=1)
assert (betax^2+1) % px == 0
assert (lambx^2+1) % rx == 0
print_params(px,rx,betax,lambx)

print("\nCyclo (BLS) k=21 D=1 trace = u+1")
E21a = Cyclo_kDe(k=21, D=3, e0=1, u=0xf1a1ddd7, b=1)
print(E21a)
E21a.print_parameters()
E21b = Cyclo_kDe(k=21, D=3, e0=1, u=0xffffccc1, b=1)
print(E21b)
E21b.print_parameters()
E21a_ = Cyclo_kDe(k=21, D=3, e0=1, u=-0xf1a1c083, b=1)
print(E21a_)
E21a_.print_parameters()
E21b_ = Cyclo_kDe(k=21, D=3, e0=1, u=-0xffff6fd1, b=1)
print(E21b_)
E21b_.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = polynomial_params(k=21,D=3,e0=1)
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)
