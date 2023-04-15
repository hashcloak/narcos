"""
Curves of the paper
  Updating Key Size Estimations for Pairings
eprint 2017/334
Authors: Razvan Barbulescu and Sylvain Duquesne

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

print("128-bit security level")
# Curve BN
print("\nk=12, D=3, BN curve")
k=12
D=3
u1=2^114+2^101-2^14-1
a1=0
b1=-4
E1 = BN(u1,b=b1)
print(E1)
E1.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bn.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (2*X+1)/(6*X^2+2*X) == lambx(X)
assert -(6*X^2+4*X+1)/(2*X+1) == lambx(X)

# Curve BLS12-1
print("\nk=12, D=3, BLS12 curve")
k=12
D=3
e0=1
u2=-2^77+2^50+2^33
a2=0
b2=4
E2 = BLS12(u2,b=b2)
print(E2)
E2.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bls12.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
lambx_ = -lambx-1
assert (lambx_^2+lambx_+1) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)

#second curve BLS12
u3=-2^77-2^59+2^9
b3=4
E3 = BLS12(u3,b=b3)
print(E3)
E3.print_parameters()

#third curve BLS12
u4=-2^77-2^71-2^64+2^37+2^35+2^22-2^5
b4=-2
E4 = BLS12(u4,b=b4)
print(E4)
E4.print_parameters()

"""
QQx<x> := PolynomialRing(QQ);
rx := 36*Factorization(Evaluate(CyclotomicPolynomial(12), 6*x^2))[2][1];
//rx := CyclotomicPolynomial(12);
K<w> := NumberField(rx);
rr := Roots(x^2+x+1, K);
l1 := QQx ! Eltseq(rr[1][1]);
l2 := QQx ! Eltseq(rr[2][1]);
l1; // -x^2
// 36*x^3 - 18*x^2 - 6*x - 2
l2; // x^2-1
// 36*x^3 + 18*x^2 + 6*x + 1
if LeadingCoefficient(l1) lt 0 then
    l_ := l1;
    l1 := l2;
    l2 := l_;
end if;
M := Matrix(QQx, 2, 2, [rx, 0, -l1, 1]);
R := LLL(M);
//R;
R := 6*R; R;
//[      1     x^2]
//[x^2 - 1      -1]
//[       -2*x - 1     6*x^2 + 2*x]
//[6*x^2 + 4*x + 1         2*x + 1]
a0 := R[1][1]; a1 := R[1][2];
a2 := R[2][1]; a3 := R[2][2];
assert ((a0 + a1*l1) mod rx) eq 0;
(a0 + a1*l1) div rx; // 1
// 6*x - 1
assert ((a2 + a3*l1) mod rx) eq 0;
(a2 + a3*l1) div rx; // -1
// 2
M := Matrix(QQx, 2, 2, [rx, 0, -l2, 1]);
R := LLL(M);
// R;
R := 6*R; R;
//[     -1 x^2 - 1]
//[    x^2       1]
//[        2*x + 1 6*x^2 + 4*x + 1]
//[    6*x^2 + 2*x        -2*x - 1]
b0 := R[1][1]; b1 := R[1][2];
b2 := R[2][1]; b3 := R[2][2];
assert ((b0 + b1*l2) mod rx) eq 0;
(b0 + b1*l2) div rx; // -1
// -6*x - 1
assert ((b2 + b3*l2) mod rx) eq 0;
(b2 + b3*l2) div rx; // 0
// 2
"""
# KSS16 curves
print("\nk=16, D=4, KSS16")
k=16
D=4
u5=-2**34+2**27-2**23+2**20-2**11+1
a5=1;b5=0
E5 = KSS16(u=u5,a=a5)
print(E5)
E5.print_parameters()

u6=2**35-2**32-2**18+2**8+1
a6=1;b6=0
E6 = KSS16(u=u6,a=a6)
print(E6)
E6.print_parameters()

px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.kss16.polynomial_params()
assert (betax^2+1) % px == 0
assert (lambx^2+1) % rx == 0
print_params(px,rx,betax,lambx)
Kr.<X> = NumberField(rx)
assert (X^4+24)/7 == lambx(X)

print("\nKSS18 curve k=18 D=3")
u7=2^44+2^22-2^9+2
a7=0
b7=3
E7= KSS18(u=u7,b=b7)
print(E7)
E7.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.kss18.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

print("192-bit security level")

print("\nKSS18 curve k=18 D=3")
u8=-2**85-2**31-2**26+2**6
a8=0
b8=2
E8= KSS18(u=u8,b=b8)
print(E8)
E8.print_parameters()

# BLS24 curve
print("\nk=24 D=3 BLS24 curve")
k=24
D=3
u9=-2**56-2**43+2**9-2**6
a9=0
b9=-2
E9 = BLS24(u=u9,b=b9)
print(E9)
E9.print_parameters()
px,rx,tx,cx,yx,betax,lambx,D = tnfs.curve.bls24.polynomial_params()
assert (betax^2+betax+1) % px == 0
assert (lambx^2+lambx+1) % rx == 0
print_params(px,rx,betax,lambx)

print("\n256-bit security level\n")
# KSS18
u10=2**186-2**75-2**22+2**4
b10=-2
E10= KSS18(u=u10,b=b10)
print(E10)
E10.print_parameters()

# BLS24
u11=-2**103-2**101+2**68+2**50
b11=-2
E11 = BLS24(u=u11,b=b11)
print(E11)
E11.print_parameters()

