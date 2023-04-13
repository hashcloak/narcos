from narcos_generation_script import test_finite_field_nfs, test_finite_field_tnfs
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.all import EllipticCurve, power_mod, GF


q = ZZ(90924519629997939991666979336023169281404607245509258933367291716054347256010286951411976973635670301538812567627538336231454957950025633906465864404045313853)
r = ZZ(115792089210356248762697446949407573530086143415290314195533631308867097853951)
k = 6
tr = ZZ(9459927709287013385885240241424446084893356409782777942811485615308356970762977)
D = 3
j = hilbert_class_polynomial(-D).roots()[0][0]

# TODO: Check Twist Security
# Since D = 3 and 6 | K, we can use a sextic twist where its degree is d = 6
# E = EllipticCurve(GF(q), j=j)
# Find i such that i is a quadratic and cubic non-residue
""" i = None
for a in range(1, q):
    if power_mod(a, (q-1)//2, q) == 1 or power_mod(a, (q-1)//3, q) == 1:
        continue # a is a quadratic or cubic residue, skip to next a
    if power_mod(a, (q-1)//6, q) != 1:
        i = a
        break
twist = E.sextic_twist(i) """

def get_estimated_bits_of_security(q, r, k):
    bitsecs = [0] * 9

    bitsecs[0] = test_finite_field_nfs(q, r, k, cost=118, sieving_dim=2, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    bitsecs[1] = test_finite_field_nfs(q, r, k, cost=117, sieving_dim=3, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    bitsecs[2] = test_finite_field_nfs(q, r, k, cost=120, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    bitsecs[3] = test_finite_field_nfs(q, r, k, cost=135, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    bitsecs[4] = test_finite_field_tnfs(q, r, k, cost=156, samples=100000, conj=True, max_coeff=2, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    bitsecs[5] = test_finite_field_tnfs(q, r, k, cost=130, samples=100000, jouxlercier=True, max_coeff=1, deg_f=3, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    bitsecs[6] = test_finite_field_tnfs(q, r, k, cost=118, samples=100000, jouxlercier=True, max_coeff=1, deg_f=4, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    bitsecs[7] = test_finite_field_tnfs(q, r, k, cost=115, samples=100000, jouxlercier=True, max_coeff=1, deg_f=5, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    bitsecs[8] = test_finite_field_tnfs(q, r, k, cost=115, samples=100000, jouxlercier=True, max_coeff=1, deg_f=6, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    return min(bitsecs)

print("------------------------------------------")
print("Estimated bits of security for E: ", get_estimated_bits_of_security(q,r,k))