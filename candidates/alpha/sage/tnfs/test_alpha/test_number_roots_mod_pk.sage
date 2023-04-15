from tnfs.alpha.number_roots_mod_pk import no_roots_mod_pk

# tests

ZZx.<x> = ZZ[]
f = x^5 + 12*x^3 + 12*x^2 - 11*x + 8

for p in [2,3,5]:
    print("p={}".format(p))
    for k in range(1,10):
        n_pk = no_roots_mod_pk(f, p, k)
        print("{:4d} roots mod {}^{}".format(n_pk,p,k))

#f = ZZx ([randint(-12,12) for i in range(5)] + [2*3*5])
#while not (f.content() == 1 and f.is_irreducible() and (f.discriminant() % (16*81*125)) == 0):
#    f = ZZx ([randint(-12,12) for i in range(5)] + [2*3*5])

f = 30*x^5 + 12*x^4 - 3*x^3 + 2*x^2 - 3*x + 12
#disc(f) = 59155942680000 = 2^6 * 3^4 * 5^4 * 18258007

disc_f = f.discriminant()
print("f = {}".format(f))
print("disc(f) = {} = {}".format(disc_f, disc_f.factor()))

for p in [2,3,5]:
    print("p={}".format(p))
    for k in range(1,10):
        n_pk = no_roots_mod_pk(f, p, k)
        print("{:4d} roots mod {}^{}".format(n_pk,p,k))

# pb of definition of projective roots

# examples from Kopp, Randall, Rojas, Zhu
# https://arxiv.org/abs/1808.10531
# Table p. 2

f = (x-1234)^3*(x-7193)^4*(x-2030)^12
disc_f = f.discriminant()
print("f = {}".format(f))
print("f = {}".format(f.factor()))
print("disc(f) = {}".format(disc_f))

p = 123456791
k=1
n_pk = no_roots_mod_pk(f, p, k)
print("{:4d} roots mod {}^{}".format(n_pk,p,k))

k=23
n_pk = no_roots_mod_pk(f, p, k)
print("{:4d} roots mod {}^{}".format(n_pk,p,k))

ans = 83524650739763670783591272793501499347381420700990366689774050080031654011699848668752654473531540039924209209663876325122031629580404523246324540823308088725469492593973

ans == n_pk

# this function is much faster

#Example 1.6
f = x^10 - 10*x + 738
disc_f = f.discriminant()
print("f = {}".format(f))
print("disc(f) = {} = {}".format(disc_f, disc_f.factor()))

for p in [2,3,5,7]:
    print("p={}".format(p))
    for k in range(1,8):
        n_pk = no_roots_mod_pk(f, p, k)    
        print("{:4d} roots mod {}^{}".format(n_pk,p,k))

