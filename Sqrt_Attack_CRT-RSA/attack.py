from Crypto.Util.number import getPrime
import random as rd
from sage.all import *

# start parameter generation

def gen_key():
	while True:
		p = getPrime(512)
		q = getPrime(512)
		if GCD(p-1, q-1) == 2:
			return p, q

def get_clue(p, q, BITS):
	while True:
		d_p = rd.randint(1, 1 << BITS)
		d_q = rd.randint(1, q - 1)
		if GCD(d_p, p - 1) == 1 and GCD(d_q, q - 1) == 1:
			d = crt(d_p, d_q, p - 1, q - 1)
			e = int(inverse_mod(d, (p - 1) * (q - 1)))
			return e

def gen_param():
	p, q = gen_key()
	e = get_clue(p, q, 36)
	n = p * q
	return n, e 

print("[+] Generating Parameters")
sys.setrecursionlimit(10 ** 6)
N, e = gen_param()

# start solving

RR = Zmod(N)
P = PolynomialRing(RR, 'x')
x = P.gen()
m = 53 # choose randomly
M = 1 << 18
bb = RR(m) ** (e * M)

POLYS = [0] * (4 * M)
REM = [0] * (4 * M)
ans = [0] * M

# divide & conquer
def calc_poly(L, R):
	if L >= R:
		return 1
	if L + 1 == R:
		return (bb ** L) * x - m
	f1 = calc_poly(L, (L+R)//2)
	f2 = calc_poly((L+R)//2, R)
	return f1 * f2

# segment tree build
def build(index, L, R):
	if L >= R:
		return
	if L + 1 == R:
		POLYS[index] = x - (RR(m) ** (e * L))
		return
	build(index << 1, L, (L + R) // 2)
	build(index << 1 | 1, (L + R) // 2, R)
	POLYS[index] = POLYS[index << 1] * POLYS[index << 1 | 1]

# multipoint evaluation
def calc(index, L, R):
	if L >= R:
		return
	if L + 1 == R:
		ans[L] = REM[index].coefficients()[0]
		return
	REM[index << 1] = REM[index] % POLYS[index << 1]
	REM[index << 1 | 1] = REM[index] % POLYS[index << 1 | 1]
	calc(index << 1, L, (L + R) // 2)
	calc(index << 1 | 1, (L + R) // 2, R)

print("[+] Calculate Poly")
f = calc_poly(0, M)
print("[+] Build Segtree")
build(1, 0, M)
REM[1] = f % POLYS[1]
print("[+] Get Remainders")
calc(1, 0, M)

print("[+] Begin Search")
for i in range(0, M):
	TT = GCD(int(ans[i]), int(N))
	if TT != 1 and TT != N:
		p = int(TT)
		q = N // p
		print(p, q)
		assert p * q == N
		break