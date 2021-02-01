from Crypto.Util.number import inverse, getPrime
from Crypto.Util.number import long_to_bytes, inverse, bytes_to_long, isPrime, GCD
import random

cnt = 0

p = getPrime(512)
q = getPrime(512)
n = p * q
phi = (p - 1) * (q - 1)

k = 128 # byte length
B = 1 << (8 * k - 8) # bound

e = 65537
d = inverse(e, phi)

def oracle(c): # Oracle
	global cnt
	cnt += 1
	ptxt = pow(c, d, n)
	if ptxt < B:
		return True
	else:
		return False

def Manger_Attack(c):
	f1 = 2
	while True:
		val = (pow(f1, e, n) * c) % n
		if oracle(val):
			f1 = 2 * f1
		else:
			break
	f12 = f1 // 2
	f2 = ((n + B) // B) * f12
	while True:
		val = (pow(f2, e, n) * c) % n
		if oracle(val):
			break
		else:
			f2 += f12
	m_min = (n + f2 - 1) // f2
	m_max = (n + B) // f2
	# note the ERRATA from https://github.com/GDSSecurity/mangers-oracle
	while m_min < m_max:
		f_tmp = (2 * B) // (m_max - m_min)
		I = (f_tmp * m_min) // n 
		f3 = (I * n + m_min - 1) // m_min
		val = (pow(f3, e, n) * c) % n
		if oracle(val):
			m_max = (I * n + B) // f3
		else:
			m_min = (I * n + B + f3 - 1) // f3
	return m_min

# test case
for _ in range(100):
	cnt = 0
	m = random.randrange(1, B)
	c = pow(m, e, n)
	recover = Manger_Attack(c)
	assert m == recover
	print("[+] Recovered in {} Oracle Calls!".format(cnt))
