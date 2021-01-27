# extremely rough showcase of Diophantine Argument of Knowledge
# we assume honest verifiers. also, the algorithm is interactive
# we could (?) add fiat-shamir heuristic to make it non-interactive
# this code does not really care about reducing proof sizes as long as they are doable

import random
from Crypto.Util.number import inverse, getPrime, isPrime
# support functions


def gpow(a, b, c): # a^b mod c
	if b >= 0:
		return pow(a, b, c)
	else:
		return pow(inverse(a, c), -b, c)
	
def genSafePrime(nbits, evade = 0): # generate safe prime not equal to evade
	while True:
		p = getPrime(nbits)
		if p != evade and isPrime((p-1)//2):
			return p
	
def get_gen(p): # get element of order (p-1) // 2
	while True:
		u = random.randrange(0, p)
		if u != 1 and pow(u, (p-1)//2, p) == 1:
			return u

def CRT(a, m, b, n): # we assume gcd(m, n) = 1
	mul = ((b - a) * inverse(m%n, n)) % n
	ret = (mul * m + a) % (m * n)
	ret = (ret + m * n) % (m * n)
	assert ret % m == a and ret % n == b
	return ret

def sqrt(n): # integer sqrt with binary search
	if n == 0:
		return 0
	lef, rig = 1, 2
	while rig ** 2 < n:
		rig = rig << 1
	while lef <= rig:
		mid = (lef + rig) // 2
		if mid ** 2 <= n:
			best, lef = mid, mid + 1
		else:
			rig = mid - 1
	return best

def tonelli(tar, x): # tonelli-shanks algorithm
	tar %= x
	if pow(tar, (x-1) // 2, x) != 1:
		return -1
	S, tp, cc, sz = 0, x, x-1, 0
	if tar == 0:
		return 0
	if x % 4 == 3:
		return pow(tar, (x+1) // 4, x)
	while cc % 2 == 0:
		cc, S = cc // 2, S + 1
	Q, z = cc, 2
	while pow(z, (x-1)//2, x) == 1:
		z += 1
	M, c, t, R = S, pow(z, Q, x), pow(tar, Q, x), pow(tar, (Q+1) // 2, x)
	while True:
		if t == 0:
			return 0
		if t == 1:
			return R % x
		sx, tem = 0, t
		while tem != 1:
			sx, tem = sx + 1, (tem * tem) % x
		b = pow(c, 1<<(M-sx-1), x)
		M, c = sx, (b * b) % x
		t, R = (t * c) % x, (R * b) % x

def abs(x):
	if x >= 0:
		return x
	else:
		return -x 

def prime_two_square(p):
	u = tonelli(p-1, p)
	if u <= p // 2:
		u = p - u
	a, b = p, u
	l = sqrt(p)
	while b >= l:
		r = a % b 
		a, b = b, r 
	c = sqrt(p - b * b)
	return b, c

def brute_four_square(x):
	sq = sqrt(x) + 1
	for i in range(0, sq):
		for j in range(i, sq):
			if i * i + j * j > x:
				break
			for k in range(j, sq):
				if i * i + j * j + k * k > x:
					break
				for l in range(k, sq):
					if i * i + j * j + k * k + l * l == x:
						return [i, j, k, l]

def four_square(x):
	if x < 0:
		ret = [random.randrange(0, 1 << 512) for i in range(0, 4)]
		return ret
	if x == 0:
		return [0, 0, 0, 0]
	ex = 0
	while x % 2 == 0:
		ex, x = ex + 1, x // 2
	if ex != 1 and ex % 2 == 1:
		ret = four_square(2 * x)
		for i in range(0, 4):
			ret[i] = ret[i] << ((ex - 1) // 2)
		return ret
	if ex % 2 == 0:
		ret = four_square(2 * x)
		for i in range(2, 3):
			if ret[i] % 2 == ret[0] % 2:
				ret[i], ret[1] = ret[1], ret[i]
		fin = []
		s = (1 << (ex // 2))
		fin.append((ret[0] + ret[1]) * s // 2)
		fin.append(abs(ret[0] - ret[1]) * s // 2)
		fin.append((ret[2] + ret[3]) * s // 2)
		fin.append(abs(ret[2] - ret[3]) * s // 2)
		return fin
	if ex == 1:
		num = 2 * x
		if num <= 500:
			return brute_four_squre(num)
		else:
			while True:
				w1 = random.randint(0, sqrt(num))
				w2 = random.randint(0, sqrt(num - w1 * w1))
				p = num - w1 * w1 - w2 * w2
				if p % 4 == 1 and isPrime(p):
					w3, w4 = prime_two_square(p)
					return [w1, w2, w3, w4]


class Verifier:
	def set_prover(self, A):
		self.prover = A
	
	def calc(self, x, r):
		return (gpow(self.g, x, self.N) * gpow(self.h, r, self.N)) % self.N

	def gen_param(self, nbits = 512, test = False): # generate parameters
		self.T = (1 << 1024)
		if test == False:
			self.nbits = nbits
			self.P = genSafePrime(nbits)
			self.Q = genSafePrime(nbits, self.P)
		else:
			self.nbits, nbits = 512, 512
			self.P = 10212507832725305011944042025376841364170191919743201627456683549189969962823106284511698786491897605084563006187459969219868985667022529399803242095049903
			self.Q = 12402356134669021912727815960862954169906553067507512563261563050098281345734222140795387822679489071990409871271079215357414720073401705325543079611345367
		self.k = 2 * nbits
		self.B = 2 * nbits
		self.C = (1 << (nbits // 2))
		self.N = self.P * self.Q
		
		self.p = (self.P - 1) // 2
		self.q = (self.Q - 1) // 2
		
		gp = get_gen(self.P)
		gq = get_gen(self.Q)
		self.h = CRT(gp, self.P, gq, self.Q)
		
		self.alpha = random.randrange(0, 1 << (self.B + self.k))
		self.g = pow(self.h, self.alpha, self.N)
		
		
	def send_param(self):
		return self.T, self.nbits, self.N, self.k, self.B, self.C, self.g, self.h

	def recv_commit(self, commit):
		self.commit = commit
		
	def select_challenge(self):
		self.challenge = random.randrange(0, self.C)
		return self.challenge
	
	def check_open(self, commit, x, r):
		return commit == self.calc(x, r)

	def verify_dlogh(self, z):
		if len(self.commit) != 2:
			return False
		LHS = gpow(self.h, z, self.N)
		RHS = (self.commit[1] * pow(self.commit[0], self.challenge, self.N)) % self.N
		if LHS == RHS:
			return True
		else:
			return False

	def verify_open(self, u, v):
		if len(self.commit) != 2:
			return False
		LHS = (pow(self.g, u, self.N) * pow(self.h, v, self.N)) % self.N
		RHS = (self.commit[1] * pow(self.commit[0], self.challenge, self.N)) % self.N
		if LHS == RHS:
			return True
		else:
			return False

	def verify_mul(self, u, v_1, v_2):
		if len(self.commit) != 5:
			return False
		LHS_1 = self.calc(u, v_1)
		RHS_1 = (self.commit[3] * pow(self.commit[1], self.challenge, self.N)) % self.N
		LHS_2 = (pow(self.commit[0], u, self.N) * gpow(self.h, v_2, self.N)) % self.N
		RHS_2 = (self.commit[4] * pow(self.commit[2], self.challenge, self.N)) % self.N
		if LHS_1 == RHS_1 and LHS_2 == RHS_2:
			return True 
		else:
			return False
	
class Prover:
	def set_verifier(self, V):
		self.verifier = V
		
	def recv_param(self):
		self.T, self.nbits, self.N, self.k, self.B, self.C, self.g, self.h = self.verifier.send_param()
	
	def send_commit(self, commit):
		self.verifier.recv_commit(commit)
	
	def recv_challenge(self):
		return self.verifier.select_challenge()
   
	def calc(self, x, r):
		return (gpow(self.g, x, self.N) * gpow(self.h, r, self.N)) % self.N

	def commit(self, x):
		r = random.randrange(0, 1 << (self.B + self.k))
		commit = self.calc(x, r)
		return r, commit

	def open(self, commit, x, r, display = True):
		if display == True:
			print("[+] Alice opening commitment...")
		result = self.verifier.check_open(commit, x, r)
		return result

	def display(self, result):
		if result == True:
			print("[+] Proof Success!")
		else:
			print("[-] Proof Failed..")

	def prove_dlogh(self, target, r, display = True):
		if display:
			print("[+] Alice wants to prove that she knows the dlog!")
		v = random.randrange(0, (1 << (self.B + 2 * self.k)))
		a = pow(self.h, v, self.N)
		self.send_commit([target, a])
		c = self.recv_challenge()
		result = self.verifier.verify_dlogh(v + c * r)
		if display:
			self.display(result)
		return result

	def prove_open(self, commit, x, r, display = True):
		if display:
			print("[+] Alice wants to prove that she knows how to open!")
		y = random.randrange(0, self.T * self.C * (1 << self.k))
		s = random.randrange(0, self.C * (1 << (self.B + 2 * self.k)))
		d = self.calc(y, s)
		if display:
			print("[+] Sending commit and values, Receiving challenge...")
		self.send_commit([commit, d])
		c = self.recv_challenge()
		u = y + c * x
		v = s + c * r
		if display:
			print("[+] Sending answer, waiting for verification...")
		result = self.verifier.verify_open(u, v)
		if display:
			self.display(result)
		return result

	def trial_open(self):
		commit = self.calc(525, 3751)
		result = self.prove_open(commit, 525, 3751)
		assert result == True
		commit = self.calc(528, 3751)
		result = self.prove_open(commit, 525, 3751)
		assert result == False

	def prove_sum(self, commits, xs, rs, display = True):
		if len(commits) != 3 or len(xs) != 3 or len(rs) != 3:
			print("[-] prove_sum : incorrect array size")
			return
		if display:
			print("[+] Alice wants to prove that sum of two values equal to another!")
		for i in range(0, 3):
			result = self.prove_open(commits[i], xs[i], rs[i], False)
			if result == False:
				return False
		target = (commits[2] * inverse(commits[0] * commits[1], self.N)) % self.N
		true_r = rs[2] - rs[0] - rs[1]
		if display:
			print("[+] Proving dlog knowledge...")
		result = self.prove_dlogh(target, true_r)
		if display:
			self.display(result)
		return result

	def trial_sum(self):
		xs = [1378, 67286, 1378+67286]
		rs = [100, 300, 1678]
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_sum(commits, xs, rs)
		assert result == True
		xs[2] = 67265
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_sum(commits, xs, rs)
		assert result == False
	
	def recv_sum(self, commits, xs, rs):
		x = xs[0] + xs[1]
		r, commit = self.commit(x)
		result = self.prove_sum([commits[0], commits[1], commit], [xs[0], xs[1], x], [rs[0], rs[1], r])
		assert result == True
		return x, r, commit

	def recv_sub(self, commits, xs, rs):
		x = xs[0] - xs[1]
		r, commit = self.commit(x)
		result = self.prove_sum([commits[1], commit, commits[0]], [xs[1], x, xs[0]], [rs[1], r, rs[0]])
		assert result == True
		return x, r, commit

	def prove_mul(self, commits, xs, rs, display = True):
		if len(commits) != 3 or len(xs) != 3 or len(rs) != 3:
			print("[-] prove_sum : incorrect array size")
			return
		if display:
			print("[+] Alice wants to prove that product of two values equal to another!")
		for i in range(0, 3):
			result = self.prove_open(commits[i], xs[i], rs[i], False)
			if result == False:
				return False
		y = random.randrange(0, self.T * self.C * (1 << self.k))
		s_1 = random.randrange(0, self.C * (1 << (self.B + 2 * self.k)))
		s_2 = random.randrange(0, self.C * self.T * (1 << (self.B + 2 * self.k)))
		d_1 = self.calc(y, s_1)
		d_2 = (pow(commits[0], y, self.N) * pow(self.h, s_2, self.N)) % self.N
		if display:
			print("[+] Sending commit and values, Receiving challenge...")
		self.send_commit([commits[0], commits[1], commits[2], d_1, d_2])
		c = self.recv_challenge()
		u = y + c * xs[1]
		v_1 = s_1 + c * rs[1]
		v_2 = s_2 + c * (rs[2] - xs[1] * rs[0])
		if display:
			print("[+] Sending answer, waiting for verification...")
		result = self.verifier.verify_mul(u, v_1, v_2)
		if display:
			self.display(result)
		return result

	def trial_mul(self):
		xs = [1378, 67286, 1378*67286]
		rs = [100, 300, 1678]
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_mul(commits, xs, rs)
		assert result == True
		xs[2] = 67265
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_mul(commits, xs, rs)
		assert result == False

	def recv_mul(self, commits, xs, rs):
		x = xs[0] * xs[1]
		r, commit = self.commit(x)
		result = self.prove_mul([commits[0], commits[1], commit], [xs[0], xs[1], x], [rs[0], rs[1], r])
		assert result == True
		return x, r, commit

	def recv_sq(self, commit, x, r):
		return self.recv_mul([commit, commit], [x, x], [r, r])

	def prove_sumsq(self, S, display = True):
		if display:
			print("[+] Preparing commits...")
		num = len(S)
		xs, rs, commits = [], [], []
		for i in range(0, num):
			r, commit = self.commit(S[i])
			xs.append(S[i]), rs.append(r), commits.append(commit)
		for i in range(0, num):
			x, r, commit = self.recv_sq(commits[i], xs[i], rs[i])
			xs.append(x), rs.append(r), commits.append(commit)
		x, r, commit = self.recv_sum([commits[num], commits[num+1]], [xs[num], xs[num+1]], [rs[num], rs[num+1]])
		xs.append(x), rs.append(r), commits.append(commit)
		for i in range(2, num):
			x, r, commit = self.recv_sum([commits[num+i], commits[2*num-2+i]], [xs[num+i], xs[2*num-2+i]], [rs[num+i], rs[2*num-2+i]])
			xs.append(x), rs.append(r), commits.append(commit)
		return xs[3*num-2], rs[3*num-2], commits[3*num-2]

	def prove_nonneg(self, commit, x, r, display = True): # I want to prove that x >= 0
		if display:
			print("[+] Alice wants to prove that her value is nonnegative!")
		self.prove_open(commit, x, r, False)
		if display:
			print("[+] Writing x as sum of four squares...")
		S = four_square(x)
		x_final, r_final, commit_final = self.prove_sumsq(S)
		target = (commit_final * inverse(commit, self.N)) % self.N 
		new_r = (r_final - r)
		result = True
		if display:
			print("[+] Finally, prove the summation result is equal to original x.")
		if self.prove_dlogh(target, new_r, False) == False:
			result = False
		if display:
			self.display(result)
		return result

	def trial_nonneg(self):
		x = 178716
		r, commit = self.commit(x)
		result = self.prove_nonneg(commit, x, r)
		assert result == True
		x = -47817
		r, commit = self.commit(x)
		result = self.prove_nonneg(commit, x, r)
		assert result == False

	def prove_minus(commits, xs, rs, display = True):
		if len(commits) != 3 or len(xs) != 3 or len(rs) != 3:
			print("[-] prove_sum : incorrect array size")
			return
		if display:
			print("[+] Alice wants to prove a - b = c")
		for i in range(0, 3):
			result = self.prove_open(commits[i], xs[i], rs[i], False)
			if result == False:
				return False
		for i in range(0, 3):
			result = self.prove_nonneg(commits[i], xs[i], rs[i], False)
			if result == False:
				return False
		# a - b = c <=> a - b >= 0 and a = b + c or a <= b and c = 0
		# VAL = ((a - b - x^2 - y^2 - z^2 - w^2)^2 + (a - b - c)^2) * ((b - a - xx^2 - yy^2 - zz^2 - ww^2)^2 + c^2)
		SEQ_1 = four_square(xs[0] - xs[1])
		SEQ_2 = four_square(xs[1] - xs[0])
		SEQ_1_x, SEQ_1_r, SEQ_1_commit = self.prove_sumsq(SEQ_1)
		SEQ_2_x, SEQ_2_r, SEQ_2_commit = self.prove_sumsq(SEQ_1)
		# preparing first equation
		0_1_x, 0_1_r, 0_1_commit = self.recv_sub([commits[0], commits[1]], [xs[0], xs[1]], [rs[0], rs[1]])
		




	def trial_minus(self):
		xs = [1378, 67286, 0]
		rs = [100, 300, 1678]
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_minus(commits, xs, rs)
		assert result == True
		xs[1] , xs[2] = 500, 878
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_minus(commits, xs, rs)
		assert result == True
		xs[1] , xs[2] = 1377, 0
		commits = [self.calc(xs[i], rs[i]) for i in range(0, 3)]
		result = self.prove_minus(commits, xs, rs)
		assert result == False

# begin main protocol setup
print("[+] Starting Protocol...")
Alice = Prover()
Bob = Verifier()
Alice.set_verifier(Bob)
Bob.set_prover(Alice)

print("[+] Generating & Sending Parameters")
Bob.gen_param(512, True)
Alice.recv_param()

print("[*] Parameter Sent. Originally, Bob now has to prove knowledge of alpha.")
print("[*] This is skipped in this program for the sake of simplicity.")

print("[*] Proving Alice knows how to open.")
Alice.trial_open()

print("[*] Proving Alice has x_1 + x_2 = x_3.")
Alice.trial_sum()

print("[*] Proving Alice has x_1 * x_2 = x_3.")
Alice.trial_mul()

print("[*] Proving Alice has x >= 0")
Alice.trial_nonneg()

print("[*] Bounded Arithmetic ShowCase : a - b = c")

# Alice.trial_minus()

print("[*] Bounded Arithmetic ShowCase : MSB(a, k) = b")
# Alice.trial_MSB()
