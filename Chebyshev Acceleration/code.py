import sys
import numpy as np
from scipy.linalg import eigh

n = 2000
A = np.random.randn(n, n)
A = A @ A.T + 30 * np.eye(n)
b = np.random.randn(n, 1)

opt = np.linalg.solve(A, b)

mu = eigh(A, eigvals_only = True, eigvals = (0, 0))[0]
L = eigh(A, eigvals_only = True, eigvals = (n-1, n-1))[0]

print(mu, L)

# begin Gradient Descent

print("[+] Begin Gradient Descent")
gamma = 2 / (mu + L)
x = np.zeros((n, 1))
for i in range(5000):
	x = x - gamma * (A @ x - b)
	if i % 100 == 0:
		print(i, np.linalg.norm(x - opt))

# begin Chebyshev Method

print("[+] Begin Chebyshev Method")
delta = 0
x_0 = np.zeros((n, 1))
x_1 = x_0 - 2 / (L + mu) * (A @ x_0 - b)
for i in range(500):
	delta = (L - mu) / (2 * (L + mu) - delta * (L - mu))
	x_nxt = x_1 - 4 * delta / (L - mu) * (A @ x_1 - b) + (1 - 2 * delta * (L + mu) / (L - mu)) * (x_0 - x_1)
	x_0, x_1 = x_1, x_nxt
	if i % 10 == 0:
		print(i, np.linalg.norm(x_1 - opt))