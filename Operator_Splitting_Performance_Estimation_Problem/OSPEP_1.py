import cvxpy as cp 
import matplotlib.pyplot as plt 
import numpy as np
import math

def THM_4_1_scaled(theta, mu, beta):
	assert 0 < theta < 2 and mu > 0 and beta > 0

	if mu * beta - mu + beta < 0 and theta <= 2 * (beta + 1) * (mu - beta - mu * beta) / (mu + mu * beta - beta - beta * beta - 2 * mu * beta * beta):
		return 1, abs(1 - theta * beta / (beta + 1))
	if mu * beta - mu - beta > 0 and theta <= 2 * (mu * mu + beta * beta + mu * beta + mu + beta - mu * mu * beta * beta) / (mu * mu + beta * beta + mu * mu * beta + mu * beta * beta + mu + beta - 2 * mu * mu * beta * beta):
		return 2, abs(1 - theta * (1 + mu * beta) / (1 + mu) / (1 + beta))
	if theta >= 2 * (mu * beta + mu + beta) / (2 * mu * beta + mu + beta):
		return 3, abs(1 - theta)
	if mu * beta + mu - beta < 0 and theta <= 2 * (mu + 1) * (beta - mu - mu * beta) / (beta + mu * beta - mu - mu * mu - 2 * mu * mu * beta):
		return 4, abs(1 - theta * mu / (mu + 1))
	
	ret = (2 - theta) / 4
	ret = ret * ((2 - theta) * mu * (beta + 1) + theta * beta * (1 - mu))
	ret = ret * ((2 - theta) * beta * (mu + 1) + theta * mu * (1 - beta))
	ret = ret / mu / beta
	ret = ret / (2 * mu * beta * (1 - theta) + (2 - theta) * (mu + beta + 1))
	return 5, math.sqrt(ret)

def THM_4_1_unscaled(alpha, theta, mu, beta):
	return THM_4_1_scaled(theta, mu * alpha, beta / alpha)

def gen_M_O(theta):
	ret = [[1, theta, -theta], 
	       [theta, theta * theta, -theta * theta], 
	       [-theta, -theta * theta, theta * theta]]
	return np.array(ret)

def gen_M_I():
	ret = np.zeros((3, 3))
	ret[0, 0] = 1
	return ret

def gen_M_mu_A(alpha, mu):
	ret = [[0, -1/2, 0], 
	       [-1/2, -(1+alpha*mu), 1],
	       [0, 1, 0]]
	return np.array(ret)

def gen_M_beta_B(alpha, beta):
	ret = [[-beta/alpha, 0, beta/alpha + 1/2], 
	       [0, 0, 0],
	       [beta/alpha + 1/2, 0, -beta/alpha-1]]
	return np.array(ret)


def Primal(alpha, theta, mu, beta):
	G = cp.Variable((3, 3), symmetric = True)
	
	M_O = gen_M_O(theta)
	M_I = gen_M_I()
	M_mu_A = gen_M_mu_A(alpha, mu)
	M_beta_B = gen_M_beta_B(alpha, beta)

	constraints = []
	constraints.append(G >> 0)
	constraints.append(cp.trace(M_I @ G) == 1)
	constraints.append(cp.trace(M_mu_A @ G) >= 0)
	constraints.append(cp.trace(M_beta_B @ G) >= 0)

	objective = cp.Maximize(cp.trace(M_O @ G))

	problem = cp.Problem(objective, constraints)
	problem.solve()

	return math.sqrt(problem.value)

def Dual(alpha, theta, mu, beta):
	lambda_mu_A = cp.Variable()
	lambda_beta_B = cp.Variable()
	rho_sq = cp.Variable()

	M_O = gen_M_O(theta)
	M_I = gen_M_I()
	M_mu_A = gen_M_mu_A(alpha, mu)
	M_beta_B = gen_M_beta_B(alpha, beta)

	constraints = []
	constraints.append(lambda_mu_A >= 0)
	constraints.append(lambda_beta_B >= 0)
	constraints.append(-M_O - lambda_mu_A * M_mu_A - lambda_beta_B * M_beta_B + rho_sq * M_I >> 0)

	objective = cp.Minimize(rho_sq)

	problem = cp.Problem(objective, constraints)
	problem.solve()

	return math.sqrt(problem.value)

eps = 0.0005
checked = [0, 0, 0, 0, 0]

for _ in range(300):
	alpha = np.random.uniform(low = 0.5, high = 2.0)
	theta = np.random.uniform(low = 0.05, high = 1.95)
	mu = np.random.uniform(low = 0.1, high = 3.9)
	beta = np.random.uniform(low = 0.1, high = 3.9)

	assert alpha > 0 and theta > 0
	assert mu > 0 and beta > 0 

	val_1 = Primal(alpha, theta, mu, beta)
	val_2 = Dual(alpha, theta, mu, beta)
	whi, val_3 = THM_4_1_unscaled(alpha, theta, mu, beta)

	checked[whi - 1] += 1

	if abs(val_2 - val_1) > eps or abs(val_3 - val_1) > eps or abs(val_3 - val_2) > eps:
		print(val_1, val_2, val_3)
		print("Failed at Case", whi)
		print("Parameters:", alpha, theta, mu, beta)
		break

print("Checked Cases")
for i in range(5):
	print("Case #{} Verified {} Times!".format(i + 1, checked[i]))

print("Finished!")