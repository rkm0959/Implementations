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

def gen_Base():
	ret = np.zeros((4, 4))
	ret[0, 3] = 1
	ret[3, 0] = 1
	ret[3, 3] = 1
	return ret

def gen_M_theta():
	ret = np.zeros((4, 4))
	ret[1, 3] = 1
	ret[2, 3] = -1
	ret[3, 1] = 1
	ret[3, 2] = -1
	return ret

def gen_M_I():
	ret = np.zeros((4, 4))
	ret[0, 0] = 1
	return ret

def gen_M_mu_A(alpha, mu):
	ret = [[0, -1/2, 0, 0], 
	       [-1/2, -(1+alpha*mu), 1, 0],
	       [0, 1, 0, 0],
	       [0, 0, 0, 0]]
	return np.array(ret)

def gen_M_beta_B(alpha, beta):
	ret = [[-beta/alpha, 0, beta/alpha + 1/2, 0], 
	       [0, 0, 0, 0],
	       [beta/alpha + 1/2, 0, -beta/alpha-1, 0], 
	       [0, 0, 0, 0]]
	return np.array(ret)


def opt_val(alpha, mu, beta, retrieve = False):
	lambda_mu_A = cp.Variable()
	lambda_beta_B = cp.Variable()
	theta = cp.Variable()
	rho_sq = cp.Variable()

	Base = gen_Base()
	M_theta = gen_M_theta()
	M_I = gen_M_I()
	M_mu_A = gen_M_mu_A(alpha, mu)
	M_beta_B = gen_M_beta_B(alpha, beta)

	constraints = []
	constraints.append(lambda_mu_A >= 0)
	constraints.append(lambda_beta_B >= 0)
	constraints.append(theta >= 0)
	constraints.append(Base + theta * M_theta - lambda_mu_A * M_mu_A - lambda_beta_B * M_beta_B + rho_sq * M_I >> 0)

	objective = cp.Minimize(rho_sq)

	problem = cp.Problem(objective, constraints)
	problem.solve()

	if retrieve == True:
		return math.sqrt(problem.value), theta.value
	else:
		return math.sqrt(problem.value)

mu = 0.53
beta = 1.35

alpha_L = 0.05
alpha_R = 3.95

print("Begin Ternary Search")
while alpha_R - alpha_L >= 0.0001:
	left = alpha_L + (alpha_R - alpha_L) / 3
	right = alpha_R - (alpha_R - alpha_L) / 3
	if opt_val(left, mu, beta) < opt_val(right, mu, beta):
		alpha_R = right
	else:
		alpha_L = left

opt_alpha = (alpha_L + alpha_R) / 2
opt_rho, opt_theta = opt_val(opt_alpha, mu, beta, True)

print("Optimal alpha, theta", opt_alpha, opt_theta)
print("Optimal contraction factor", opt_rho)
print("Check with Theorem 4.1.", THM_4_1_unscaled(opt_alpha, opt_theta, mu, beta))

print("Begin Graphing")
X = np.linspace(0.25, 3.75, 50)
Y = np.array([opt_val(alpha, mu, beta) for alpha in X])

plt.plot(X, Y)
plt.show()

plt.savefig("Graph.PNG")

print("Finished!")