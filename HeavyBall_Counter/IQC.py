import numpy as np 
import cvxpy as cp
import warnings
import matplotlib.pyplot as plt
from tqdm import tqdm

warnings.filterwarnings("ignore")

def dynamic(r, s, L):
    A = np.matrix([[1 + r, -r], [1, 0]])
    B = np.matrix([[-s], [0]])
    C = np.matrix([1, 0])
    D = np.matrix([0])
    return A, B, C, D 

def IQC_weighted_off_by_one(mu, L):
    A_phi = np.matrix([0])
    B_y_phi = np.matrix([-L])
    B_u_phi = np.matrix([1])
    return A_phi, B_y_phi, B_u_phi

def attempt(rho_1, r, s, mu, L):
    A, B, C, D = dynamic(r, s, L)
    A_phi, B_y_phi, B_u_phi = IQC_weighted_off_by_one(mu, L)
    A_hat = np.block([[A, np.zeros((2, 1))], [B_y_phi @ C, A_phi]])
    B_hat = np.block([[B], [B_u_phi]])

    P = cp.Variable((3, 3), PSD = True)
    lam_1 = cp.Variable()  
    lam_2 = cp.Variable()
    Target1 = cp.bmat([[A_hat.transpose() @ P @ A_hat - rho_1 * P, 
                        A_hat.transpose() @ P @ B_hat],
                        [B_hat.transpose() @ P @ A_hat, B_hat.transpose() @ P @ B_hat]
    ])
    
    MM = np.zeros((4, 4))
    CC = np.zeros((4, 4))
    MM[0, 0] = - 2 * mu * L 
    MM[0, 3] = mu + L 
    MM[3, 0] = mu + L
    MM[3, 3] = - 2
    CC[0, 2] = - mu 
    CC[2, 0] = - mu 
    CC[2, 3] = 1
    CC[3, 2] = 1
    
    constraints = []
    constraints.append(P >> np.eye(3))
    constraints.append(P == P.T)
    constraints.append(lam_1 >= 0)
    constraints.append(lam_2 >= 0)
    constraints.append(rho_1 * lam_2 >= lam_1)
    constraints.append(Target1 + (CC * lam_1 + MM * lam_2) << 0)
    prob = cp.Problem(cp.Minimize(1), constraints = constraints)
    prob.solve()

    if prob.status == "infeasible" or prob.status == "infeasible_inaccurate":
        return False # failure
    else:
        return True
    
def get_result(r, s, mu, L):
    lef, rig, eps, best = 0.5, 5, 10 ** -8, 10
    while rig - lef > eps:
        mid = (lef + rig) / 2
        if attempt(mid, r, s, mu, L) == True:
            best, rig = mid, mid 
        else:
            lef = mid 
    return best

x = []
y = []
z = []

for i in tqdm(range(0, 100)):
    kappa = 2 + 0.3 * i
    mu = 1 / kappa
    L = 1

    r = ((np.sqrt(kappa) - 1) / (np.sqrt(kappa) + 1)) ** 2
    s = 4 / ((np.sqrt(L) + np.sqrt(mu)) ** 2)

    x.append(kappa)
    y.append(get_result(r, s, mu, L))
    z.append(1)

plt.plot(x, y, 'b')
plt.plot(x, z, 'r')
plt.savefig("resultHeavy.PNG")
