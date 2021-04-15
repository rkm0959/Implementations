from sklearn.datasets import load_svmlight_file
import numpy as np
import scipy.sparse
from tqdm import tqdm
import matplotlib.pyplot as plt 

print("[+] Loading Data and Converting")

data, whi = load_svmlight_file("real-sim")

# Catalyst-SAGA with Criterion C3 vs non-accelerated SAGA

n, m = data.get_shape()
CYC = 200

values = data.data
indices = data.indices
indptr = data.indptr

# preproccess : unit length

a = []
for i in range(0, n):
    row = []
    tot = 0
    assert whi[i] == -1 or whi[i] == 1
    for j in range(indptr[i], indptr[i+1]):
        tot += values[j] * values[j]
    dv = np.sqrt(tot)
    if abs(dv) < 10 ** -6:
        a.append(row)
        continue
    for j in range(indptr[i], indptr[i+1]):
        row.append((indices[j], values[j] / dv))
    a.append(row)


L = 0.25 
reg = 1 / (256 * n)

# alpha_k -> alpha_{k+1}
def update_alpha(alpha, q):
    A = 1
    B = alpha * alpha - q
    C = -alpha * alpha
    x_1 = (-B + np.sqrt(B * B - 4 * A * C)) / (2 * A)
    x_2 = (-B - np.sqrt(B * B - 4 * A * C)) / (2 * A)   
    if 0 < x_1 < 1:
        return x_1
    return x_2

# dot product <a_i, x>
def dot_prod(i, x):
    ret = 0
    for loc, val in a[i]:
        ret += val * x[loc]
    return ret 

# gradient of f_i(x) but a bit convoluted...
def grad1(i, x, y, rat):
    val = rat * dot_prod(i, x) + dot_prod(i, y)
    T = np.exp(-whi[i] * val)
    return - whi[i] * T / (1 + T)

# gradient of f_i(x)
def grad2(i, x, rat):
    val = rat * dot_prod(i, x)
    T = np.exp(-whi[i] * val)
    return - whi[i] * T / (1 + T)

# 1/n * f_i(x)
def func(i, x):
    val = dot_prod(i, x)
    return 1 / n * np.log(1 + np.exp(-whi[i] * val))

# reg / 2 ||x||^2 + sum 1/n * f_i(x)
def f(x):
    ret = reg / 2 * (np.linalg.norm(x) ** 2)
    for i in range(0, n):
        ret += func(i, x)
    return ret 

# f(x) + kappa / 2 ||x-y||^2
def h(x, y, kappa):
    ret = f(x)
    ret += kappa / 2 * (np.linalg.norm(x - y) ** 2)
    return ret

# prediction accuracy
def prediction(x):
    correct = 0
    for i in range(0, n):
        val = dot_prod(i, x)
        if np.sign(val) == whi[i]:
            correct += 1
    return correct / n 

# Inner Iteration in Catalyst-SAGA 
def SAGA_InnerIteration(REG, yk, diff, step):
    d = np.zeros(m)
    y = [0] * n
    S = []
    V = [0] * m
    x, z = yk, yk
    rat, cnt = 1, 0

    # one pass strategy (C3)
    for i in range(0, n):
        cnt += 1
        for loc, val in a[i]:
            z[loc] -= d[loc] * ((S[i-1] if i >= 1 else 0) - (S[V[loc]-1] if V[loc] >= 1 else 0))
            V[loc] = i
        vv = grad1(i, z, diff, rat)
        for loc, val in a[i]:
            d[loc] += val * (vv - y[i])
        rat = rat * (1 - step * REG)
        for loc, val in a[i]:
            z[loc] -= step / rat * (cnt-1) / cnt * (vv - y[i]) * val
        cc = S[i-1] if i >= 1 else 0
        S.append(cc + step / rat / cnt)
        y[i] = vv
    
    # final updates
    for i in range(0, m):
        z[i] -= d[i] * (S[n-1] - (S[V[i]-1] if V[i] >= 1 else 0))
        x[i] = rat * z[i]
    
    return x + diff

def CatalystSAGA():
    print("[+] Setting Parameters")

    x = np.zeros(m)
    y = np.zeros(m)

    kappa = (L - reg) / (n + 1) - reg
    q = reg / (reg + kappa) 
    alpha = np.sqrt(q)
    gamma = (reg + kappa) / (2 * ((reg + kappa) * n + L + kappa))

    cand_1 = x
    cand_2 = x

    NUM_ITER = []
    FUNC = []
    PRED = []
    
    print("[+] Begin Catalyst_SAGA")

    for i in tqdm(range(0, CYC + 1)):
        NUM_ITER.append(i)
        FUNC.append(f(x))
        PRED.append(prediction(x))
        if i == CYC:
            break
        # choose initial value
        val_1 = h(cand_1, y, kappa)
        val_2 = h(cand_2, y, kappa)
        # run inner iteration of Catalyst_SAGA
        if val_1 < val_2:
            st = kappa * cand_1 / (kappa + reg)
            x_new = SAGA_InnerIteration(reg + kappa, cand_1 - st, st, gamma)
        else:
            st = kappa * cand_2 / (kappa + reg)
            x_new = SAGA_InnerIteration(reg + kappa, cand_2 - st, st, gamma) 
        # update alpha
        alpha_new = update_alpha(alpha, q)
        # extrapolation & updates
        beta = (alpha * (1 - alpha)) / (alpha * alpha + alpha_new)
        cand_2 = - kappa / (kappa + reg) * y
        y = x_new + beta * (x_new - x)
        cand_2 += kappa / (kappa + reg) * y
        x, alpha = x_new, alpha_new
        cand_2 += x
        cand_1 = x

    return NUM_ITER, FUNC, PRED

def SAGA():
    print("[+] Setting Parameters")

    x = np.zeros(m)
    y = np.zeros(m)
    z = np.zeros(m)
    step = reg / (2 * (reg * n + L))

    d = np.zeros(m)
    y = [0] * n
    S = []
    V = [0] * m
    rat, cnt = 1, 0

    NUM_ITER = []
    FUNC = []
    PRED = []

    print("[+] Begin SAGA")

    NUM_ITER.append(0)
    FUNC.append(f(x))
    PRED.append(prediction(x))
    for i in tqdm(range(0, CYC*n)):
        if i < n:
            cnt += 1
        for loc, val in a[i % n]:
            z[loc] -= d[loc] * ((S[i-1] if i >= 1 else 0) - (S[V[loc]-1] if V[loc] >= 1 else 0))
            V[loc] = i
        vv = grad2(i % n, z, rat)
        for loc, val in a[i % n]:
            d[loc] += val * (vv - y[i % n])
        rat = rat * (1 - step * reg)
        for loc, val in a[i % n]:
            z[loc] -= step / rat * (cnt-1) / cnt * (vv - y[i % n]) * val
        cc = S[i-1] if i >= 1 else 0
        S.append(cc + step / rat / cnt)
        y[i % n] = vv
        
        # one pass complete
        if i % n == n - 1:
            # full update on x
            for j in range(0, m):
                z[j] -= d[j] * (S[i] - (S[V[j]-1] if V[j] >= 1 else 0))
                x[j] = rat * z[j]
                V[j] = i + 1
            NUM_ITER.append((i+1)//n)
            FUNC.append(f(x))
            PRED.append(prediction(x))
    
    return NUM_ITER, FUNC, PRED

X1, FUNC_1, PRED_1 = CatalystSAGA()
X2, FUNC_2, PRED_2 = SAGA()

assert X1 == X2 

plt.xlabel('Number of Passes')
plt.ylabel('Function Value')
plt.plot(X1, FUNC_1, color = 'r', label = 'CatalystSAGA')
plt.plot(X1, FUNC_2, color = 'b', label = 'SAGA')
plt.legend()
plt.savefig("result_real-sim_loss.PNG")

plt.clf()

plt.xlabel('Number of Passes')
plt.ylabel('Prediction Accuracy')
plt.plot(X1, PRED_1, color = 'r', label = 'CatalystSAGA')
plt.plot(X1, PRED_2, color = 'b', label = 'SAGA')
plt.legend()
plt.savefig("result_real-sim_pred.PNG")