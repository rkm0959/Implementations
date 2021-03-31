import matplotlib.pyplot as plt
import numpy as np

def grad(x):
    if x < 1:
        return 25 * x
    if 1 <= x < 2:
        return x + 24
    if x >= 2:
        return 25 * x - 24

mu = 1
L = 25
kappa = 25

alpha = 4 / ((np.sqrt(L) + np.sqrt(mu)) ** 2)
beta = ((np.sqrt(kappa) - 1) / (np.sqrt(kappa) + 1)) ** 2

x_0 = - 973 / 300
x_1 = x_0 - 2 / (L + mu) * grad(x_0)
x = [x_0, x_1]

for i in range(2, 30):
    x_2 = x_1 - alpha * grad(x_1) + beta * (x_1 - x_0)
    x.append(x_2)
    x_0, x_1 = x_1, x_2

plt.plot(x)
plt.savefig("counterexample.PNG")
