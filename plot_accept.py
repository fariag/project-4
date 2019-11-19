from  matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm

data = np.loadtxt("output_accept1_1_spin=20.txt")
x = data[:, 0]
y = data[:, 1]

data1 = np.loadtxt("output_accept1_2_spin=20.txt")
x1 = data1[:, 0]
y1 = data1[:, 1]

plt.figure(1)

plt.semilogx(x, y, c="b")
plt.semilogx(x1, y1, c="r")
plt.grid()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Accepted Configuration / (MC cycles*spins^2)")
plt.title("Accepted Configuration, Temperature=1.0")
plt.legend(["Random config", "Ordered config"])
plt.savefig("4c_mcs_t=1")
plt.show()
