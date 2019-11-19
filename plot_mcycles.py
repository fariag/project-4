from  matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm

# data = np.loadtxt("output_energy_2_r_spin=20.txt")
# x = data[:, 0]
# y = data[:, 1]

data1 = np.loadtxt("output_energy_spin=2.txt")
x1 = data1[:, 0]
y1 = data1[:, 1]

plt.figure(1)

# plt.semilogx(x, y, c="b")
plt.semilogx(x1, y1, c="r")
plt.grid()
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Mean energy")
plt.title("Mean energy for L=2, Temperature=1.0")
plt.legend(["Ordered config"])
plt.savefig("4bPlot_energy_2")
plt.show()
