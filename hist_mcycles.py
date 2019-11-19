from  matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm



data = np.loadtxt("Poutput_1_spin=20.txt")
y1 = data
y1_l = y1[10000:]
data2 = np.loadtxt("output_2.4_spin=20.txt")
y2 = data2
y2_l = y2[10000:]
plt.figure(1)
weights1 = np.ones_like(y1_l)/len(y1_l)
plt.hist(y1_l, bins=10, weights=weights1)
plt.grid()
plt.xlabel("Energy levels")
plt.ylabel("P(E)")
plt.title("Probability for L=20 at temperature=1.0")
plt.savefig("4dPlot_1")

plt.figure(2)
weights2 = np.ones_like(y2_l)/len(y2_l)
plt.hist(y2_l, bins=20, weights=weights2)

plt.grid()
plt.xlabel("Energy levels")
plt.ylabel("P(E)")
plt.title("Probability for L=20 at temperature=2.4")
plt.savefig("4dPlot_2")
plt.show()
