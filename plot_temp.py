from  matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm

data100 = np.loadtxt("data_l=100.txt")
x = data100[:, 0]
y = data100[:, 1]

data80 = np.loadtxt("data_l=80.txt")
x1 = data80[:, 0]
y1 = data80[:, 1]

data60 = np.loadtxt("data_l=60.txt")
x2 = data60[:, 0]
y2 = data60[:, 1]

data40 = np.loadtxt("data_l=40.txt")
x3 = data40[:, 0]
y3 = data40[:, 1]

plt.figure(1) # Mean Energy
plt.plot(x, y, c="r")
plt.plot(x1, y1, c="g")
plt.plot(x2, y2, c="b")
plt.plot(x3, y3, c="c")
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("Mean Energy")
plt.legend(["L=100", "L=80", "L=60", "L=40"])
plt.title("Mean Energy over Temperature ")
plt.savefig("4e_Energy_plot")

y = data100[:, 2]
idx = np.argmax(y)
print(x[idx])
y1 = data80[:, 2]

y2 = data60[:, 2]

y3 = data40[:, 2]

plt.figure(2) # Specific Heat
plt.plot(x[idx], y[idx], "o")
plt.plot(x, y, c="r")
plt.plot(x1, y1, c="g")
plt.plot(x2, y2, c="b")
plt.plot(x3, y3, c="c")
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("Specific Heat")
plt.legend(["L=100", "L=80", "L=60", "L=40"])
plt.title("Specific Heat over Temperature ")
plt.savefig("4e_Cv_plot")

y = data100[:, 3]

y1 = data80[:, 3]

y2 = data60[:, 3]

y3 = data40[:, 3]

plt.figure(3) # Magnetization

plt.plot(x, y, c="r")
plt.plot(x1, y1, c="g")
plt.plot(x2, y2, c="b")
plt.plot(x3, y3, c="c")
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("Magnetization")
plt.legend(["L=100", "L=80", "L=60", "L=40"])
plt.title("Magnetization over Temperature ")
plt.savefig("4e_mag_plot")

y = data100[:, 4]

y1 = data80[:, 4]

y2 = data60[:, 4]

y3 = data40[:, 4]

plt.figure(4) # Susceptibility

plt.plot(x, y, c="r")
plt.plot(x1, y1, c="g")
plt.plot(x2, y2, c="b")
plt.plot(x3, y3, c="c")
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")
plt.legend(["L=100", "L=80", "L=60", "L=40"])
plt.title("Susceptibility over Temperature ")
plt.savefig("4e_sus_plot")

y = data100[:, 5]

y1 = data80[:, 5]

y2 = data60[:, 5]

y3 = data40[:, 5]

plt.figure(5) # |Magnetization|

plt.plot(x, y, c="r")
plt.plot(x1, y1, c="g")
plt.plot(x2, y2, c="b")
plt.plot(x3, y3, c="c")
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("|Magnetization|")
plt.legend(["L=100", "L=80", "L=60", "L=40"])
plt.title("|Magnetization| over Temperature ")
plt.savefig("4e_absmag_plot")
plt.show()
