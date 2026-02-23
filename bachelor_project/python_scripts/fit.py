import sys
sys.path.insert(0,'..')
from prec import pprec
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

x, dx, y, dy = np.loadtxt("vinv", unpack=True)

dx = (dx**2 + .03**2*x**2)**.5
dy = (dy**2 + .03**2*y**2)**.5

xx = np.linspace(np.min(x), np.max(x), 10)

plt.minorticks_on()
plt.grid()

plt.errorbar(x, y, dy, dx, fmt="k.")

f = lambda x, a, b: a + b * x
p0 = [.0, 10.]
##plt.plot(xx, f(xx, *p0), "k")

sigma = dy
popt, pcov = curve_fit(f, x, y, p0, sigma, absolute_sigma=False)

yy = f(xx, *popt)
yr = (y - f(x, *popt)) / sigma

print("p:\n{}".format(popt))
print("cov:\n{}".format(pcov))
print("dp:\n{}".format(np.sqrt(pcov.diagonal())))
print("chisq: {}".format(sum(yr**2)))
print("ndof: {}".format(x.size-len(p0)))

plt.plot(xx, yy, color="k")

plt.xlabel(r"$V_\mathrm{in}$")
plt.ylabel(r"$V_\mathrm{out}$")

#plt.figure()
#plt.minorticks_on()
#plt.grid()
#plt.errorbar(x0, yr0, fmt="--b.")
#plt.errorbar(x, yr, fmt="--r.")
#plt.plot(x, np.full(x.size, 0.0), "k")

#for i in range(x.size): print("{}\t&\t{} \\\\".format(pprec(x[i], dx[i]), pprec(y[i], dy[i])))
#for i in range(x.size): print(x[i], dx[i], y[i], dy[i])

plt.show()

