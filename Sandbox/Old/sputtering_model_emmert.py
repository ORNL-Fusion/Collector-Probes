import numpy as np
import scipy.integrate as integrate
import scipy.optimize
from math import *
import matplotlib.pyplot as plt

# Boltzmann constant in ev * K^-1
k = 8.617 * 10**(-5)
# Mass of deuterium in eV s^2 m^-2
M = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
# Mass of electron in eV s^2 m^-2
m = 0.511 * 10**6 / ((3*10**8)**2.0)
# Dimensionless charge variable
Z = 1.0
# Dimensionless temperature variable
tau = 1.0



def D(x):
    def integrand(t):
        return np.exp(t**2)
    int_part = integrate.quad(integrand, 0, x)
    result = np.exp(-x**2) * int_part[0]
    return result

def F(psi, Z, tau):
    result = erf(sqrt(Z*tau*psi)) + 2 / sqrt(pi*Z*tau) \
        * np.exp(-(1+Z*tau)*psi) * D(sqrt(psi))
    return result

def H(psi, Z, tau, Te, Ti, M, psi0, vx):
    v0 = sqrt(2/M * Z*k*Te*psi)
    if vx > v0:
        #print("H case 1")
        return 1 + F(psi, Z, tau)
    elif vx > 0 and vx < v0:
        #print("H case 2")
        return 1 + F(psi, Z, tau) - 2*F(psi0, Z, tau)
    else:
        #print("H case 3")
        return 1 - F(psi, Z, tau)

def f(vx, psi, Z, tau, Te, Ti, M, n0):
    psi0 = psi - M*vx**2 / (2*k*Ti*Z*tau)
    #print("psi0 = " + str(psi0))
    #print(H(psi, Z, tau, Te, Ti, M, psi0, vx))
    return (n0/Z) * sqrt(M/(2*pi*k*Ti)) * np.exp(Z*tau*psi0) \
        * H(psi, Z, tau, Te, Ti, M, psi0, vx)

def solve_for_psi1(Z, tau):
    def func (psi1, Z, tau):
        return 2.0 / sqrt(pi*Z*tau) * np.exp(-(1+Z*tau) * psi1) * D(sqrt(psi1)) \
            + erf(sqrt(Z*tau*psi1))  - 1.0
    minimum = scipy.optimize.minimize(func, x0=0.5, args=(Z, tau), method="SLSQP", bounds=(0.0, None))
    return minimum[0]

psi1 = solve_for_psi1(Z, tau)
Ti=25
Te=25
vxs = np.arange(-10, 10, 0.1) * sqrt(2*k*25 / M)

print("***************")
print("psi = 0")
print("****************")
fs = np.array([])
for v in vxs:
    tmp_f = f(v, 0, 1, 1, Te, Ti, M=M, n0=10**18)
    #print("f = " + str(tmp_f))
    if tmp_f < 0.0:
        fs = np.append(fs, 0.0)
    else:
        fs = np.append(fs, tmp_f)
fs = fs / np.max(fs)
plt.plot(vxs / np.max(vxs), fs)

print("****************")
print("psi = 0.5*psi1")
print("****************")
fs = np.array([])
for v in vxs:
    tmp_f = f(v, 0.5*psi1, 1, 1, Te, Ti, M=M, n0=10**18)
    #print("f = " + str(tmp_f))
    if tmp_f < 0.0:
        fs = np.append(fs, 0.0)
    else:
        fs = np.append(fs, tmp_f)
fs = fs / np.max(fs)
plt.plot(vxs / np.max(vxs), fs)

print("****************")
print("psi = 1.0*psi1")
print("****************")
fs = np.array([])
for v in vxs:
    tmp_f = f(v, psi1, 1, 1, Te, Ti, M=M, n0=10**18)
    #print("f = " + str(tmp_f))
    if tmp_f <= 0.0:
        fs = np.append(fs, 0.0)
    else:
        fs = np.append(fs, tmp_f)
fs = fs / np.max(fs)
plt.plot(vxs / np.max(vxs), fs)

#plt.ylim([0, 1])
#plt.xlim([0, 0.75])
plt.show()
