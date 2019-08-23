import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import special
import scipy


def integrate(fn, start, end):
    domain = end - start
    x = start + domain * np.random.rand(10000)
    return np.average(fn(x)) * domain


def normGaussian(x, mean, variance):
    frac = np.square(x - mean)/(2.0*variance)
    return np.exp(-frac) / np.sqrt(2.0 * np.pi * variance)


def csch(x):
    return 2.0 / (np.exp(x) - np.exp(-x))

# modified bessel function of the first kind I0


def factorial(n):
    result = 1
    while (n > 0):
        result *= n
        n -= 1

    return result


print(factorial(1))
print(factorial(2))
print(factorial(3))
print(factorial(4))
print(factorial(5))


def I0(z):
    sum = 1.0
    nom = 0.25*z*z
    for k in range(1, 10):
        fact = np.prod(np.arange(1.0, k+1))
        sum += np.power(nom, k) / (fact*fact)

    return sum


def I00(X):

    P1 = 1.0
    P2 = 3.5156229
    P3 = 3.0899424
    P4 = 1.2067492
    P5 = 0.2659732
    P6 = 0.360768e-1
    P7 = 0.45813e-2
    Q1 = 0.39894228
    Q2 = 0.1328592e-1
    Q3 = 0.225319e-2
    Q4 = -0.157565e-2
    Q5 = 0.916281e-2
    Q6 = -0.2057706e-1
    Q7 = 0.2635537e-1
    Q8 = -0.1647633e-1
    Q9 = 0.392377e-2

    if np.absolute(X) < 3.75:
        Y = (X/3.75)*(X/3.75)
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    else:
        AX = np.absolute(X)
        Y = 3.75/AX
        BX = np.exp(AX)/np.sqrt(AX)
        AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        return (AX*BX)


def M_scipy(v, theta_i, theta_r):
    a = csch(1.0/v) / (2.0*v)
    print("Msc::a = ", a)
    b = (np.sin(-theta_i) * np.sin(theta_r))/v
    print("Msc::b = ", b)
    print("Msc:: arg = ", np.cos(-theta_i) * np.cos(theta_r)/v)
    bss = special.iv(1, np.cos(-theta_i) * np.cos(theta_r)/v)
    return a * np.exp(b) * bss


aa = np.vectorize(I00)


def M(v, theta_i, theta_r):
    a = csch(1.0/v) / (2.0*v)
    print("M::a = ", a)

    b = (np.sin(-theta_i) * np.sin(theta_r))/v
    print("M::b = ", b)
    bss = aa(np.cos(-theta_i) * np.cos(theta_r)/v)
    return a * np.exp(b) * bss


beta = 0.05
variance = beta*beta

x = np.linspace(-.5*np.pi, .5*np.pi, 1000)


def _normGaussian(x): return normGaussian(x, 0, 0.4)


thetaR = 0.0


def _M(x): return M(variance, x, thetaR)


def _M_scipy(x): return M_scipy(variance, x, thetaR)


print("Integrating normalize gaussian: ",
      integrate(_normGaussian, -.5*np.pi, .5*np.pi))
print("Integrating M: ", integrate(_M, -.5*np.pi, .5*np.pi))
print("Integrating Mscopy: ", integrate(_M_scipy, -.5*np.pi, .5*np.pi))

print("Msc(0): ", M_scipy(variance, 0.0, 0.0))
print("M(0): ", M(variance, 0.0, 0.0))


# Plot results
# fig = plt.figure()
# plot = fig.add_subplot(2, 1, 1)
# #plot.plot(x, I0(x), label="my bessel")
# #plot.plot(x, special.iv(0, x), label="bessel (scipy)")
# plot.plot(x, special.iv(0, x) - I0(x), label="should be zero")
# #plot.plot(x, _M(x), label="M")
# #plot.plot(x, _M_scipy(x), label="M scipy")
# plot.legend()

# plot = fig.add_subplot(2, 1, 2)
# plot.plot(x, special.iv(0, x), label="bessel (scipy)")
# plot.legend()


# # Plot 3d
# # thetai = np.linspace(-.5*np.pi, .5*np.pi, 100)
# # thetar = np.linspace(-.5*np.pi, .5*np.pi, 100)
# # X, Y = np.meshgrid(thetai, thetar)

# # plot = fig.add_subplot(2, 1, 2, projection='3d')
# # plot.plot_surface(X, Y, M(variance, X**2, Y**2), cmap="viridis")

# plt.show()
