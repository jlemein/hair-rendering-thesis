from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

def sampleSphericalGaussian(variance, u):
    return variance * np.log(np.exp(1.0 / variance) - 2.0 * u * np.sinh(1.0 / variance))

def sampleThetaI(thetaCone, variance, u1, u2):
    e1 = sampleSphericalGaussian(variance, u1)
    thetaAccent = 0.5 * np.pi - thetaCone
    return np.arcsin(e1 * np.cos(thetaAccent) + np.sqrt(1.0 - e1 * e1) * np.cos(2.0 * np.pi * u2) * np.sin(thetaAccent))


# variance = 0.2
# x = np.linspace(0.0, 0.99)
# z = np.linspace(0.0, 0.99)
# X, Z = np.meshgrid(x, z)

# fig = plt.figure()

# # Plot the surface.
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Z, aa(thetaCone, variance, X**2, Z**2), cmap=cm.coolwarm, linewidth=0, antialiased=True)
# fig.colorbar(surf, shrink=0.5, aspect=5)


variance = 0.2
thetaO = -np.pi * .4
alpha = (5.0 / 180.0) * np.pi
thetaCone = -thetaO + alpha

SIZE = 100000
u1 = np.random.sample(SIZE)
u2 = np.random.sample(SIZE)

print("u1: " + str(u1))
print("u2: " + str(u2))

samples = np.empty(SIZE)
print("samples before: " + str(samples))
for idx in range(SIZE):
    samples[idx] = sampleThetaI(thetaCone, variance, u1[idx], u2[idx])

print("samples after: " + str(samples))

minThetaI = np.amin(samples)
maxThetaI = np.amax(samples)
print("lowest thetaI: " + str(minThetaI))
print("highest thetaI: " + str(maxThetaI))

plt.hist(samples, bins=90)

print("Sum: " + str(np.sum(samples)))


# Plot 2d plot
#thetaO = aa(thetaCone, variance, x, 1.0)
#thetaD = .5 * (thetaI + thetaO);
#plt.plot(x, thetaD, label="variance = 0.2")


plt.legend()
plt.show()