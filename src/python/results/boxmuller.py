from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


def BoxMuller(u1, u2):
    return np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2)

def gaussian(x, mu, sigma):
    a = 1.0 / sigma * np.sqrt(2.0 * np.pi)
    nom = (x - mu) / sigma;

    return a * np.exp(-.5 * nom * nom)


#
# Generate samples from Box muller
#
SIZE = 100000
u1 = np.random.sample(SIZE)
u2 = np.random.sample(SIZE)
samples = np.empty(SIZE)

for idx in range(SIZE):
    samples[idx] = BoxMuller(u1[idx], u2[idx])

#
# Plot a histogram for the distribution of samples (should be gaussian like)
# #plt.hist(samples, bins=90)

#
# convert histogram to a normalized function (and compare it with a reference gaussian)
#
nBins = 10
binRanges = np.linspace(-np.pi/2, np.pi/2, nBins+1)
histdata = np.histogram(samples, binRanges)[0]
print("Histdata: ", str(histdata))
print("Sum of histdata = ", str(np.sum(histdata)))

# normalize histogram data
histdata = np.true_divide(histdata, np.sum(histdata))
print("normalized Histdata: ", str(histdata))
print("sum histdata: " + str(np.sum(histdata)))
plt.plot(np.linspace(-np.pi/2, np.pi/2, nBins), histdata, label="Distribution of box muller samples")

# Plot refence gaussian
x = np.linspace(-np.pi/2, np.pi/2)
#plt.plot(x, gaussian(x, 0, 1.0), label="Reference normalized gaussian")

print("gaussian sum: " + str(np.pi*np.sum(gaussian(x, 0, 1.0))/))

plt.legend()
plt.show()