import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MultipleLocator, FuncFormatter

# Make some fake data.
uniform_vs_groundtruth = np.array([0.019034,
                                   0.009072,
                                   0.004352,
                                   0.001842,
                                   0.001004,
                                   0.000357,
                                   0.000156,
                                   0.000054,
                                   0.000017,
                                   0.000000])

deon_vs_groundtruth = np.array([0.010120,
                                0.004829,
                                0.002440,
                                0.001119,
                                0.000760,
                                0.000415,
                                0.000321,
                                0.000282,
                                0.000264,
                                0.000259])

deon_vs_deon = np.array([0.009840,
                         0.004560,
                         0.002169,
                         0.000860,
                         0.000493,
                         0.000153,
                         0.000058,
                         0.000022,
                         0.000007,
                         0.000000])

#samples_per_pixel = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
samples_per_pixel = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Create plots with pre-defined labels.
fig, ax = plt.subplots()
x = np.logspace(0, np.log10(512), 10)
ax.plot(samples_per_pixel, uniform_vs_groundtruth,
        label='Uniform Sampling', marker='x')
ax.plot(samples_per_pixel, deon_vs_groundtruth,
        label='Importance Sampling', marker='x')
# ax.plot(samples_per_pixel, deon_vs_deon,
#       label='Importance Sampling (relative)', marker='x')

plt.xscale('log')
plt.xlabel("Samples per pixel")
plt.ylabel("Variance (compared with uniform sampling 512 spp)")
plt.grid(True)
# ax.semilogx(range(512))
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.xaxis.set_minor_formatter(ScalarFormatter())
# ax.xaxis.set_major_locator(plt.MaxNLocator(12))
ax.xaxis.set_major_locator(MultipleLocator(1))
# Format the ticklabel to be 2 raised to the power of `x`
ax.xaxis.set_major_formatter(FuncFormatter(
    lambda x, pos: int(np.power(2, x-1))))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)

legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')

plt.show()
