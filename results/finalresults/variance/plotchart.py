import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MultipleLocator, FuncFormatter

# Make some fake data.
uniform_vs_groundtruth = np.array([0.004720, 0.001962, 0.000909, 0.000459,
                                   0.000275, 0.000173, 0.000130, 0.000105, 0.000003, 0.000000])
deon_vs_groundtruth = np.array([0.002616, 0.001177, 0.000689, 0.000481,
                                0.000406, 0.000371, 0.000355, 0.000346, 0.000257, 0.000257])
deon_vs_deon = np.array([0.002374, 0.000930, 0.000432, 0.000223,
                         0.000146, 0.000109, 0.000092, 0.000084, 0.000001, 0.000000])

#samples_per_pixel = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
samples_per_pixel = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Create plots with pre-defined labels.
fig, ax = plt.subplots()
x = np.logspace(0, np.log10(512), 10)
ax.plot(samples_per_pixel, uniform_vs_groundtruth,
        label='Uniform Sampling', marker='x')
ax.plot(samples_per_pixel, deon_vs_groundtruth,
        label='Importance Sampling', marker='x')
ax.plot(samples_per_pixel, deon_vs_deon,
        label='Importance Sampling (relative)', marker='x')

plt.xscale('log')
plt.xlabel("Samples per pixel")
plt.ylabel("Variance compared to reference (ground truth)")
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
