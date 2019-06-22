
from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
import sys

def sample(u0, u1, beta):
    print("-2.0 * log(", u0, ") = ", -2*np.log(u0))
    return beta * np.sqrt(-2*np.log(u0)) * np.cos(2.0 * np.pi * u1)

x = np.linspace(0.0, 1.0)
y = np.linspace(0.0, 1.0)

data = []
for u1 in range(0, 50):
    row = []
    for u2 in range(0, 50):
        row.append(sample(u1/50.0, u2/50.0, 1.0));
    data.append(row);


print("TEST: ", sample(0.99, 0.01, 1.0))
print(data);
X, Y = np.meshgrid(x, y)

fig = plt.figure()
plot = fig.add_subplot(2, 1, 1, projection='3d')
#plot.set_title(title)
plot.plot_surface(X, Y, np.array(data), cmap="viridis");
plot.set_ylabel(r'$\theta_i$')
plot.set_xlabel(r'$\phi_i$')

plt.show()


    # axs1.set_title(r'Response for different $\theta_i$, given $\phi_i$')
    # axs1.plot(thetas, data[:, 17], color='red', label=r'$\phi_i=-\frac{2}{3}\pi$');
    # axs1.plot(thetas, data[:, 34], color='blue', label=r'$\phi_i=-\frac{1}{3}\pi$');
    # axs1.plot(thetas, data[:, 50], color='green', label=r'$\phi_i=0$');
    # axs1.plot(thetas, data[:, 66], "--", color='blue', label=r'$\phi_i=-\frac{1}{3}\pi$');
    # axs1.plot(thetas, data[:, 83], "--", color='red', label=r'$\phi_i=-\frac{2}{3}\pi$');
    # axs1.set_ylabel(r'response ($y$)')
    # axs1.set_xlabel(r'$\theta_i$')
    # axs1.legend()