import numpy as np
import matplotlib.pyplot as plt
import marschner

# Fixing random state for reproducibility
np.random.seed(19680801)

mar = marschner.Marschner()

# Compute areas and colors
N = 100
thetai2 = np.linspace(0.0, 0.0, num=N)
thetai1 = np.linspace(-np.pi/4, -np.pi/4, num=N)
thetai3 = np.linspace(np.pi/4, np.pi/4, num=N)
thetar = np.linspace(-.5*np.pi, .5*np.pi, num=N)
# marschner_mr = mar.M_r(thetai, thetar)
# marschner_mtt = mar.M_tt(thetai, thetar)
# marschner_mtrt = mar.M_trt(thetai, thetar)


fig, axs = plt.subplots(1, 3, subplot_kw=dict(projection='polar'))
fig.suptitle("Marschner reflection")
axs[0].plot(thetar, mar.M_r(thetai1, thetar), label="Mr")
axs[0].plot(thetar, mar.M_r(thetai1, thetar), label="Mtt")
axs[0].plot(thetar, mar.M_r(thetai1, thetar), label="Mtrt")
axs[0].arrow(-np.pi/4, 5, 0, -5, alpha=0.5, width=0.015,
             edgecolor='black', facecolor='green', lw=2, zorder=5)
axs[0].set_xlim(thetar.min(), thetar.max())

axs[1].plot(thetar, mar.M_r(thetai2, thetar), label="Mr")
axs[1].plot(thetar, mar.M_tt(thetai2, thetar), label="Mtt")
axs[1].plot(thetar, mar.M_trt(thetai2, thetar), label="Mtrt")
axs[1].set_xlim(thetar.min(), thetar.max())
axs[1].arrow(0, 5, 0, -5, alpha=0.5, width=0.015,
             edgecolor='black', facecolor='green', lw=2, zorder=5)

axs[2].plot(thetar, mar.M_r(thetai3, thetar), label="Mr")
axs[2].plot(thetar, mar.M_tt(thetai3, thetar), label="Mtt")
axs[2].plot(thetar, mar.M_trt(thetai3, thetar), label="Mtrt")
axs[2].set_xlim(thetar.min(), thetar.max())
axs[2].arrow(np.pi/4, 5, 0, -5, alpha=0.5, width=0.015,
             edgecolor='black', facecolor='green', lw=2, zorder=5)

# c = ax.scatter(thetar, r, c=colors, s=area, cmap='hsv', alpha=0.75)
plt.show()
