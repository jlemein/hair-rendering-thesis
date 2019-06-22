from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
import sys

def readData(fileName):
    y = []
    with open(fileName) as f:
        lines = f.read().splitlines()
        for line in lines:
            l = line.split(" ")
            # thetas.append(float(l[0]))
            # phis.append(float(l[1]))
            y.append(float(l[0]))
    return y

def readAbAf(fileName):
    thetaD = []
    ab = []
    af = []
    with open(fileName) as f:
        lines = f.read().splitlines()
        for line in lines:
            l = line.split(" ")
            thetaD.append(float(l[0]))
            ab.append(float(l[1]))
            af.append(float(l[2]))
    return thetaD, ab, af


def readData2D(fileName):
    data = []
    with open(fileName) as f:
        dimensionLine = f.readline().split(" ")
        sizeX = int(dimensionLine[0])
        sizeY = int(dimensionLine[1])
        print("Dimensions: ", sizeX, sizeY)

        data = []
        for x in range(0, sizeX):
            row = []
            for y in range(0, sizeY):
                row.append(float(f.readline()))
            data.append(row)

    thetas = np.linspace(-.5*np.pi, .5*np.pi, sizeX)
    phis = np.linspace(-np.pi, np.pi, sizeY)
    X, Y = np.meshgrid(phis, thetas)
    return X, Y, np.array(data)

def getTheta(data, rowSize, theta):
    row = round(rowSize * (theta + .5*np.pi) / np.pi)
    print("Row: ", row)
    return data[row, :]

# y0 = readData("thetaI_0.data")
# y1 = readData("thetaI_1.data")
# y2 = readData("thetaI_2.data")
# y4 = readData("thetaI_4.data")

# phi0 = readData("phiI_0.data")
# phi1 = readData("phiI_1.data")
# phi2 = readData("phiI_2.data")
# phi4 = readData("phiI_4.data")




# thetas = np.linspace(-.5*np.pi, .5*np.pi, len(y1))
# phis = np.linspace(-np.pi, np.pi, len(phi0))



def plotFigures(fig, colIndex, columnCount, thetas, phis, data, title):
    axs0 = fig.add_subplot(3, columnCount, colIndex, projection='3d')
    axs0.set_title(title)
    axs0.plot_surface(x, y, data, cmap="viridis");
    axs0.set_ylabel(r'$\theta_i$')
    axs0.set_xlabel(r'$\phi_i$')


    axs1 = fig.add_subplot(3, columnCount, colIndex+columnCount)
    axs1.set_title(r'Response for different $\theta_i$, given $\phi_i$')
    axs1.plot(thetas, data[:, 17], color='red', label=r'$\phi_i=-\frac{2}{3}\pi$');
    axs1.plot(thetas, data[:, 34], color='blue', label=r'$\phi_i=-\frac{1}{3}\pi$');
    axs1.plot(thetas, data[:, 50], color='green', label=r'$\phi_i=0$');
    axs1.plot(thetas, data[:, 66], "--", color='blue', label=r'$\phi_i=-\frac{1}{3}\pi$');
    axs1.plot(thetas, data[:, 83], "--", color='red', label=r'$\phi_i=-\frac{2}{3}\pi$');
    axs1.set_ylabel(r'response ($y$)')
    axs1.set_xlabel(r'$\theta_i$')
    axs1.legend()

    axs3 = fig.add_subplot(3, columnCount, colIndex+2*columnCount)
    axs3.set_title(r'Response for different $\phi_i$, given $\theta_i$')
    axs3.plot(phis, data[17, :], color='red', label=r'$\theta_i=-\frac{2}{6}\pi$');
    axs3.plot(phis, data[34, :], color='blue', label=r'$\theta_i=-\frac{1}{6}\pi$');
    axs3.plot(phis, data[50, :], color='green', label=r'$\theta_i=0$');
    axs3.plot(phis, data[66, :], "--", color='blue', label=r'$\theta_i=\frac{1}{6}\pi$');
    axs3.plot(phis, data[83, :], "--", color='red', label=r'$\theta_i=\frac{2}{6}\pi$');
    axs3.set_ylabel(r'response($y$)')
    axs3.set_xlabel(r'$\phi_i$')
    axs3.legend()


#x, y, data = readData2D("test.data")
x, y, blonde0 = readData2D("blonde0.data")
x, y, blonde1 = readData2D("blonde1.data")
x, y, blonde2 = readData2D("blonde2.data")

x, y, brunette0 = readData2D("brunette0.data")
x, y, brunette1 = readData2D("brunette1.data")
x, y, brunette2 = readData2D("brunette2.data")

thetas = np.linspace(-.5*np.pi, .5*np.pi, 100)
phis = np.linspace(-np.pi, np.pi, 100)

fig = plt.figure(1)
fig.canvas.set_window_title('Blonde hair') 
plotFigures(fig, 1, 3, thetas, phis, blonde0, r"Scatter count ($n=0$)")
plotFigures(fig, 2, 3, thetas, phis, blonde1, r"Scatter count ($n=1$)")
plotFigures(fig, 3, 3, thetas, phis, blonde2, r"Scatter count ($n=2$)")
fig.savefig("blonde_hair.pdf", dpi=150)

fig = plt.figure(2)
fig.canvas.set_window_title('Brown hair') 
plotFigures(fig, 1, 3, thetas, phis, brunette0, r"Scatter count ($n=0$)")
plotFigures(fig, 2, 3, thetas, phis, brunette1, r"Scatter count ($n=1$)")
plotFigures(fig, 3, 3, thetas, phis, brunette2, r"Scatter count ($n=2$)")
fig.savefig("brown_hair.pdf", dpi=150)

fig = plt.figure(3)
fig.canvas.set_window_title('Absorption after scatter amount')
ax = fig.add_subplot(2, 1, 1)
ax.plot(thetas, brunette0[:, 50], "-", color="brown", label='n=0')
ax.plot(thetas, brunette1[:, 50], "--", color="brown", label='n=1')
ax.plot(thetas, brunette2[:, 50], ":", color="brown", label='n=2')
ax.plot(thetas, blonde0[:, 50], "-", color="orange", label='n=0')
ax.plot(thetas, blonde1[:, 50], '--', color="orange", label='n=1')
ax.plot(thetas, blonde2[:, 50], ':', color="orange", label='n=2')
ax.legend()
ax = fig.add_subplot(2, 1, 2)
ax.plot(phis, brunette0[50, :], "-", color="brown", label='n=0')
ax.plot(phis, brunette1[50, :], "--", color="brown", label='n=1')
ax.plot(phis, brunette2[50, :], ":", color="brown", label='n=2')
ax.plot(phis, blonde0[50, :], "-", color="orange", label='n=0')
ax.plot(phis, blonde1[50, :], '--', color="orange", label='n=1')
ax.plot(phis, blonde2[50, :], ':', color="orange", label='n=2')
ax.legend()
fig.savefig("absorption_after_scatter_amount.pdf", dpi=150)

#
# Displaying average forward and backward scattering attentuation
#
thetaD, blonde_ab, blonde_af = readAbAf("blonde_abaf.data")
thetaD, brunette_ab, brunette_af = readAbAf("brunette_abaf.data")
thetaD, white_ab, white_af = readAbAf("white_abaf.data")
fig = plt.figure(4)
fig.canvas.set_window_title('Average forward/backward scattering attenuation')
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'Transmittance')
ax.set_xlabel(r'Difference angle ($\theta_d$)')
ax.plot(thetaD, blonde_ab, "-", color="orange", label=r"$a_b$ Blonde hair")
ax.plot(thetaD, blonde_af, "-", color="yellow", label=r"$a_f$ Blonde hair")
ax.plot(thetaD, brunette_ab, "-", color="brown", label=r"$a_b$ Brown hair")
ax.plot(thetaD, brunette_af, "-", color="red", label=r"$a_f$ Brown hair")
ax.plot(thetaD, white_ab, "-", color="cyan", label=r"$a_b$ No absorption")
ax.plot(thetaD, white_af, "-", color="blue", label=r"$a_f$ No absorption")
ax.plot(thetaD, np.cos(np.array(thetaD)), "-", color="black", label=r"$\cos \theta_d$")
ax.legend()
fig.savefig("absorption_after_scatter_amount.pdf", dpi=150)



plt.show()