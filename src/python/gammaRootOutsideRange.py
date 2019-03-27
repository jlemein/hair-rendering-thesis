import numpy as np
import matplotlib.pyplot as plt
import marschner

a = -0.360782
c = 0.670584
etaPerp = 1.553825
phi = -1.609246
d = 2.0 * np.pi - phi


def evaluateCubic(a, b, c, d, x):
    return a*x*x*x + b*x*x + c*x + d

def evaluateApproxPhi(p, etaPerp, gamma, phi):
    c = np.arcsin(1.0/etaPerp)
    d = (p%2)*np.pi - phi;
    while (d > np.pi): d -= 2.0 * np.pi;
    while (d < -np.pi): d += 2.0 * np.pi;
    return gamma * (6.0 * p * c / np.pi - 2.0) - gamma * gamma * gamma * 8.0*p*c/(np.pi*np.pi*np.pi) + d

def evaluateApprox(p, etaPerp, gamma):
    c = np.arcsin(1.0/etaPerp)
    return gamma * (6.0 * p * c / np.pi - 2.0) - gamma * gamma * gamma * 8.0*p*c/(np.pi*np.pi*np.pi) 


def evaluateReal(p, etaPerp, gamma):
    gammat = np.arcsin(np.sin(gamma)/etaPerp)
    return 2.0*p*gammat - 2.0*gamma + (p % 2)*np.pi

gammas = np.linspace(-np.pi, np.pi, 100)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
fig.suptitle(r'Incident angle $\gamma_i$ (in radians)', fontsize=16)

ax1.axhline(y=0)
ax1.axvline(x=-.5*np.pi, linestyle='dashed')
ax1.axvline(x=.5*np.pi, linestyle='dashed')
ax1.set_title("Scattering for R (p = 0)")

ax2.axhline(y=0)
ax2.axvline(x=-.5*np.pi, linestyle='dashed')
ax2.axvline(x=.5*np.pi, linestyle='dashed')
ax2.set_title("Scattering for TT (p = 1)")

ax3.axhline(y=0)
ax3.axvline(x=-.5*np.pi, linestyle='dashed')
ax3.axvline(x=.5*np.pi, linestyle='dashed')
ax3.set_title("Scattering for TRT (p = 2)")
# ax1.set_ylabel(r'Exitant angle $\phi_i$ (in radians)')

# R
curveRMin = evaluateApprox(0, etaPerp, gammas) + np.pi
curveRMax = evaluateApprox(0, etaPerp, gammas) - np.pi
curveR0 = evaluateApprox(0, etaPerp, gammas)

curveTTMin = evaluateApproxPhi(1, etaPerp, gammas, np.pi)
curveTTMax = evaluateApproxPhi(1, etaPerp, gammas, - np.pi)
curveTT0 = evaluateApproxPhi(1, etaPerp, gammas, 0)

curveTRTMin = evaluateApproxPhi(2, etaPerp, gammas, np.pi)
curveTRTMax = evaluateApproxPhi(2, etaPerp, gammas, -np.pi)
curveTRT0 = evaluateApproxPhi(2, etaPerp, gammas, 0)

ax1.fill_between(gammas, curveRMin, curveRMax, facecolor='salmon')
ax2.fill_between(gammas, curveTTMin, curveTTMax, facecolor='ivory')
ax3.fill_between(gammas, curveTRTMin, curveTRTMax, facecolor='lightblue')

# R
ax1.plot(gammas, curveR0, color='red',
         label=r'$R (\phi = 0)$', linestyle='dashed')
ax1.plot(gammas, curveRMin, color='red')
ax1.plot(gammas, curveRMax, color='red')


# TT
ax2.plot(gammas, curveTT0, color='orange',
         label=r'$TT (\phi = 0)$', linestyle='dashed')
ax2.plot(gammas, curveTTMin, color='orange')
ax2.plot(gammas, curveTTMax, color='orange')

# TRT
ax3.plot(gammas, curveTRT0, color='blue',
         label=r'$TRT (\phi = 0)$', linestyle='dashed')
ax3.plot(gammas, curveTRTMin, color='blue')
ax3.plot(gammas, curveTRTMax, color='blue')


plt.show()
