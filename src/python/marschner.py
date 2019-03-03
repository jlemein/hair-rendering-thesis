# THis is mostly a copy of the C++ code in order to test functions and make visualizations

import numpy as np


def Gaussian(width, x):
    a = 1.0 / (width * np.sqrt(2.0 * np.pi))
    c = width
    # width of the curve is beta(might also be 0.5 * sigma)
    nom = np.square(x)
    den = 2.0 * np.square(width)
    return a * np.exp(-nom / den)


def HalfAngle(a, b):
    return 0.5 * (a + b)


class Marschner:
    def __init__(self):
        self.alphaR = np.radians(-7.5)
        self.alphaTT = -.5 * self.alphaR
        self.alphaTRT = -1.5*self.alphaR

        self.betaR = np.radians(7.5)
        self.betaTT = .5 * self.betaR
        self.betaTRT = 2.0 * self.betaR

    def M_r(self, theta_i, theta_r):
        theta_h = HalfAngle(theta_i, theta_r)
        return Gaussian(self.betaR, theta_h - self.alphaR)

    def M_tt(self, theta_i, theta_r):
        theta_h = HalfAngle(theta_i, theta_r)
        return Gaussian(self.betaTT, theta_h - self.alphaTT)

    def M_trt(self, theta_i, theta_r):
        theta_h = HalfAngle(theta_i, theta_r)
        return Gaussian(self.betaTRT, theta_h - self.alphaTRT)
