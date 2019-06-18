import numpy as np
import matplotlib.pyplot as plt

def balanceHeuristic(x, currentMIS, listMIS):
    currentSamples, currentPdf = currentMIS
    ns = len(currentSamples)
    ps = currentPdf(x) if callable(currentPdf) else currentPdf
    
    denom = .0;
    for (samples, pdf) in listMIS:
        ni = len(samples)
        pi = pdf(x) if callable(pdf) else pdf
        denom += ni*pi 
        
    return ns*ps / denom

def integrate(fn, samples, pdf):
    sum = .0
    
    for x in samples:
        p = pdf(x) if callable(pdf) else pdf
        sum += fn(x) / p

    return sum / len(samples)

def integrateMIS(fn, listMIS):
    result = .0
    
    for (samples, pdf) in listMIS:
        sumComponent = .0
        for x in samples:
            p = pdf(x) if callable(pdf) else pdf 
            w = balanceHeuristic(x, (samples, pdf), listMIS)
            sumComponent += fn(x) * w / p

        print("Sum component: ", sumComponent/len(samples))
        result += sumComponent / len(samples)

    return result

def integrateMIS2(fn, samples1, pdf1, samples2, pdf2):
    sum = .0;

    for x in samples1:
        sum += fn(x) * len(samples1) / (len(samples1) * pdf1(x) + len(samples2) * pdf2(x))
    
    for x in samples2:
        sum += fn(x) * len(samples2) / (len(samples1) * pdf1(x) + len(samples2) * pdf2(x))

    return sum


def fancyFn(x):
    if x < 0:
        return 0

    if x >= 0 and x <= np.pi:
        return np.sin(x)

    a = 0.15 * (x - np.pi)
    return a*a;

def f2(x):
    a = (x - np.pi)
    return a * a

# pdf
def p(x):
    return np.power(x-np.pi, 2) / 117.87

# CDF of p(x) = (x - pi)^2
def P(x):
    return (1.0/117.87) * (np.power(x - np.pi, 3)/3.0 + 10.3354)

def invP(x):
    a = 1.0/117.87
    b = 10.3354
    return np.pi + np.cbrt(3*x/a - 3*b)

def p2(x):
    return (2.0*np.pi - x) / (2.0*np.pi*np.pi)

def P2(x):
    return (2.0*np.pi*x - .5*x*x) / (2.0*np.pi*np.pi)

def invP2(x):
    return 2.0*np.pi - np.sqrt(4.0*np.pi*np.pi - x*4.0*np.pi*np.pi)

def pdf(x):
    return np.power(x-np.pi, 3) / 3.0

vFancyFn = np.vectorize(fancyFn)
#vFancyFn = np.vectorize(p)


def sampleLinear(a, b, nSamples):
    return np.linspace(a, b, nSamples)

def sampleUniform(a, b, nSamples):
    return a + (b-a) * np.random.uniform(0.0, 1.0, nSamples)

def sampleFromInversion(inversedFn, a, b, nSamples):
    return inversedFn(np.random.uniform(0.0, 1.0, nSamples), a, b)

def importanceSamples(invertedCdf, nSamples):
    return invertedCdf(np.random.uniform(0.0, 1.0, nSamples))


# x: canonical random variable between 0 and 1
# a: lower bound
# b: upper bound
def invCdfUniform(x, a, b):
    return (b-a) * x + a


samplesLinear = sampleLinear(0, 10, 100000)
samplesUniform = sampleUniform(0, 10, 10000)
samplesFromInversion = sampleFromInversion(invCdfUniform, 0, 10, 1000)
samplesImportance = importanceSamples(invP, 100000)

samplesImportance1 = sampleUniform(0, 5, 1000)#importanceSamples(invP, 10)
samplesImportance2 = sampleUniform(5, 10, 1000)#importanceSamples(invP2, 20)


# print("Samples for 2 MIS: ", samplesImportance1, samplesImportance2)

print("Integrated linearly gives: ", integrate(vFancyFn, samplesLinear, lambda x: 0.1 ))
# print("Integrated uniformly gives: ", integrate(vFancyFn, samplesUniform, 1/10.0))
# print("Integrated from inversion gives: ", integrate(vFancyFn, samplesFromInversion, 0.1))
# print("Integrated with importance sampling: ", integrate(vFancyFn, samplesImportance, p))
print("Integrated with MIS: ", integrateMIS(vFancyFn, [(samplesImportance1, 0.2), (samplesImportance2, 0.2)]))
print("Integrated with MIS2: ", integrateMIS2(vFancyFn, 
    samplesImportance1, lambda x: 0.2 if x >= 0 and x <= 5 else 0, 
    samplesImportance2, lambda x: 0.2 if x >= 5 and x <= 10 else 0))

# Fixing random state for reproducibility
np.random.seed(19680801)

# Compute areas and colors
N = 100
x = np.linspace(0.0, 10.0, num=N)
y = vFancyFn(x)



fig, axs = plt.subplots(1, 2)
fig.suptitle("Marschner reflection")
axs[0].plot(x, y, label="y")
axs[0].plot(x, p(x), label="pdf")
axs[0].plot(x, P(x), label="CDF")
axs[0].set_xlim(x.min(), x.max())

xx = np.linspace(0.0, 1.0, 100)
x2 = np.linspace(0.0, 1.0, 100)
f = np.vectorize(invP2)
axs[1].plot(xx, f(xx), label="invCDF")

plt.show()
