
#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "sampler.h"
#include "materials/marschner.cpp"
#include <atomic>

#include <vector>

using namespace pbrt;

class DepressedCubic {
public:
    Float a, b, c;

    DepressedCubic(Float a, Float b, Float c) : a(a), b(b), c(c) {
    }

    Float operator()(Float x) {
        return a * x * x * x + b * x + c;
    }
};

// Set up sample phis for testing
static const int SAMPLE_SIZE = 100;

static Float EvaluateCubic(Float a, Float b, Float c, Float x) {
    return a * x * x * x + b * x + c;
}

/**
 This test assures that when we wrap gamma around its boundaries, that
 the root h = sin(gamma) still holds.
 **/
TEST(Marschner, RangeBoundGammaInversed) {
    for (Float gammaI = -Pi; gammaI < Pi; gammaI += 0.01) {
        EXPECT_NEAR(sin(gammaI), sin(RangeBoundGammaInversed(gammaI)), 1e-5);
    }
}

TEST(Marschner, EtaEccentricityIdentity) {
    EXPECT_FLOAT_EQ(1.55, EtaEccentricity(1.0, 1.55, 0.0));
    EXPECT_FLOAT_EQ(1.55, EtaEccentricity(1.0, 1.55, -.5 * Pi));
    EXPECT_FLOAT_EQ(1.55, EtaEccentricity(1.0, 1.55, .5 * Pi));
}

TEST(Marschner, EtaEccentricity) {
    EXPECT_GT(EtaEccentricity(0.9, 1.55, 0.0), 1.0);
    EXPECT_GT(EtaEccentricity(0.9, 1.55, -.5 * Pi), 1.0);
    EXPECT_GT(EtaEccentricity(0.9, 1.55, .5 * Pi), 1.0);
}


//
//TEST(Marschner, RangeBoundGamma) {
//    EXPECT_FLOAT_EQ(.5 * Pi - 0.1, RangeBoundGamma(.5 * Pi - 0.1));
//    EXPECT_FLOAT_EQ(.5 * Pi - 0.1, RangeBoundGamma(.5 * Pi + 0.1));
//    EXPECT_FLOAT_EQ(.5 * Pi - 0.1, RangeBoundGamma(-1.5 * Pi - 0.1));
//    EXPECT_FLOAT_EQ(.5 * Pi - 0.1, RangeBoundGamma(-1.5 * Pi + 0.1));
//
//    EXPECT_FLOAT_EQ(-.5 * Pi + 0.1, RangeBoundGamma(-.5 * Pi - 0.1));
//    EXPECT_FLOAT_EQ(-.5 * Pi + 0.1, RangeBoundGamma(-.5 * Pi + 0.1));
//    EXPECT_FLOAT_EQ(-.5 * Pi + 0.1, RangeBoundGamma(1.5 * Pi + 0.1));
//    EXPECT_FLOAT_EQ(-.5 * Pi + 0.1, RangeBoundGamma(1.5 * Pi - 0.1));
//
//    EXPECT_FLOAT_EQ(-.5 * Pi + 0.1, RangeBoundGamma(-.5 * Pi - 0.1));
//    EXPECT_FLOAT_EQ(-.5 * Pi + 0.1, RangeBoundGamma(-.5 * Pi + 0.1));
//
//
//    EXPECT_FLOAT_EQ(.0, RangeBoundGamma(Pi));
//    EXPECT_FLOAT_EQ(.0, RangeBoundGamma(2.0 * Pi));
//    EXPECT_FLOAT_EQ(.0, RangeBoundGamma(-2.0 * Pi));
//
//    // testing with near, because of precision errors
//    EXPECT_NEAR(.0, RangeBoundGamma(3.0 * Pi), 1e-3);
//    EXPECT_NEAR(.0, RangeBoundGamma(-3.0 * Pi), 1e-3);
//    EXPECT_NEAR(.0, RangeBoundGamma(4.0 * Pi), 1e-3);
//    EXPECT_NEAR(.0, RangeBoundGamma(-4.0 * Pi), 1e-3);
//}
//

TEST(Marschner, SolveRoot) {
    Float a = 16.0;
    Float b = 0.0;
    Float c = -6.0;
    Float d = -9.0;

    // ax^3 + bx^2 -c = d
    Float root;
    int nRoots = SolveDepressedCubic(a, c, d, &root);

    EXPECT_EQ(1, nRoots);
    EXPECT_NEAR(0.0, EvaluateCubic(a, c, d, root), 1e-5);
}
//
//TEST(Marschner, SolveRootsForTT) {
//    Float eta = 1.55;
//    const Float C = asin(1.0 / eta);
//
//    for (int i = 0; i < SAMPLE_SIZE; ++i) {
//        Float phi = -Pi + (i / (Float) SAMPLE_SIZE) * 2.0 * Pi;
//
//        // make sure that phi stays between [-Pi, Pi] due to floating point
//        // precision errors
//        phi = Clamp(phi, -Pi, Pi);
//
//        // depressed cubic equation: ax^3 + cx + d = 0
//        Float a = -8.0 * C / (Pi * Pi * Pi);
//        Float c = 6.0 * C / Pi - 2.0;
//        Float d = Pi - phi;
//
//        Float root;
//        int nRoots = SolveDepressedCubic(a, c, d, &root);
//
//        // always expect 1 root for TT scattering components
//        EXPECT_EQ(1, nRoots);
//
//        // result should always be 0, we check here for smaller than 1e-5
//        Float result = EvaluateCubic(a, c, d, root);
//        EXPECT_NEAR(fabs(result), 0.0, 1e-4);
//    }
//}
//
//TEST(Marschner, SolveRootsForTT2) {
//    Float eta = 1.55;
//
//    for (int i = 0; i < SAMPLE_SIZE; ++i) {
//        Float phi = -Pi + (i / (Float) SAMPLE_SIZE) * 2.0 * Pi;
//        Float gammaI = SolveGammaRoot_TT(phi, eta);
//
//        // always expect 1 root for TT scattering components
//        EXPECT_LE(gammaI, .5 * Pi);
//        EXPECT_GE(gammaI, -.5 * Pi);
//
//        // transforming the root should roughly be equal to the phi for which
//        // we found the root.
//
//        //EXPECT_NEAR(phi, ClampPhi(Phi(1, gammaI, GammaT(gammaI, eta))), 1e-5);
//    }
//}
//
//TEST(Marschner, SolveEquationWithThreeRoots) {
//
//    Float a = 1.0;
//    Float c = -15.0;
//    Float d = -4.0;
//    DepressedCubic fn(a, c, d);
//
//    // ax^3 + bx^2 -c = d
//    Float roots[3];
//    int nRoots = SolveDepressedCubic(a, c, d, roots);
//    EXPECT_EQ(3, nRoots);
//
//    EXPECT_NEAR(fabs(fn(roots[0])), .0, 1e-4);
//    EXPECT_NEAR(fabs(fn(roots[1])), .0, 1e-4);
//    EXPECT_NEAR(fabs(fn(roots[2])), .0, 1e-4);
//}
//

//TEST(Marschner, SolveRootsForTRT) {
//    //Float eta = 1.55;
//    //const Float C = asin(1.0 / eta);
//    //const Float a = -8.0 * 2.0 * C / (Pi * Pi * Pi);
//    //const Float c = 6.0 * 2.0 * C / Pi - 2.0;
//
//    const Float a = -0.352144;
//    const Float c = 0.606644;
//    // const Float d = 0.024114;
//
//    // walk around cylinder for incoming phi values
//    for (Float phi = -Pi; phi <= Pi; phi += 0.01) {
//
//        // make sure that phi stays between [-Pi, Pi] due to floating point
//        // precision errors
//        phi = Clamp(phi, -Pi, Pi);
//
//        // depressed cubic equation: ax^3 + cx + d = 0
//        Float d = 2.0 * Pi - phi;
//
//        Float roots[3];
//        int nRoots = SolveDepressedCubic(a, c, d, roots);
//
//        printf("roots: %d\n", nRoots);
//        if (nRoots == 3) {
//            printf("phi = %f\n", phi);
//            printf("roots: %f -- %f -- %f\n\n", roots[0], roots[1], roots[2]);
//        }
//
//        // TRT could result in one or three roots
//        EXPECT_TRUE(nRoots == 1 || nRoots == 3);
//
//        for (int i = 0; i < nRoots; ++i) {
//            // expect evaluation of function to be zero
//
//            EXPECT_LT(fabs(EvaluateCubic(a, c, d, roots[i])), 1e-3);
//        }
//    }
//}

TEST(Marschner, SolveThreeRootsForTRT) {

    const Float a = -0.352144;
    const Float c = 0.606644;
    const Float d = 0.024114;

    Float roots[3];
    int nRoots = SolveDepressedCubic(a, c, d, roots);

    printf("roots: %d\n", nRoots);
    if (nRoots == 3) {
        //printf("phi = %f\n", phi);
        printf("roots: %f -- %f -- %f\n\n", roots[0], roots[1], roots[2]);
    }

    // TRT could result in one or three roots
    EXPECT_TRUE(nRoots == 3);

    for (int i = 0; i < nRoots; ++i) {
        // expect evaluation of function to be zero

        EXPECT_LT(fabs(EvaluateCubic(a, c, d, roots[i])), 1e-3);
    }
}

TEST(Marschner, SolveThreeOutOfBoundsRootsForTRT) {

    const Float etaPerp = 1.550605;
    const Float phi = 0.346892;
    const Float a = -0.361684;
    const Float c = 0.677259;
    const Float d = -phi;

    Float roots[3];
    int nRoots = SolveDepressedCubic(a, c, d, roots);

    printf("roots: %d\n", nRoots);
    if (nRoots == 3) {
        //printf("phi = %f\n", phi);
        printf("roots: %f -- %f -- %f\n\n", roots[0], roots[1], roots[2]);
    }

    // TRT could result in one or three roots
    //EXPECT_TRUE(nRoots == 3);

    for (int i = 0; i < nRoots; ++i) {
        EXPECT_LT(fabs(EvaluateCubic(a, c, d, roots[i])), 1e-3);
        EXPECT_NEAR(PhiApprox(2, roots[i], etaPerp), phi, 1e-5);
    }
}
//
//TEST(Marschner, SolveRootsForTRT_LowEta) {
//    Float eta = 2.5;
//    const Float C = asin(1.0 / eta);
//    const Float a = -8.0 * 2.0 * C / (Pi * Pi * Pi);
//    const Float c = 6.0 * 2.0 * C / Pi - 2.0;
//
//    // walk around cylinder for incoming phi values
//    for (Float phi = -Pi; phi <= Pi; phi += 0.2 * Pi) {
//
//        // make sure that phi stays between [-Pi, Pi] due to floating point
//        // precision errors
//        phi = Clamp(phi, -Pi, Pi);
//
//        // depressed cubic equation: ax^3 + cx + d = 0
//        Float d = 2.0 * Pi - phi;
//
//        Float roots[3];
//        int nRoots = SolveDepressedCubic(a, c, d, roots);
//
//        // TRT could result in one or three roots
//        EXPECT_TRUE(nRoots == 1 || nRoots == 3);
//
//        for (int i = 0; i < nRoots; ++i) {
//            // expect evaluation of function to be zero
//
//            EXPECT_LT(fabs(EvaluateCubic(a, c, d, roots[i])), 1e-4);
//        }
//    }
//}
//
//TEST(Marschner, Fresnel) {
//    Float etaPerp = 1.55;
//    Float etaPar = 1.55;
//
//    for (Float gammaI = -0.5 * Pi; gammaI <= 0.5 * Pi; gammaI += 0.01 * Pi) {
//
//        Float frDielectric = FrDielectric(cos(gammaI), 1.0, 1.55);
//        Float fresnel = Fresnel(etaPerp, etaPar, gammaI);
//
//        EXPECT_FLOAT_EQ(frDielectric, fresnel);
//    }
//}
//
//// TODO:
//// * test if dphidh(0, gamma, etaPerp) equals dphidh_r
//// * test if dphi comes close to gradient and is always above 0
//// * check if gammaRoots equals gammaRoots_TT for p = 1
//// * check solutions for caustics and see if in between is the switch for 1 solution vs 3 solutions
//
//TEST(Marschner, RelativeAzimuth) {
//
//    // test all variants where difference is 0 degrees
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, -2.0 * Pi));
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, 2.0 * Pi));
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, 0.0));
//
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(2.0 * Pi, -2.0 * Pi));
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(2.0 * Pi, 2.0 * Pi));
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(2.0 * Pi, 0.0));
//
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(0.0, 0.0));
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(0.0, 2.0 * Pi));
//    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(0.0, -2.0 * Pi));
//
//    // test difference of 180 degrees (-pi)
//    EXPECT_FLOAT_EQ(-0.99 * Pi, RelativeAzimuth(0.0, -0.99 * Pi));
//    EXPECT_FLOAT_EQ(0.99 * Pi, RelativeAzimuth(0.0, 0.99 * Pi));
//    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(0.0, -Pi)));
//    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(0.0, Pi)));
//
//    EXPECT_FLOAT_EQ(0.99 * Pi, RelativeAzimuth(-0.99 * Pi, 0.0));
//    EXPECT_FLOAT_EQ(-0.99 * Pi, RelativeAzimuth(0.99 * Pi, 0.0));
//    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-Pi, 0.0)));
//    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(Pi, 0.0)));
//
//    // test difference of 90 degrees (pi/2)
//    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(.0, .5 * Pi));
//    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(2 * Pi, .5 * Pi));
//    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-2 * Pi, .5 * Pi));
//
//    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.0, -.5 * Pi));
//    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(2 * Pi, -.5 * Pi));
//    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(-2 * Pi, -.5 * Pi));
//
//    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.5 * Pi, .0));
//    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.5 * Pi, 2 * Pi));
//    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.5 * Pi, -2 * Pi));
//
//    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-.5 * Pi, 0.0));
//    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-.5 * Pi, 2 * Pi));
//    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-.5 * Pi, -2 * Pi));
//
//    // test wrapping around (tests are failing because Pi is not precise enough
//    // but this is no problem for us, since wrapping does not occur)
//    //    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, 4.0 * Pi));
//    //    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(4.0 * Pi, -2.0 * Pi));
//    //    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-Pi, 4.0 * Pi)));
//    //    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(4.0 * Pi, -Pi)));
//}

double _integrateMonteCarloFront(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    double V = 2.0 * Pi;
    double sum = .0;

    for (int i = 0; i < nSamples; ++i) {
        //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
        Vector3f wr = SampleFrontHemisphere(Point2f(distribution(generator), distribution(generator)));
        sum += f(wr, wi);
    }

    return V / static_cast<double> (nSamples) * sum;
}

double _integrateMonteCarloBack(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    double V = 2.0 * Pi;
    double sum = .0;

    for (int i = 0; i < nSamples; ++i) {
        //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
        Vector3f wr = SampleBackHemisphere(Point2f(distribution(generator), distribution(generator)));
        sum += f(wr, wi);
    }

    return V / static_cast<double> (nSamples) * sum;
}

Vector3f UniformSampleSphere(const Point2f &u) {
    Float z = 1 - 2 * u[0];
    Float r = std::sqrt(std::max((Float) 0, (Float) 1 - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

double _integrateMonteCarlo(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    double V = 4.0 * Pi;
    double sum = .0;

    for (int i = 0; i < nSamples; ++i) {
        //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
        Vector3f wr = UniformSampleSphere(Point2f(distribution(generator), distribution(generator)));
        sum += f(wr, wi);
    }

    return V / static_cast<double> (nSamples) * sum;
}

double _integrateMonteCarlo(std::function<double(const Float) > f, int nSamples) {

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(.0, 1.0);

    double V = Pi;
    double sum = .0;

    for (int i = 0; i < nSamples; ++i) {
        //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
        Float x = -.5 * Pi + distribution(generator) * Pi;
        sum += f(x);
    }

    return V / static_cast<double> (nSamples) * sum;
}

//alphaR: 0.0523599 alphaTT: -0.0261799 alphaTRT: -0.0785398 betaR: 0.244346 betaTT: 0.139626 betaTRT: 0.383972
//causticFadeRange: 0.3 causticIntensityLimit: 0.5 causticWidth: 0.174533
//eccentricity: 0.9 eta: 1.55 glintScaleFactor: 0.4 hairRadius: 1
//sigmaA: 0.432, 0.612, 0.98
const Float eta = 1.55;
const Float eccentricity = 0.9;
const Vector3f alpha = Vector3f(0.0523599, -0.0261799, -0.0785398);
const Vector3f alphaSquared = Vector3f(alpha.x*alpha.x, alpha.y*alpha.y, alpha.z*alpha.z);
const Vector3f alphaSqrt = Vector3f(sqrt(alpha.x), sqrt(alpha.y), sqrt(alpha.z));

const Vector3f beta = Vector3f(0.244346, 0.139626, 0.383972);
const Vector3f betaSquared = Vector3f(beta.x*beta.x, beta.y*beta.y, beta.z*beta.z);
const Vector3f betaSqrt = Vector3f(sqrt(beta.x), sqrt(beta.y), sqrt(beta.z));

const Float glintScale = 0.4;
const Float causticWidth = 0.174533;
const Float causticFade = 0.3;
const Float causticIntensityLimit = 0.5;
const Float hairRadius = 1.0;
const Float sigmaARgb[3] = {.0, .0, .0};
const Spectrum sigmaA = Spectrum::FromRGB(sigmaARgb);
const SurfaceInteraction si = SurfaceInteraction();

MarschnerBSDF* marschner = new MarschnerBSDF(si, alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2], hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);
MarschnerBSDF* marschnerSquared = new MarschnerBSDF(si, alphaSquared[0], alphaSquared[1], alphaSquared[2], betaSquared[0], betaSquared[1], betaSquared[2], hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);
MarschnerBSDF* marschnerSqrt = new MarschnerBSDF(si, alphaSqrt[0], alphaSqrt[1], alphaSqrt[2], betaSqrt[0], betaSqrt[1], betaSqrt[2], hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);

TEST(Marschner, BsdfShouldBeEnergyConservant) {
    MyRandomSampler sampler(0.0, 1.0);

    auto fn = [&](const Vector3f wr, const Vector3f wi) {
        Float cosAngle = Dot(wr, wi);
        return marschner->f(wr, wi).y(); // * fabs(cosAngle);
    };

    auto fnSquared = [&](const Vector3f wr, const Vector3f wi) {
        Float cosAngle = Dot(wr, wi);
        return marschnerSquared->f(wr, wi).y() * fabs(cosAngle);
    };

    auto fnSqrt = [&](const Vector3f wr, const Vector3f wi) {
        Float cosAngle = Dot(wr, wi);
        return marschnerSqrt->f(wr, wi).y() * fabs(cosAngle);
    };

    auto fnSphereSurface = [&](const Vector3f wr, const Vector3f wi) {
        return 1.0;
    };

    Vector3f wii = SampleBackHemisphere(0.0, sampler.next());

    Float sphere = _integrateMonteCarlo(wii, fn, 10000);

    printf("integrated BSDF energy: %f\n", sphere);

    EXPECT_LE(sphere, 1.2);
    EXPECT_GE(sphere, .8);
}

TEST(Marschner, NormalizedGaussianSumsToOne) {
    MyRandomSampler sampler(0.0, 1.0);
    auto fn = [&](const Float thetaH) {
        return Gaussian(beta.x, thetaH);
    };

    Float area = _integrateMonteCarlo(fn, 10000);
    printf("area under normalized gaussian graph: %f\n", area);
    EXPECT_NEAR(1.0, area, 0.03);
}
