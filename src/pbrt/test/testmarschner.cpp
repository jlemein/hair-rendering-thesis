
#include "tests/gtest/gtest.h"
#include "pbrt.h"
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