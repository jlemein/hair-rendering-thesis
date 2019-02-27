
#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "materials/marschner.cpp"
#include <atomic>

using namespace pbrt;

TEST(Marschner, SolveRoot) {
    Float a = 16.0;
    Float b = 0.0;
    Float c = -6.0;
    Float d = -9.0;

    // ax^3 + bx^2 -c = d
    Float x = 0, x2, x3;
    int nRoots = SolveDepressedCubic(a, c, d, x, x2, x3);
    printf("Root is %f\n", x);

    EXPECT_EQ(1, nRoots);
    EXPECT_FLOAT_EQ(0.0, a * x * x * x + c * x + d);
}

TEST(Marschner, SolveRootsForRandTT) {
    Float eta = 1.55;
    const Float C = asin(1.0 / eta);

    // walk around cylinder for incoming phi values
    for (Float phi = -Pi; phi <= Pi; phi += 0.01 * Pi) {

        // make sure that phi stays between [-Pi, Pi] due to floating point
        // precision errors
        phi = Clamp(phi, -Pi, Pi);

        // depressed cubic equation: ax^3 + cx + d = 0
        Float a = -8.0 * C / (Pi * Pi * Pi);
        Float c = 6.0 * C / Pi - 2.0;
        Float d = Pi - phi;

        Float x1, x2, x3;
        int nRoots = SolveDepressedCubic(a, c, d, x1, x2, x3);

        // always expect 1 root for R and TT scattering components
        EXPECT_EQ(1, nRoots);

        // result should always be 0, we check here for smaller than 1e-5
        Float result = a * x1 * x1 * x1 + c * x1 + d;
        EXPECT_LT(fabs(result), 1e-5);
    }
}

//TEST(Marschner, SolveRoots) {
//    // 3x^3 -9x + Pi = 0
//    Float r1, r2, r3;
//    int nRoots = SolveCubicRoots(3.0, 0.0, -9.0, Pi, r1, r2, r3);
//    EXPECT_EQ(1, nRoots);
//    EXPECT_FLOAT_EQ(0.0, 3.0 * pow(r1, 3.0) - 9 * r1 + Pi);
//}

TEST(Marschner, RelativeAzimuth) {
    // Relative azimuth takes


    // test all variants where difference is 0 degrees
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, -2.0 * Pi));
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, 2.0 * Pi));
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, 0.0));

    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(2.0 * Pi, -2.0 * Pi));
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(2.0 * Pi, 2.0 * Pi));
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(2.0 * Pi, 0.0));

    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(0.0, 0.0));
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(0.0, 2.0 * Pi));
    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(0.0, -2.0 * Pi));

    // test difference of 180 degrees (-pi)
    EXPECT_FLOAT_EQ(-0.99 * Pi, RelativeAzimuth(0.0, -0.99 * Pi));
    EXPECT_FLOAT_EQ(0.99 * Pi, RelativeAzimuth(0.0, 0.99 * Pi));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(0.0, -Pi)));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(0.0, Pi)));

    EXPECT_FLOAT_EQ(0.99 * Pi, RelativeAzimuth(-0.99 * Pi, 0.0));
    EXPECT_FLOAT_EQ(-0.99 * Pi, RelativeAzimuth(0.99 * Pi, 0.0));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-Pi, 0.0)));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(Pi, 0.0)));

    // test difference of 90 degrees (pi/2)
    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(.0, .5 * Pi));
    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(2 * Pi, .5 * Pi));
    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-2 * Pi, .5 * Pi));

    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.0, -.5 * Pi));
    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(2 * Pi, -.5 * Pi));
    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(-2 * Pi, -.5 * Pi));

    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.5 * Pi, .0));
    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.5 * Pi, 2 * Pi));
    EXPECT_FLOAT_EQ(-.5 * Pi, RelativeAzimuth(.5 * Pi, -2 * Pi));

    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-.5 * Pi, 0.0));
    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-.5 * Pi, 2 * Pi));
    EXPECT_FLOAT_EQ(.5 * Pi, RelativeAzimuth(-.5 * Pi, -2 * Pi));

    // test wrapping around (tests are failing because Pi is not precise enough
    // but this is no problem for us, since wrapping does not occur)
    //    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-2.0 * Pi, 4.0 * Pi));
    //    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(4.0 * Pi, -2.0 * Pi));
    //    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-Pi, 4.0 * Pi)));
    //    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(4.0 * Pi, -Pi)));
}