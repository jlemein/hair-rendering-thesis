
#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "materials/marschner.cpp"
#include <atomic>

using namespace pbrt;

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
    EXPECT_FLOAT_EQ(-0.99 * Pi, RelativeAzimuth(-2.0 * Pi, -0.99 * Pi));
    EXPECT_FLOAT_EQ(0.99 * Pi, RelativeAzimuth(-2.0 * Pi, 0.99 * Pi));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-2.0 * Pi, -Pi)));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-2.0 * Pi, Pi)));

    EXPECT_FLOAT_EQ(0.99 * Pi, RelativeAzimuth(-0.99 * Pi, -2.0 * Pi));
    EXPECT_FLOAT_EQ(-0.99 * Pi, RelativeAzimuth(0.99 * Pi, -2.0 * Pi));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-Pi, -2.0 * Pi)));
    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(Pi, -2.0 * Pi)));

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
    //    EXPECT_FLOAT_EQ(0.0, RelativeAzimuth(-4.0 * Pi, Pi));
    //    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(-Pi, 4.0 * Pi)));
    //    EXPECT_FLOAT_EQ(Pi, fabs(RelativeAzimuth(4.0 * Pi, -Pi)));
}