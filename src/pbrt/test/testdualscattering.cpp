
#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "materials/dualscattering.cpp"
#include <atomic>

#include <vector>

using namespace pbrt;

TEST(DualScattering, SampleFrontHemisphere) {
    std::vector<Vector3f> frontVectors;

    frontVectors.push_back(SampleFrontHemisphere(-.5 * Pi, 0.0));
    frontVectors.push_back(SampleFrontHemisphere(-.5 * Pi, 0.5));
    frontVectors.push_back(SampleFrontHemisphere(-.5 * Pi, 1.0));
    frontVectors.push_back(SampleFrontHemisphere(.0, 0.0));
    frontVectors.push_back(SampleFrontHemisphere(.0, 0.5));
    frontVectors.push_back(SampleFrontHemisphere(.0, 1.0));
    frontVectors.push_back(SampleFrontHemisphere(.5 * Pi, 0.0));
    frontVectors.push_back(SampleFrontHemisphere(.5 * Pi, 0.5));
    frontVectors.push_back(SampleFrontHemisphere(.5 * Pi, 1.0));
    frontVectors.push_back(SampleFrontHemisphere(Point2f(0.0, 0.0)));
    frontVectors.push_back(SampleFrontHemisphere(Point2f(0.5, 0.5)));
    frontVectors.push_back(SampleFrontHemisphere(Point2f(0.1, 0.9)));
    frontVectors.push_back(SampleFrontHemisphere(Point2f(0.9, 0.1)));
    frontVectors.push_back(SampleFrontHemisphere(Point2f(1.0, 1.0)));

    for (const auto& w : frontVectors) {
        EXPECT_FLOAT_EQ(w.Length(), 1.0);
        EXPECT_GE(w.x, 0.0);
    }
}

TEST(DualScattering, SampleBackHemisphere) {
    std::vector<Vector3f> backVectors;

    backVectors.push_back(SampleBackHemisphere(-.5 * Pi, 0.0));
    backVectors.push_back(SampleBackHemisphere(-.5 * Pi, 0.5));
    backVectors.push_back(SampleBackHemisphere(-.5 * Pi, 1.0));
    backVectors.push_back(SampleBackHemisphere(.0, 0.0));
    backVectors.push_back(SampleBackHemisphere(.0, 0.5));
    backVectors.push_back(SampleBackHemisphere(.0, 1.0));
    backVectors.push_back(SampleBackHemisphere(.5 * Pi, 0.0));
    backVectors.push_back(SampleBackHemisphere(.5 * Pi, 0.5));
    backVectors.push_back(SampleBackHemisphere(.5 * Pi, 1.0));
    backVectors.push_back(SampleBackHemisphere(Point2f(0.0, 0.0)));
    backVectors.push_back(SampleBackHemisphere(Point2f(0.5, 0.5)));
    backVectors.push_back(SampleBackHemisphere(Point2f(0.1, 0.9)));
    backVectors.push_back(SampleBackHemisphere(Point2f(0.9, 0.1)));
    backVectors.push_back(SampleBackHemisphere(Point2f(1.0, 1.0)));

    for (const auto& w : backVectors) {
        EXPECT_FLOAT_EQ(w.Length(), 1.0);
        EXPECT_LE(w.x, 0.0);
    }
}
