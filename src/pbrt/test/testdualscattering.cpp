
#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "materials/dualscattering.cpp"
#include "materials/dualscatteringlookup.h"
#include <atomic>

#include <vector>

using namespace pbrt;

const Float eta = 1.55;
const Float eccentricity = 1.0;
const Vector3f alpha = Vector3f(1.0, 1.0, 1.0);
const Vector3f beta = Vector3f(1.0, 1.0, 1.0);
const Float glintScale = 0.5;
const Float causticWidth = 1.5;
const Float causticFade = 1.0;
const Float causticIntensityLimit = 1.0;
const Float hairRadius = 1.0;
const Float sigmaARgb[3] = {0.4, 0.5, 0.1};
const Spectrum sigmaA = Spectrum::FromRGB(sigmaARgb);
const SurfaceInteraction si = SurfaceInteraction();

//TEST(DualScattering, AverageBackwardScatteringAttenuation) {
//
//    MarschnerBSDF* marschner = new MarschnerBSDF(si, alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2], hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);
//    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, eta, marschner, alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2], "unnamed.vdb");
//    const DualScatteringLookup* lookup = DualScatteringLookup::Get(dualScattering);
//
//    const int SAMPLES = 10;
//    for (int i = 0; i < SAMPLES; ++i) {
//        Float thetaD = -.5 * Pi + (static_cast<Float> (i) / SAMPLES) * Pi;
//        Spectrum response = dualScattering->AverageBackwardScatteringAlpha(thetaD);
//        Float rgb[3];
//        response.ToRGB(rgb);
//        printf("thetaD: %f -- rgb: %f %f %f\n", thetaD, rgb[0], rgb[1], rgb[2]);
//
//        EXPECT_GE(rgb[0], 0.0);
//        EXPECT_GE(rgb[1], 0.0);
//        EXPECT_GE(rgb[2], 0.0);
//    }
//}

//TEST(DualScattering, AverageForwardScatteringAttenuation) {
//
//    Spectrum theta0 = dualScattering->AverageForwardScatteringAttenuation(0.0);
//}

TEST(DualScattering, NumberGeneration) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    Float avg = 0.0,
            max = std::numeric_limits<Float>::min(),
            min = std::numeric_limits<Float>::max();
    for (int i = 0; i < 10000; ++i) {
        Float r = dis(gen);
        avg += r;
        min = std::min(r, min);
        max = std::max(r, max);
        EXPECT_LE(r, 1.0);
        EXPECT_GE(r, 0.0);
    }

    avg /= 10000.0;
    printf("Min: %f, Max: %f, Average: %f\n", min, max, avg);
    EXPECT_NEAR(min, .0, 0.01);
    EXPECT_NEAR(max, 1.0, 0.01);
    EXPECT_NEAR(avg, .5, 0.02);
}

//TEST(DualScattering, RandomSamplerNumberGeneration) {
//    const int N_SAMPLES = 10000;
//    RandomSampler sampler(N_SAMPLES);
//
//    Float avg = 0.0,
//            max = std::numeric_limits<Float>::min(),
//            min = std::numeric_limits<Float>::max();
//    for (auto i = 0; i < N_SAMPLES; ++i) {
//        Float r = sampler.Get1D()
//                avg += r;
//        min = std::min(r, min);
//        max = std::max(r, max);
//        EXPECT_LE(r, 1.0);
//        EXPECT_GE(r, 0.0);
//    }
//
//    avg /= 10000.0;
//    printf("Min: %f, Max: %f, Average: %f\n", min, max, avg);
//    EXPECT_NEAR(min, .0, 0.01);
//    EXPECT_NEAR(max, 1.0, 0.01);
//    EXPECT_NEAR(avg, .5, 0.02);
//}

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
//
//TEST(DualScattering, SampleBackHemisphere) {
//    std::vector<Vector3f> backVectors;
//
//    backVectors.push_back(SampleBackHemisphere(-.5 * Pi, 0.0));
//    backVectors.push_back(SampleBackHemisphere(-.5 * Pi, 0.5));
//    backVectors.push_back(SampleBackHemisphere(-.5 * Pi, 1.0));
//    backVectors.push_back(SampleBackHemisphere(.0, 0.0));
//    backVectors.push_back(SampleBackHemisphere(.0, 0.5));
//    backVectors.push_back(SampleBackHemisphere(.0, 1.0));
//    backVectors.push_back(SampleBackHemisphere(.5 * Pi, 0.0));
//    backVectors.push_back(SampleBackHemisphere(.5 * Pi, 0.5));
//    backVectors.push_back(SampleBackHemisphere(.5 * Pi, 1.0));
//    backVectors.push_back(SampleBackHemisphere(Point2f(0.0, 0.0)));
//    backVectors.push_back(SampleBackHemisphere(Point2f(0.5, 0.5)));
//    backVectors.push_back(SampleBackHemisphere(Point2f(0.1, 0.9)));
//    backVectors.push_back(SampleBackHemisphere(Point2f(0.9, 0.1)));
//    backVectors.push_back(SampleBackHemisphere(Point2f(1.0, 1.0)));
//
//    for (const auto& w : backVectors) {
//        EXPECT_FLOAT_EQ(w.Length(), 1.0);
//        EXPECT_LE(w.x, 0.0);
//    }
//}
