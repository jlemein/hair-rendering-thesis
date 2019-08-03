
#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "materials/dualscattering.cpp"
#include "materials/dualscatteringlookup.h"
#include <atomic>

#include <vector>

using namespace pbrt;

static double _integrateMonteCarloFront(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {

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

static double _integrateMonteCarloBack(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {

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

static Vector3f MyUniformSampleSphere(const Point2f &u) {
    Float z = 1 - 2 * u[0];
    Float r = std::sqrt(std::max((Float) 0, (Float) 1 - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

static double _integrateMonteCarlo(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // total area for unit sphere
    double V = 4.0 * Pi;

    double sum = .0;

    for (int i = 0; i < nSamples; ++i) {
        //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
        Vector3f wr = MyUniformSampleSphere(Point2f(distribution(generator), distribution(generator)));
        sum += f(wr, wi);
    }

    return V / static_cast<double> (nSamples) * sum;
}

static double _integrateMonteCarlo(std::function<double(const Float) > f, int nSamples) {

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

const Float eta = 1.55;
const Float eccentricity = 0.9;
const Vector3f alpha = Vector3f(0.0523599, -0.0261799, -0.0785398);
const Vector3f beta = Vector3f(0.244346, 0.139626, 0.383972);

const Float glintScale = 0.4;
const Float causticWidth = 0.174533;
const Float causticFade = 0.3;
const Float causticIntensityLimit = 0.5;
const Float hairRadius = 1.0;
const Float sigmaARgb[3] = {.0, .0, .0};
const Spectrum sigmaA = Spectrum::FromRGB(sigmaARgb);
SurfaceInteraction si = SurfaceInteraction();

//MarschnerBSDF* marschner = new MarschnerBSDF(si, alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2], hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);


//const Float eta = 1.55;
//const Float eccentricity = 1.0;
//const Vector3f alpha = Vector3f(1.0, 1.0, 1.0);
//const Vector3f beta = Vector3f(5.0, 2.5, 10.0);
//const Float glintScale = 0.5;
//const Float causticWidth = 1.5;
//const Float causticFade = 1.0;
//const Float causticIntensityLimit = 1.0;
//const Float hairRadius = 1.0;
//const Float sigmaARgb[3] = {0.4, 0.5, 0.1};
//const Spectrum sigmaA = Spectrum::FromRGB(sigmaARgb);
//SurfaceInteraction si = SurfaceInteraction();

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

void compareValues(MyRandomSampler& sampler, MarschnerBSDF* marschner, DualScatteringBSDF* dualScattering) {

    Vector3f wi = SampleBackHemisphere(0.0, sampler.next());
    Vector3f wo = SampleFrontHemisphere(0.0, sampler.next());

    Float rgbM[3];
    Float rgbD[3];
    marschner->f(wo, wi).ToRGB(rgbM);
    dualScattering->f(wo, wi).ToRGB(rgbD);

    printf("comparing marschner to dualscattering: [ %f %f %f ] vs [ %f %f %f ]\n",
            rgbM[0], rgbM[1], rgbM[2], rgbD[0], rgbD[1], rgbD[2]);
}

TEST(Dualscattering, CompareValues) {
    si.shading.dpdu = Vector3f(1.0, 0.0, 0.0);
    MarschnerBSDF* marschner = new MarschnerBSDF(si,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);

    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, (Scene*) 0,
            eta, marschner,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            1.0 /*0.7*/, 1.0 /*0.7*/, 0.0, "unnamed.vdb", true);

    MyRandomSampler sampler(0.0, 1.0);

    compareValues(sampler, marschner, dualScattering);
    compareValues(sampler, marschner, dualScattering);

}

TEST(ImportanceSampling, Prefer_TT_selection_when_h_close_to_zero) {
    MyRandomSampler sampler(0.0, 1.0);
    Point2f u = Point2f(sampler.next(), sampler.next());
    Vector3f wi = UniformSampleSphere(u);
    Float thetaI, phiI;
    ToSphericalCoords(wi, thetaI, phiI);

    Float h = 0.01;
    Float etaT = 1.55; //this->mEta;
    Float gammaI = SafeASin(h);
    Float gammaT = SafeASin(h / etaT);
    Float cosGammaI = cos(gammaI);

    // TODO: is this correct? To use thetaI instead of thetaR
    Float sinThetaI = sin(thetaI);
    Float sinThetaT = sinThetaI / etaT;
    Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    Float rgb[3] = {0.05, 0.15, 0.66};
    Spectrum sigA = Spectrum::FromRGB(rgb);


    Float wR = AttenuationSpec(0, cosGammaI, gammaT, cosThetaT, sigA, etaT);
    Float wTT = AttenuationSpec(1, cosGammaI, gammaT, cosThetaT, sigA, etaT);
    Float wTRT = AttenuationSpec(2, cosGammaI, gammaT, cosThetaT, sigA, etaT);
    Float wSum = wR + wTT + wTRT;
    wR /= wSum;
    wTT /= wSum;
    wTRT /= wSum;

    //printf("Probability R: %f proc, TT: %f proc, TRT: %f proc\n", wR, wTT, wTRT);

    CHECK_GT(wTT, wR);
    CHECK_GT(wTT, wTRT);

}

TEST(ImportanceSampling, Prefer_R_selection_when_h_close_to_one) {
    MyRandomSampler sampler(0.0, 1.0);
    Point2f u = Point2f(sampler.next(), sampler.next());
    Vector3f wi = UniformSampleSphere(u);
    Float thetaI, phiI;
    ToSphericalCoords(wi, thetaI, phiI);

    // Given a glancing angle, then the R lobe should have highest change of being selected
    Float h = 0.99;
    Float etaT = 1.55; //this->mEta;
    Float gammaI = SafeASin(h);
    Float gammaT = SafeASin(h / etaT);
    Float cosGammaI = cos(gammaI);

    // TODO: is this correct? To use thetaI instead of thetaR
    Float sinThetaI = sin(thetaI);
    Float sinThetaT = sinThetaI / etaT;
    Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    Float rgb[3] = {0.05, 0.15, 0.66};
    Spectrum sigA = Spectrum::FromRGB(rgb);


    Float wR = AttenuationSpec(0, cosGammaI, gammaT, cosThetaT, sigA, etaT);
    Float wTT = AttenuationSpec(1, cosGammaI, gammaT, cosThetaT, sigA, etaT);
    Float wTRT = AttenuationSpec(2, cosGammaI, gammaT, cosThetaT, sigA, etaT);
    Float wSum = wR + wTT + wTRT;
    wR /= wSum;
    wTT /= wSum;
    wTRT /= wSum;

    //printf("Probability R: %f proc, TT: %f proc, TRT: %f proc\n", wR, wTT, wTRT);

    CHECK_GT(wR, wTT);
    CHECK_GT(wR, wTRT);

}

TEST(DualScattering, ImportanceSamplingPdfIntegratesOne) {
    si.shading.dpdu = Vector3f(1.0, 0.0, 0.0);
    MarschnerBSDF* marschner = new MarschnerBSDF(si,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);

    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, (Scene*) 0,
            eta, marschner,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            1.0 /*0.7*/, 1.0 /*0.7*/, 0.0, "unnamed.vdb", true);

    const int SAMPLES_PHI = 100;
    const int SAMPLES_THETA = 100;

    // give a fixed incoming direction
    const Vector3f wi(.0f, 1.f, .0f);

    Float sum = .0f, reference = .0f;

    for (int i = 0; i < SAMPLES_THETA; ++i) {
        Float theta = (i / 100.0) * Pi - .5 * Pi;

        for (int j = 0; j < SAMPLES_PHI; ++j) {
            Float phi = (j / 100.0) * 2.0 * Pi;
            const Vector3f wo = FromSphericalCoords(theta, phi);

            //sum += bsdf->Pdf(wo, wi) * cos(theta);
            Float dPhi = 2.0 * Pi / SAMPLES_PHI;
            Float dTheta = Pi / SAMPLES_THETA;
            Float dA = dTheta * dPhi * cos(theta);

            //Float uniformPdf = dualScattering->Pdf(wo, wi);
            Float pdf = dualScattering->DEonPdf(wo, wi);
            ASSERT_GE(pdf, 0.0f);

            sum += dA * pdf;
            reference += dA * Inv4Pi;
        }
    }

    CHECK_NEAR(reference, 1.0, 0.001);
    CHECK_NEAR(sum, 1.0, 1e-2);

}
/*
TEST(Dualscattering, BsdfShouldBeEnergyConservant) {
    MyRandomSampler sampler(0.0, 1.0);

    si.shading.dpdu = Vector3f(1.0, 0.0, 0.0);
    MarschnerBSDF* marschner = new MarschnerBSDF(si,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);

    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, (Scene*) 0,
            eta, marschner,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            0.7, 0.7, 0.0, "unnamed.vdb");


    Spectrum sigmaNoAbsorption = Spectrum(.0);

    auto fnSphere = [&](const Vector3f wr, const Vector3f wi){
        return 1.0;
    };

    auto fnDualscattering = [&](const Vector3f wr, const Vector3f wi){
        Float cosAngle = Dot(wr, wi);
        return dualScattering->f(wr, wi).y(); // * fabs(cosAngle);
    };

    auto fnMarschner = [&](const Vector3f wr, const Vector3f wi){
        Float cosAngle = Dot(wr, wi);
        return marschner->f(wr, wi).y(); // * fabs(cosAngle);
    };

    //    auto fnSquared = [&](const Vector3f wr, const Vector3f wi){
    //        Float cosAngle = Dot(wr, wi);
    //        return marschnerSquared->f(wr, wi).y() * fabs(cosAngle);
    //    };
    //
    //    auto fnSqrt = [&](const Vector3f wr, const Vector3f wi){
    //        Float cosAngle = Dot(wr, wi);
    //        return marschnerSqrt->f(wr, wi).y() * fabs(cosAngle);
    //    };

    //    auto fnSphereSurface = [&](const Vector3f wr, const Vector3f wi){
    //        return 1.0;
    //    };

    Vector3f wii = SampleBackHemisphere(0.0, sampler.next());

    Float integratedSphere = _integrateMonteCarlo(wii, fnSphere, 10000);
    Float integratedMarschner = _integrateMonteCarlo(wii, fnMarschner, 10000);
    Float integratedDualscattering = _integrateMonteCarlo(wii, fnDualscattering, 10000);

    printf("Reference sphere integrated: %f (should be: %f)\n", integratedSphere, 4.0 * Pi);
    printf("Marschner integrated: %f (should be: %f)\n", integratedMarschner, 1.0);
    printf("DualScattering integrated: %f\n (should be: %f)", integratedDualscattering, 1.0);

    // expect that the integration are near the value (with an allowed offset of 0.05)
    EXPECT_NEAR(integratedSphere, 4.0 * Pi, 0.05);
    //    EXPECT_NEAR(integratedMarschner, 1.0, 0.05);
    //    EXPECT_LE(integratedDualscattering, 1.2);
    //    EXPECT_GE(integratedDualscattering, .8);
}
 */

//TEST(DualScattering, Pdf_MarschnerR_Must_integrate_to_1) {
//    si.shading.dpdu = Vector3f(1.0, 0.0, 0.0);
//
//    MarschnerBSDF* marschner = new MarschnerBSDF(si,
//            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
//            hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);
//
//    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, (Scene*) 0,
//            eta, marschner,
//            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
//            0.7, 0.7, 0.0, "unnamed.vdb");
//
//    const int SAMPLES = 1000;
//    const Vector3f wo = FromSphericalCoords(0.0, 0.0);
//    Float summedPdf = .0;
//
//    for (int i = 0; i < SAMPLES; ++i) {
//        Float theta = -.5 * Pi + Pi * (i / (SAMPLES - 1.0));
//        Float circumference = 2.0 * Pi * cos(theta) / SAMPLES;
//        Float s = .0;
//
//        for (int j = 0; j < SAMPLES; ++j) {
//
//            Float phi = -Pi + 2.0 * Pi * (j / (SAMPLES - 1.0));
//
//            Vector3f wi = FromSphericalCoords(theta, phi);
//            s += dualScattering->PdfMarschnerR(wo, wi);
//        }
//
//        summedPdf += s * circumference;
//    }
//
//    CHECK_NEAR(summedPdf / SAMPLES, 1.0, 1e-5);
//}


// is this test relevant?
/*
TEST(DualScattering, Pdf_MarschnerR_Must_integrate_to_10) {
    //    si.shading.dpdu = Vector3f(1.0, 0.0, 0.0);
    //
    //    MarschnerBSDF* marschner = new MarschnerBSDF(si,
    //            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
    //            hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);
    //
    //    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, (Scene*) 0,
    //            eta, marschner,
    //            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
    //            0.7, 0.7, 0.0, "unnamed.vdb");
    //
    //    const int SAMPLES = 100000;
    //    const Vector3f wo = FromSphericalCoords(1.0, 0.0);
    //    Float summedPdf = .0;
    //
    //    for (int i = 0; i < SAMPLES; ++i) {
    //        Vector3f wi = UniformSampleSphere(Point2f(dis(gen), dis(gen)));
    //        summedPdf += dualScattering->PdfMarschnerR(wo, wi);
    //    }
    //
    //    CHECK_NEAR(4.0 * Pi * summedPdf / SAMPLES, 1.0, 1e-5);
}

TEST(DualScattering, AverageBackwardScatteringAttenuation) {

    MarschnerBSDF* marschner = new MarschnerBSDF(si,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            hairRadius, eta, sigmaA, eccentricity, glintScale, causticWidth, causticFade, causticIntensityLimit);

    DualScatteringBSDF* dualScattering = new DualScatteringBSDF(si, (Scene*) 0,
            eta, marschner,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            0.7, 0.7, -1.0, "unnamed.vdb");

    const DualScatteringLookup* lookup = DualScatteringLookup::Get(dualScattering);

    const int SAMPLES = 10;
    for (int i = 0; i < SAMPLES; ++i) {
        Float thetaD = -.5 * Pi + (static_cast<Float> (i) / SAMPLES) * Pi;
        Spectrum response = dualScattering->AverageBackwardScatteringAlpha(thetaD);
        Float rgb[3];
        response.ToRGB(rgb);
        printf("thetaD: %f -- rgb: %f %f %f\n", thetaD, rgb[0], rgb[1], rgb[2]);

        EXPECT_GE(rgb[0], 0.0);
        EXPECT_GE(rgb[1], 0.0);
        EXPECT_GE(rgb[2], 0.0);
    }
}

//TEST(DualScattering, AverageForwardScatteringAttenuation) {
//
//    Spectrum theta0 = dualScattering->AverageForwardScatteringAttenuation(0.0);
//}

TEST(DualScattering, NumberGeneration) {


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

TEST(DualScattering, FromToSphericalCoordsVaryingTheta) {
    Float theta, phi, theta2, phi2;

    for (int i = 0; i < 10; ++i) {
        // we explicitly don't test theta= -/+ Pi, because phi can then be anything
        theta = -.49 * Pi + .98 * Pi * (i / 9.0);
        phi = 0.0;
        Vector3f w = FromSphericalCoords(theta, phi);
        ToSphericalCoords(w, theta2, phi2);

        //printf("w: %f %f %f (phi original: %f) --> returned phi: %f\n", w.x, w.y, w.z, phi, phi2);

        EXPECT_NEAR(theta, theta2, 1e-5);
        EXPECT_NEAR(phi, phi2, 1e-5);
    }
}

TEST(DualScattering, FromToSphericalCoordsVaryingPhi) {
    Float theta, phi, theta2, phi2;

    for (int i = 0; i < 10; ++i) {
        theta = 0.0;
        // we explicitly don't test phi = -/+ Pi, because of numeric issues, it can be 2Pi apart, which is still correct
        phi = -.99 * Pi + 1.98 * Pi * (i / 9.0);
        Vector3f w = FromSphericalCoords(theta, phi);
        ToSphericalCoords(w, theta2, phi2);

        EXPECT_NEAR(theta, theta2, 1e-5);
        EXPECT_NEAR(phi, phi2, 1e-5);
    }
}
//

TEST(DualScattering, SampleBackHemisphereToSphericalCoordsMustBeSame) {
    Float theta = 0.0;
    Float theta2, phi2;

    for (Float sample = 0.01; sample < .98; sample += 0.09) {
        theta = -.49 * Pi + sample * Pi;
        const Vector3f wBack = SampleBackHemisphere(theta, sample);

        ToSphericalCoords(wBack, theta2, phi2);

        EXPECT_FLOAT_EQ(theta, theta2);
        EXPECT_GE(phi2, -Pi);
        EXPECT_LE(phi2, Pi);
    }
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
        Float theta, phi;
        ToSphericalCoords(w, theta, phi);

        EXPECT_FLOAT_EQ(w.Length(), 1.0);
        EXPECT_GE(fabs(phi), .5 * Pi);
        EXPECT_LE(w.z, 0.0);
    }
}
//

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
        Float theta, phi;
        ToSphericalCoords(w, theta, phi);

        EXPECT_FLOAT_EQ(w.Length(), 1.0);
        EXPECT_LE(fabs(phi), .5 * Pi);
        EXPECT_GE(w.z, 0.0);
    }
}
 */