#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>

#include "pbrt.h"
#include "materials/marschner.h"
#include "materials/dualscattering.h"
#include "materials/dualscatteringlookup.h"
#include "materials/hairutil.h"
#include "sampling.h"
using namespace pbrt;

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);

const Float eta = 1.55;

Float r() {
    return distribution(generator);
}

void plotResponse(DualScatteringBSDF* bsdf, std::string fileName, int nSamples = 1000);
void plotAbAf(DualScatteringBSDF* bsdf, std::string fileName, int nSamples = 1000);

void createThetaIData(DualScatteringBSDF* bsdf, std::string outFileName, int nSamples = 1000);
void createPhiData(DualScatteringBSDF* bsdf, std::string outFileName, int nSamples = 1000);

DualScatteringBSDF* createDualScattering(Spectrum sigmaA, Float scatterCount,
        Float alphaR = -7.5,
        Float betaR = 7.5, Float betaTT = 3.725, Float betaTRT = 15.0) {
    Float alphaTT = -.5 * alphaR;
    Float alphaTRT = -1.5 * alphaR;
    Float alpha[3] = {Radians(alphaR), Radians(alphaTT), Radians(alphaTRT)};
    Float beta[3] = {Radians(betaR), Radians(betaTT), Radians(betaTRT)};

    Float hairRadius = 1.0;
    Float eccentricity = 0.9;

    Float glintScaleFactor = 0.4;
    Float causticWidth = 1.5;
    Float causticFade = 0.3;
    Float causticIntensityLimit = 0.5;

    SurfaceInteraction si = SurfaceInteraction();
    si.p = Point3f(0.0, 0.0, 0.0);
    si.shading.n = Normal3f(.0, 1.0, .0);
    si.n = Normal3f(.0, 1.0, .0);
    si.shading.dpdu = Vector3f(1.0, .0, .0);
    si.shading.dpdv = Vector3f(.0, .0, 1.0);

    MarschnerBSDF* marschner = new MarschnerBSDF(si,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            hairRadius, eta, sigmaA, eccentricity,
            glintScaleFactor, causticWidth, causticFade, causticIntensityLimit);

    Float df = 0.7, db = 0.7;

    return new DualScatteringBSDF(si, (Scene*) 0,
            eta, marschner,
            alpha[0], alpha[1], alpha[2], beta[0], beta[1], beta[2],
            df, db, scatterCount, "unnamed.vdb", true);
}

const Float rgbBrunette[3] = {0.432, 0.612, 0.98};
const Float rgbBlonde[3] = {0.15, 0.20, 0.30};
const Float rgbWhite[3] = {0.0, 0.0, 0.0};

const Spectrum BLONDE_HAIR = Spectrum::FromRGB(rgbBlonde);
const Spectrum BRUNETTE_HAIR = Spectrum::FromRGB(rgbBrunette);
const Spectrum WHITE_HAIR = Spectrum::FromRGB(rgbWhite);

int main(int argc, char** argv) {

    DualScatteringBSDF* blonde0 = createDualScattering(BLONDE_HAIR, 0);
    plotResponse(blonde0, "blonde0.data", 1000);
    plotAbAf(blonde0, "blonde_abaf.data");
    //    DualScatteringLookup::Reset();

    DualScatteringBSDF* blonde1 = createDualScattering(BLONDE_HAIR, 1);
    plotResponse(blonde1, "blonde1.data", 1000);
    //    DualScatteringLookup::Reset();

    DualScatteringBSDF* blonde2 = createDualScattering(BLONDE_HAIR, 2);
    plotResponse(blonde2, "blonde2.data", 1000);
    DualScatteringLookup::Reset();

    DualScatteringBSDF* brunette0 = createDualScattering(BRUNETTE_HAIR, 0, 3.0, 14, 8, 22);
    plotResponse(brunette0, "brunette0.data", 1000);
    plotAbAf(brunette0, "brunette_abaf.data");
    //    DualScatteringLookup::Reset();

    DualScatteringBSDF* brunette1 = createDualScattering(BRUNETTE_HAIR, 1, 3.0, 14, 8, 22);
    plotResponse(brunette1, "brunette1.data", 1000);
    //    DualScatteringLookup::Reset();

    DualScatteringBSDF* brunette2 = createDualScattering(BRUNETTE_HAIR, 2, 3.0, 14, 8, 22);
    plotResponse(brunette2, "brunette2.data", 1000);
    DualScatteringLookup::Reset();

    DualScatteringBSDF* white = createDualScattering(WHITE_HAIR, 0, 3.0, 14, 8, 22);
    plotAbAf(white, "white_abaf.data");
    DualScatteringLookup::Reset();

    //    int sampleCount = argc >= 2 ? std::stoi(std::string(argv[1])) : 1000;

    //    createThetaIData(dualScattering0, "thetaI_0.data", sampleCount);
    //    createThetaIData(dualScattering1, "thetaI_1.data", sampleCount);
    //    createThetaIData(dualScattering2, "thetaI_2.data", sampleCount);
    //    createThetaIData(dualScattering4, "thetaI_4.data", sampleCount);
    //
    //    createPhiData(dualScattering0, "phiI_0.data", sampleCount);
    //    createPhiData(dualScattering1, "phiI_1.data", sampleCount);
    //    createPhiData(dualScattering2, "phiI_2.data", sampleCount);
    //    createPhiData(dualScattering4, "phiI_4.data", sampleCount);


    return 0;
}

void plotAbAf(DualScatteringBSDF* bsdf, std::string fileName, int nSamples) {
    std::ofstream out(fileName.c_str());
    if (out.fail()) {
        std::cout << "Failed to open filename " << fileName << " for writing ab and af data\n";
        return;
    }

    for (int i = 0; i < nSamples; ++i) {
        Float thetaD = -.5 * Pi + (i / (nSamples - 1.0)) * Pi;
        Float af = bsdf->AverageForwardScatteringAttenuation(thetaD).y();
        Float ab = bsdf->AverageBackwardScatteringAttenuation(thetaD).y();
        out << thetaD << " " << ab << " " << af << std::endl;
    }
    out.close();
}

void plotResponse(DualScatteringBSDF* bsdf, std::string fileName, int nSamples) {
    std::vector<Vector3f> samplesAroundSphere;

    std::ofstream out(fileName.c_str());
    if (out.fail()) {
        std::cout << "Failed to open file " << fileName << " for writing\n";
        return;
    }

    for (int n = 0; n < nSamples; ++n) {
        const Point2f u = Point2f(r(), r());
        samplesAroundSphere.push_back(UniformSampleSphere(u));
    }

    out << 100 << " " << 100 << std::endl;

    const int THETA_SAMPLES = 100;
    const int PHI_SAMPLES = 100;

    for (int i = 0; i < THETA_SAMPLES; ++i) {
        for (int j = 0; j < PHI_SAMPLES; ++j) {
            Float thetaI = -.5 * Pi + i * Pi / (THETA_SAMPLES - 1.0);
            Float phiI = -Pi + j * 2.0 * Pi / (PHI_SAMPLES - 1.0);
            Vector3f wi = FromSphericalCoords(thetaI, phiI);

            Spectrum integral = .0;
            for (Vector3f wo : samplesAroundSphere) {
                integral += bsdf->f(wo, wi, 0, 0) * cos(thetaI);
            }

            integral *= 4.0 * Pi / samplesAroundSphere.size();
            std::cout << thetaI << " " << phiI << " " << integral.y() << std::endl;
            out << integral.y() << std::endl;
        }
    }
}

void createThetaIData(DualScatteringBSDF* bsdf, std::string outFileName, int nSamples) {
    std::vector<Vector3f> samples;

    std::ofstream out(outFileName.c_str());
    if (out.fail()) {
        std::cout << "Failed to open file " << outFileName << " for writing\n";
        return;
    }

    for (int n = 0; n < nSamples; ++n) {
        const Point2f u = Point2f(r(), r());
        samples.push_back(UniformSampleSphere(u));
    }

    for (int i = 0; i < nSamples; ++i) {
        Spectrum sum = .0;
        Float thetaI = -.5 * Pi + i * Pi / (nSamples - 1.0);

        for (int j = 0; j < nSamples; ++j) {
            Float phiI = -Pi + j * 2.0 * Pi / (nSamples - 1.0);
            Vector3f wi = FromSphericalCoords(thetaI, phiI);

            for (Vector3f wo : samples) {
                sum += bsdf->f(wo, wi, 0, 0) * cos(thetaI);
            }
        }

        // normalize
        std::cout << "dividing by sample count : " << nSamples * samples.size() << std::endl;
        sum = 4.0 * Pi * sum / static_cast<Float> (nSamples * samples.size());
        std::cout << "thetaI: " << thetaI << " -- r: " << sum.y() << std::endl;
        out << sum.y() << std::endl;
    }

    out.close();
}

void createPhiData(DualScatteringBSDF* bsdf, std::string outFileName, int nSamples) {
    std::vector<Vector3f> samples;

    std::ofstream out(outFileName.c_str());
    if (out.fail()) {
        std::cout << "Failed to open file " << outFileName << " for writing\n";
        return;
    }

    for (int n = 0; n < nSamples; ++n) {
        const Point2f u = Point2f(r(), r());
        samples.push_back(UniformSampleSphere(u));
    }

    for (int i = 0; i < nSamples; ++i) {
        Spectrum sum = .0;
        Float phiI = -Pi + i * 2.0 * Pi / (nSamples - 1.0);

        for (int j = 0; j < nSamples; ++j) {
            Float thetaI = -.5 * Pi + j * Pi / (nSamples - 1.0);

            Vector3f wi = FromSphericalCoords(thetaI, phiI);

            for (Vector3f wo : samples) {
                sum += bsdf->f(wo, wi, 0, 0) * cos(thetaI);
            }
        }

        // normalize
        sum = 4.0 * Pi * sum / static_cast<Float> (nSamples * samples.size());
        std::cout << "phiI: " << phiI << " -- r: " << sum.y() << std::endl;
        out << sum.y() << std::endl;
    }

    out.close();
}

