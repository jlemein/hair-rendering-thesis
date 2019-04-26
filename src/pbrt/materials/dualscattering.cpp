
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// materials/plastic.cpp*
#include "dualscattering.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "materials/marschner.h"
#include "hairutil.h"
#include "hair.h"
#include "samplers/random.h"
#include "sampling.h"
#include <random>
#include <mutex>

namespace pbrt {

    // DualscatteringMaterial Method Definitions

    void DualscatteringMaterial::ComputeScatteringFunctions(
            SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
            bool allowMultipleLobes) const {

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);

        MarschnerBSDF* marschnerBSDF = mMarschnerMaterial->CreateMarschnerBSDF(si, arena);
        DualScatteringBSDF* dualScatteringBSDF =
                ARENA_ALLOC(arena, DualScatteringBSDF)(*si, mEta, marschnerBSDF,
                mAlphaR, mAlphaTT, mAlphaTRT, mBetaR, mBetaTT, mBetaTRT);


        //static DualScatteringLookup* lookup = DualScatteringLookup::Get(arena, dualScatteringBSDF);
        //ARENA_ALLOC(arena, DualScatteringLookup)(dualScatteringBSDF);

        si->bsdf->Add(dualScatteringBSDF);
    }

    DualscatteringMaterial *CreateDualscatteringMaterial(const TextureParams &mp) {
        MarschnerMaterial* marschnerMaterial = CreateMarschnerMaterial(mp);

        // TODO: Check if dual scattering really wants to square the widths
        //        marschnerMaterial->mBr = Sqr(marschnerMaterial->mBr);
        //        marschnerMaterial->mBtt = Sqr(marschnerMaterial->mBtrt);
        //        marschnerMaterial->mBtrt = Sqr(marschnerMaterial->mBtrt);

        Float eta = 1.55;

        return new DualscatteringMaterial(eta, marschnerMaterial,
                marschnerMaterial->mAr, marschnerMaterial->mAtt, marschnerMaterial->mAtrt,
                marschnerMaterial->mBr, marschnerMaterial->mBtt, marschnerMaterial->mBtrt);
    }

    static Float Transmittance() {

    }

    static Float Variance() {

    }

    /*******************************
     * DualScatteringBSDF
     *******************************/

    DualScatteringBSDF::DualScatteringBSDF(const SurfaceInteraction& si, Float eta,
            MarschnerBSDF* marschnerBSDF,
            Float alphaR, Float alphaTT, Float alphaTRT,
            Float betaR, Float betaTT, Float betaTRT)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
    mEta(eta),
    mMarschnerBSDF(marschnerBSDF),
    mPosition(si.p),
    mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT),
    mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT), mLookup(DualScatteringLookup::Get(this)) {

    };

    /**
     *
     * @param eccentricity Eccentricity value (1 means circular)
     * @param etaPerp The index of refraction of the fiber
     * @param thetaH The longitudinal half angle
     * @return the adjusted index of refraction taking into account eccentricity
     */
    static Float EtaEccentricity(Float eccentricity, Float etaPerp, Float thetaH) {
        Float eta1 = 2.0 * (etaPerp - 1.0) * Sqr(eccentricity) - etaPerp + 2.0;
        Float eta2 = 2.0 * (etaPerp - 1.0) / Sqr(eccentricity) - etaPerp + 2.0;

        return 0.5 * ((eta1 + eta2) + cos(2.0 * thetaH)*(eta1 - eta2));
    }

    Spectrum DualScatteringBSDF::f(const Vector3f &wo, const Vector3f &wi) const {

        Float thetaI, phiI, thetaR, phiR;
        ToSphericalCoords(wi, thetaI, phiI);
        ToSphericalCoords(wo, thetaR, phiR);

        Float thetaD = DifferenceAngle(thetaI, thetaR);
        Float phi = RelativeAzimuth(phiI, phiR);
        Float thetaH = HalfAngle(thetaI, thetaR);
        Float phiH = HalfAngle(phiI, phiR);

        // Global multiple scattering
        // TODO: Check if wd = wi?
        const Vector3f& wd = wi;
        static GlobalScatteringInformation gsi;
        GatherGlobalScatteringInformation(wd, gsi);

        // Shading

        // compute local multiple scattering contribution
        // TODO: Is this correct?


        Spectrum Ab = BackscatteringAttenuation(thetaD);
        Spectrum deltaB = BackscatteringMean(thetaD);

        Spectrum backScatterVariance = BackscatteringVariance(thetaD);
        Spectrum forwardScatterVariance = ForwardScatteringVariance(thetaD);
        Spectrum fBack = 2.0 * Ab * Gaussian(backScatterVariance + forwardScatterVariance,
                Spectrum(thetaD + thetaR) - deltaB) / (Pi * Sqr(cos(thetaD)));

        // compute BCSDF of the fiber due to direct illumination
        //TODO: find a way to use beta squared
        Spectrum fDirectS = mMarschnerBSDF->f(wo, wi);
        Spectrum FDirect = gsi.directIlluminationFraction * (fDirectS + mDb * fBack);
        //
        //        // Compute BCSDF of the fiber due to forward scattered illumination similarly
        //        Spectrum fScatterS = EvaluateForwardScatteredMarschner(thetaR, thetaH, thetaD, phi, forwardScatterVariance);
        //        //        //TODO: Does this Pi belong here, or is it a typo in paper??
        //        Spectrum FScatter = (gsi.transmittance - gsi.directIlluminationFraction) * mDf * (fScatterS + Pi * mDb * fBack);

        // combine direct and forward scattered components
        //return (FDirect + FScatter) * cos(thetaI);
        return FDirect * cos(thetaI);
    }

    void DualScatteringBSDF::GatherGlobalScatteringInformation(const Vector3f& wd, GlobalScatteringInformation& gsi) const {
        //if (this->FindScatteringCount(wd) == 0) {
        gsi.directIlluminationFraction = 1.0;
        gsi.transmittance = 1.0;
        gsi.variance = 0.0;
        //        } else {
        //            gsi.directIlluminationFraction = 0.0;
        //            gsi.transmittance = ForwardTransmittance();
        //            gsi.variance = Variance();
        //        }
    }

    int DualScatteringBSDF::FindScatteringCount(const Vector3f& wd) const {
        //return VoxelGrid.interpolateP(mPosition, wd);
        return 2;
    }

    Spectrum DualScatteringBSDF::EvaluateForwardScatteredMarschner(Float thetaR, Float thetaH, Float thetaD,
            Float phi, Float forwardScatteredVariance) const {

        // TODO: check if we can remove eta perp calculation and let marschner
        // be responsible for it
        Float etaPerp, etaPar;
        ToBravais(mEta, thetaR, etaPerp, etaPar);

        // take into account eccentricity (by again adjusting the eta)
        Float eccentricity = mMarschnerBSDF->getEccentricity();
        if (eccentricity != 1.0) {
            etaPerp = EtaEccentricity(eccentricity, etaPerp, thetaH);
        }

        // Very similar to Marschner, but MG_r, MG_tt and MG_trt are used
        // to incorporate the forward scattering variance
        Float sinThetaR = sin(thetaR);
        Float sinThetaT = sinThetaR / etaPerp;
        Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        Spectrum result = (
                MG_r(2.0 * thetaH, forwardScatteredVariance) * mMarschnerBSDF->N_r(phi, etaPerp)
                + MG_tt(2.0 * thetaH, forwardScatteredVariance) * mMarschnerBSDF->N_tt(phi, etaPerp, cosThetaT)
                + MG_trt(2.0 * thetaH, forwardScatteredVariance) * mMarschnerBSDF->N_trt(phi, etaPerp, cosThetaT)
                ) / CosineSquared(thetaD);

        return result;
    }

    Spectrum DualScatteringBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const {

        return mMarschnerBSDF->Sample_f(wo, wi, sample, pdf, sampledType);
    }

    Float DualScatteringBSDF::MG_r(Float theta, Float forwardScatteringVariance) const {
        // WATCHOUT: do not square betaR, already squared at material creation
        return Gaussian(mBetaR + forwardScatteringVariance, theta - mAlphaR);
    }

    Float DualScatteringBSDF::MG_tt(Float theta, Float forwardScatteringVariance) const {
        // WATCHOUT: do not square betaTT, already squared at material creation
        return Gaussian(mBetaTT + forwardScatteringVariance, theta - mAlphaTT);
    }

    Float DualScatteringBSDF::MG_trt(Float theta, Float forwardScatteringVariance) const {
        // WATCHOUT: do not square betaTRT, already squared at material creation
        return Gaussian(mBetaTRT + forwardScatteringVariance, theta - mAlphaTRT);
    }

    /**
     * Equation 5 in dual scattering approximation
     */
    Spectrum DualScatteringBSDF::ForwardScatteringTransmittance() const {
        //        Float af = mMarschnerBSDF->f(wo, wi) * cos(thetaD);
        //        return df * pow(af)
        return Spectrum(.0);
    }

    Float DualScatteringBSDF::ForwardScatteringVariance(Float theta) const {
        return 1.0;
    }

    DualScatteringLookup* DualScatteringLookup::instance = 0;

    std::mutex myMutex;

    const DualScatteringLookup& DualScatteringLookup::Get(const DualScatteringBSDF* bsdf) {

        // We are locking here, because the singleton requires initialization
        std::lock_guard<std::mutex> myLock(myMutex);

        if (DualScatteringLookup::instance == 0) {
            instance = new DualScatteringLookup(bsdf);
            instance->Init();
        }

        return *instance;
    }

    void DualScatteringLookup::Init() {
        this->PrecomputeAverageScatteringAttenuation();
    }

    void DualScatteringLookup::PrecomputeAverageScatteringAttenuation() {
        printf("PRECOMPUTING FOR LOOKUP TABLE SIZE: %d", LOOKUP_TABLE_SIZE);
        for (int i = 0; i < LOOKUP_TABLE_SIZE; ++i) {
            Float ratio = i / static_cast<Float> (LOOKUP_TABLE_SIZE - 1);
            Float thetaD = (-.5 + ratio) * Pi;

            mAverageBackwardScatteringAttenuation.push_back(mDualScatteringBSDF->AverageBackwardScatteringAttenuation(thetaD));
            mAverageForwardScatteringAttenuation.push_back(mDualScatteringBSDF->AverageForwardScatteringAttenuation(thetaD));

            mAverageForwardScatteringAlpha.push_back(mDualScatteringBSDF->AverageForwardScatteringAlpha(thetaD));
            mAverageBackwardScatteringAlpha.push_back(mDualScatteringBSDF->AverageBackwardScatteringAlpha(thetaD));

            mAverageForwardScatteringBeta.push_back(mDualScatteringBSDF->AverageForwardScatteringBeta(thetaD));
            mAverageBackwardScatteringBeta.push_back(mDualScatteringBSDF->AverageBackwardScatteringBeta(thetaD));
        }

        printf(" [DONE]\n");
    }

    Spectrum Lookup(Float value, Float min, Float max, const std::vector<Spectrum>& lookupTable) {
        Float range = max - min;
        int lookupIndex = round((lookupTable.size() - 1) * (value - min) / range);

        CHECK_GE(lookupIndex, 0);
        CHECK_LT(lookupIndex, lookupTable.size());

        return lookupTable[lookupIndex];
    }

    Spectrum DualScatteringLookup::AverageBackwardScatteringAttenuation(Float thetaD) const {
        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageBackwardScatteringAttenuation);
    }

    Spectrum DualScatteringLookup::AverageForwardScatteringAttenuation(Float thetaD) const {
        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageForwardScatteringAttenuation);
    }

    Spectrum DualScatteringLookup::AverageBackwardScatteringAlpha(Float thetaD) const {
        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageBackwardScatteringAlpha);
    }

    Spectrum DualScatteringLookup::AverageBackwardScatteringBeta(Float thetaD) const {
        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageBackwardScatteringBeta);
    }

    Spectrum DualScatteringLookup::AverageForwardScatteringAlpha(Float thetaD) const {
        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageForwardScatteringAlpha);
    }

    Spectrum DualScatteringLookup::AverageForwardScatteringBeta(Float thetaD) const {
        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageForwardScatteringBeta);
    }

    static Vector3f SampleFrontHemisphere(const Point2f& uv) {
        Vector3f w = UniformSampleSphere(uv);
        if (w.x < 0.0) {
            w.x *= -1.0;
        }
        return w;
    }

    /**
     * Samples in the front hemisphere for a given theta
     * @param theta Fixed theta to use
     * @param u random number between 0 and 1 for phi generation
     * @return
     */
    static Vector3f SampleFrontHemisphere(Float theta, Float u) {
        Float phi = (-1.0 + 2.0 * u) * Pi;
        Vector3f w = FromSphericalCoords(theta, phi);

        if (w.x < 0.0) {
            w.x *= -1.0;
        }
        return w;
    }

    static Vector3f SampleBackHemisphere(const Point2f& uv) {
        Vector3f w = UniformSampleSphere(uv);
        if (w.x > 0.0) {
            w.x *= -1.0;
        }
        return w;
    }

    static Vector3f SampleBackHemisphere(Float theta, Float u) {
        Float phi = (-1.0 + 2.0 * u) * Pi;
        Vector3f w = FromSphericalCoords(theta, phi);
        if (w.x > 0.0) {
            w.x *= -1.0;
        }
        return w;
    }

    Spectrum DualScatteringBSDF::AverageForwardScatteringAttenuation(Float thetaD) const {
        const int SAMPLES = 1000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        Spectrum sum(.0);
        for (int i = 0; i < SAMPLES; ++i) {

            const Vector3f wi = SampleBackHemisphere(thetaD, dis(gen));
            const Vector3f woForward = SampleFrontHemisphere(Point2f(dis(gen), dis(gen)));

            sum += mMarschnerBSDF->f(woForward, wi);
        }

        // TODO: is dividing by .5*Pi correct
        //return sum / (.5 * Pi);
        return sum / static_cast<Float> (SAMPLES);
    }

    Spectrum DualScatteringBSDF::AverageBackwardScatteringAttenuation(Float thetaD) const {
        const int SAMPLES = 1000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        Spectrum sum(.0);
        for (int i = 0; i < SAMPLES; ++i) {
            Float phi = (-.5 + dis(gen)) * Pi;

            const Vector3f wi = SampleBackHemisphere(thetaD, dis(gen));
            const Vector3f woBackward = SampleBackHemisphere(Point2f(dis(gen), dis(gen)));

            sum += mMarschnerBSDF->f(woBackward, wi);
        }

        // TODO: is dividing by .5*Pi correct
        return sum / static_cast<Float> (SAMPLES);
    }

    Spectrum DualScatteringBSDF::AverageForwardScatteringAlpha(Float thetaD) const {
        const int SAMPLES = 1000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        Spectrum sum(.0);
        Spectrum denominator(.0);

        for (int i = 0; i < SAMPLES; ++i) {
            Float phi = (-.5 + dis(gen)) * Pi;

            const Vector3f wi = SampleBackHemisphere(thetaD, dis(gen));
            const Vector3f woForward = SampleFrontHemisphere(Point2f(dis(gen), dis(gen)));

            Spectrum fR = mMarschnerBSDF->f_r(woForward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woForward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woForward, wi);

            sum += fR * mAlphaR + fTT * mAlphaTT + fTRT * mAlphaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::AverageBackwardScatteringAlpha(Float thetaD) const {
        const int SAMPLES = 1000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        Spectrum sum(.0), denominator(.0);

        for (int i = 0; i < SAMPLES; ++i) {
            Float phi = (-.5 + dis(gen)) * Pi;

            const Vector3f wi = SampleBackHemisphere(thetaD, dis(gen));
            const Vector3f woBackward = SampleBackHemisphere(Point2f(dis(gen), dis(gen)));

            Spectrum fR = mMarschnerBSDF->f_r(woBackward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woBackward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woBackward, wi);

            sum += fR * mAlphaR + fTT * mAlphaTT + fTRT * mAlphaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::AverageForwardScatteringBeta(Float thetaD) const {
        const int SAMPLES = 1000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        Spectrum sum(.0), denominator(.0);

        for (int i = 0; i < SAMPLES; ++i) {
            Float phi = (-.5 + dis(gen)) * Pi;

            const Vector3f wi = SampleBackHemisphere(thetaD, dis(gen));
            const Vector3f woForward = SampleFrontHemisphere(Point2f(dis(gen), dis(gen)));

            Spectrum fR = mMarschnerBSDF->f_r(woForward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woForward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woForward, wi);

            sum += fR * mBetaR + fTT * mBetaTT + fTRT * mBetaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::AverageBackwardScatteringBeta(Float thetaD) const {
        const int SAMPLES = 1000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        Spectrum sum(.0), denominator(.0);

        for (int i = 0; i < SAMPLES; ++i) {
            Float phi = (-.5 + dis(gen)) * Pi;

            const Vector3f wi = SampleBackHemisphere(thetaD, dis(gen));
            const Vector3f woBackward = SampleBackHemisphere(Point2f(dis(gen), dis(gen)));

            Spectrum fR = mMarschnerBSDF->f_r(woBackward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woBackward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woBackward, wi);

            sum += fR * mBetaR + fTT * mBetaTT + fTRT * mBetaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::BackscatteringAttenuation(Float thetaD) const {
        Spectrum afSquared = Sqr(mLookup.AverageForwardScatteringAttenuation(thetaD));
        Spectrum ab = mLookup.AverageBackwardScatteringAttenuation(thetaD);

        Spectrum abPow3 = Pow3(ab);
        Spectrum oneMinusAfSquared = Spectrum(1.0) - afSquared;
        Spectrum oneMinusAfSquaredPow3 = Sqr(oneMinusAfSquared) * oneMinusAfSquared;

        Spectrum A1 = ab * afSquared / oneMinusAfSquared;
        Spectrum A3 = abPow3 * afSquared / oneMinusAfSquaredPow3;

        return A1 + A3;
    }

    //! equation 16

    Spectrum DualScatteringBSDF::BackscatteringMean(Float thetaD) const {
        //TODO: how to compute alphaF and alphaB?
        // according to paper these are average forward scattering means based on BCSDF of hair fiber

        // average forward and backward scattering shifts, taken from BCSDF of the hair fiber
        Spectrum alphaF = mLookup.AverageForwardScatteringAlpha(thetaD);
        Spectrum alphaB = mLookup.AverageBackwardScatteringAlpha(thetaD);

        Spectrum ab = mLookup.AverageBackwardScatteringAttenuation(thetaD);
        Spectrum af = mLookup.AverageForwardScatteringAttenuation(thetaD);

        Spectrum oneMinSquaredAf = Spectrum(1.0) - Sqr(af);
        Spectrum oneMinSquaredAfSquared = Sqr(oneMinSquaredAf);
        Spectrum oneMinSquaredAfPow3 = Pow3(oneMinSquaredAf);

        Spectrum part1 = Spectrum(1.0) - (2.0 * Sqr(ab) / oneMinSquaredAfSquared);
        Spectrum part2 = (2.0 * oneMinSquaredAfSquared + 4.0 * Sqr(af) * Sqr(ab)) / oneMinSquaredAfPow3;

        return alphaB * part1 + alphaF * part2;
    }

    //! equation 17

    Spectrum DualScatteringBSDF::BackscatteringVariance(Float thetaD) const {
        //TODO: how to compute betaF and betaB?
        // according to paper these are average backward scattering means based on BCSDF of hair fiber
        Spectrum betaF = mLookup.AverageForwardScatteringBeta(thetaD);
        Spectrum betaB = mLookup.AverageBackwardScatteringBeta(thetaD);

        Spectrum ab = mLookup.AverageBackwardScatteringAttenuation(thetaD);
        Spectrum af = mLookup.AverageForwardScatteringAttenuation(thetaD);

        Spectrum sqrt1 = Sqrt(2.0 * Sqr(betaF) + Sqr(betaB));
        Spectrum sqrt2 = Sqrt(2.0 * Sqr(betaF) + 3.0 * Sqr(betaB));
        Spectrum abPow3 = Pow3(ab);

        Spectrum factor = Spectrum(1.0) + this->mDb * Sqr(af);
        Spectrum nom = ab * sqrt1 + abPow3 * sqrt2;
        Spectrum denom = ab + abPow3 * (2.0 * betaF + 3.0 * betaB);

        return factor * nom / denom;
    }

    Float DualScatteringBSDF::Pdf(const Vector3f &wo, const Vector3f & wi) const {

        return PiOver4;
    }

    std::string DualScatteringBSDF::ToString() const {
        return "DualScatteringBSDF";
    }

} // namespace pbrt
