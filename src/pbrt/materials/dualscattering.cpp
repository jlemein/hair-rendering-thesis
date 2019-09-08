#include "dualscattering.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "materials/marschner.h"
#include "hairutil.h"
#include "hair.h"
#include "light.h"
#include "samplers/random.h"
#include "sampling.h"
#include <random>
#include <fstream>
#include "OpenVdbReader.h"
#include "materials/dualscatteringlookup.h"
#include "scene.h"
#include "materials/hairutil.h"

#include <mutex>

#include <fstream>

namespace pbrt {

    static std::default_random_engine generator;
    static std::uniform_real_distribution<double> distribution(0.0, 0.9999999999999);

    // DualscatteringMaterial Method Definitions

    void DualscatteringMaterial::ComputeScatteringFunctions(
            SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
            bool allowMultipleLobes) const {

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);

        MarschnerBSDF* marschnerBSDF = mMarschnerMaterial->CreateMarschnerBSDF(si, arena);

        DualScatteringBSDF* dualScatteringBSDF =
                ARENA_ALLOC(arena, DualScatteringBSDF)(*si, this->scene, mEta, marschnerBSDF,
                mAlphaR, mAlphaTT, mAlphaTRT, mBetaR, mBetaTT, mBetaTRT,
                mDf, mDb, mScatterCount, mVoxelGridFileName, mUniformSampling);

        si->bsdf->Add(dualScatteringBSDF);
    }

    DualscatteringMaterial *CreateDualscatteringMaterial(const TextureParams &mp) {
        MarschnerMaterial* marschnerMaterial = CreateMarschnerMaterial(mp);

        // TODO: Check if dual scattering really wants to square the widths.
        // In marschner no squared widths are used.
        //        marschnerMaterial->mBr = Sqr(marschnerMaterial->mBr);
        //        marschnerMaterial->mBtt = Sqr(marschnerMaterial->mBtrt);
        //        marschnerMaterial->mBtrt = Sqr(marschnerMaterial->mBtrt);

        Float df = mp.FindFloat("df", 0.7);
        Float db = mp.FindFloat("db", 0.7);
        Float scatterCount = mp.FindFloat("scatterCount", -1.0);
        std::string vdbFileName = mp.FindString("vdbFileName", "voxelgrid.vdb");
        std::string samplingMethod = mp.FindString("samplingMethod", "uniform");
        bool uniformSampling = samplingMethod == "uniform";

        return new DualscatteringMaterial(marschnerMaterial->mEta, marschnerMaterial,
                marschnerMaterial->mAr, marschnerMaterial->mAtt, marschnerMaterial->mAtrt,
                marschnerMaterial->mBr, marschnerMaterial->mBtt, marschnerMaterial->mBtrt,
                df, db, scatterCount, vdbFileName, uniformSampling);
    }

    /*******************************
     * DualScatteringBSDF
     *******************************/

    DualScatteringBSDF::DualScatteringBSDF(const SurfaceInteraction& si,
            const Scene* scene,
            Float eta,
            MarschnerBSDF* marschnerBSDF,
            Float alphaR, Float alphaTT, Float alphaTRT,
            Float betaR, Float betaTT, Float betaTRT,
            Float df, Float db, Float scatterCount,
            std::string voxelGridFileName,
            bool uniformSampling)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
    mEta(eta), mDb(db), mDf(df),
    mScene(scene),
    mMarschnerBSDF(marschnerBSDF),
    mPosition(si.p),
    ns(si.shading.n),
    ng(si.n),
    ss(Normalize(si.shading.dpdu)),
    ts(Cross(ns, ss)),
    mWorldToObject(si.shape != 0 ? *si.shape->WorldToObject : Transform()),
    mObjectToWorld(si.shape != 0 ? *si.shape->ObjectToWorld : Transform()),
    mShape(si.shape),
    mObjectBound(si.shape != 0 ? si.shape->ObjectBound() : Bounds3f()), mWorldBound(si.shape != 0 ? si.shape->WorldBound() : Bounds3f()),
    mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT),
    mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT),
    mScatterCount(scatterCount),
    mVoxelGridFileName(voxelGridFileName),
    mSampleUniform(uniformSampling) {

        // lookup table initialized here, so that all members
        // of this reference are initialized
        mLookup = DualScatteringLookup::Get(this);
    };

    Vector3f DualScatteringBSDF::WorldToLocal(const Vector3f &v) const {
        return Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
    }

    Vector3f DualScatteringBSDF::LocalToWorld(const Vector3f &v) const {
        return Vector3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
                ss.y * v.x + ts.y * v.y + ns.y * v.z,
                ss.z * v.x + ts.z * v.y + ns.z * v.z);
    }

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
        return f(wo, wi, 0, 0);
    }

    Spectrum DualScatteringBSDF::Li(const Vector3f& wi, const Spectrum& Li, const Scene& scene, const VisibilityTester& visibilityTester) const {
        // Added this function to do occlusion checking in the material itself.
        // For Dualscattering, occlusion checking is not necessary, because we approximate the global illumination.
        return Li;
    }

    Spectrum DualScatteringBSDF::f(const Vector3f &woLocal, const Vector3f &wiLocal, const Scene* scene, const VisibilityTester* visibilityTester) const {

        //        const Vector3f wiWorld = LocalToWorld(wiLocal);
        //        const Vector3f woWorld = LocalToWorld(woLocal);

        const MarschnerAngles angles(woLocal, wiLocal, mEta, this->mMarschnerBSDF->getEccentricity());
        GlobalScatteringInformation gsi = GatherGlobalScatteringInformation(mScene, visibilityTester, wiLocal, angles.thetaD);

        Spectrum forwardTransmittance = gsi.transmittance;
        Spectrum Ab = BackscatteringAttenuation(angles.thetaD);
        Spectrum deltaB = BackscatteringMean(angles.thetaD);

        Spectrum backScatterVariance = BackscatteringVariance(angles.thetaD);
        Spectrum forwardScatterVariance = gsi.variance;
        CHECK_GE(forwardScatterVariance.y(), 0.0);

        Spectrum fDirectBack = 2.0 * Ab
                * Gaussian(Sqrt(backScatterVariance), -deltaB + angles.thetaH)
                / (Pi * Sqr(cos(angles.thetaD)));

        Spectrum fScatterBack = 2.0 * Ab
                * Gaussian(Sqrt(backScatterVariance + forwardScatterVariance), -deltaB + angles.thetaH)
                / (Pi * Sqr(cos(angles.thetaD)));

        Spectrum F;

        if (gsi.isDirectIlluminated) {
            Spectrum fDirectS = mMarschnerBSDF->f(woLocal, wiLocal);
            F = fDirectS + mDb * fDirectBack;
        } else {
            Spectrum fScatterS = EvaluateForwardScatteredMarschner(angles.thetaR,
                    angles.thetaH, angles.thetaD, angles.phi, forwardScatterVariance);
            F = forwardTransmittance * mDf * (fScatterS + Pi * mDb * fScatterBack);
        }

        CHECK_GE(F.y(), 0.0);
        CHECK(!F.HasNaNs());

        return F * std::fabs(cos(angles.thetaI));
    }

    std::once_flag flag1;

    GlobalScatteringInformation DualScatteringBSDF::GatherGlobalScatteringInformation(
            const Scene* scene, const VisibilityTester* visibilityTester,
            const Vector3f& wiLocal, Float thetaD) const {

        GlobalScatteringInformation gsi;

        // @ref Dual-Scattering Approximation, Zinke et al (2007), section 4.1.1.

        //        const Ray ray(mPosition, wd);
        //        bool isDirectIlluminated = visibilityTester.Unoccluded(scene);

        // Interpolates a ray through a voxel grid, where each voxel contains the hair density
        // thereby returning the approximated number of hair strands intersected.

        //        Point3f pObject = mWorldToObject(mPosition);
        //        const InterpolationResult interpolationResult = this->mLookup->getVdbReader()->interpolateToInfinity(pObject, -wd);
        //        Float scatterCount = interpolationResult.scatterCount;

        Float scatterCount;
        bool isDirectIlluminated = false;

        if (mScatterCount < 0.0) {
            Point3f pObject = mWorldToObject(mPosition);
            Vector3f wd = LocalToWorld(wiLocal);

            isDirectIlluminated = visibilityTester != 0 ?
                    visibilityTester->Unoccluded(*scene) : mScene->IntersectP(Ray(mPosition, wd));


            const InterpolationResult interpolationResult = this->mLookup->getVdbReader()->interpolateToInfinity(pObject, -wd);
            scatterCount = interpolationResult.scatterCount;
        } else {
            scatterCount = mScatterCount;

            std::call_once(flag1, [&](){
                std::cout << std::endl << "Fixed scatter count: set to " << scatterCount << std::endl;
            });
        }



        if (isDirectIlluminated || scatterCount <= .5) {
            gsi.isDirectIlluminated = true;
            gsi.transmittance = 1.0;
            gsi.variance = Spectrum(.0);
        } else {
            gsi.isDirectIlluminated = false;
            gsi.transmittance = this->ForwardScatteringTransmittance(scatterCount, thetaD);
            gsi.variance = this->ForwardScatteringVariance(scatterCount, thetaD);
        }

        return gsi;
    }

    /**
     * Equation 5 in dual scattering approximation
     * @param scatterCount the amount of intersections along the shadow ray
     * @param thetaD the average orientation (thetaD) for all intersection along the shadow ray
     */
    Spectrum DualScatteringBSDF::ForwardScatteringTransmittance(Float scatterCount, Float thetaD) const {

        Spectrum s = this->mLookup->AverageForwardScatteringAttenuation(thetaD);
        return mDf * Pow(s, scatterCount);
    }

    Spectrum DualScatteringBSDF::ForwardScatteringVariance(Float scatterCount, Float thetaD) const {

        return scatterCount * this->mLookup->AverageForwardScatteringBeta(thetaD);
    }

    Spectrum DualScatteringBSDF::EvaluateForwardScatteredMarschner(Float thetaR, Float thetaH, Float thetaD,
            Float phi, Spectrum forwardScatteredVariance) const {

        Float etaPerp, etaPar;
        ToBravais(mEta, thetaR, etaPerp, etaPar);

        Float eccentricity = mMarschnerBSDF->getEccentricity();
        if (eccentricity != 1.0) {
            etaPerp = EtaEccentricity(eccentricity, etaPerp, thetaH);
        }

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

    Spectrum DualScatteringBSDF::MG_r(Float theta, Spectrum forwardScatteringVariance) const {
        // WATCHOUT: do not square betaR, already squared at material creation

        return Gaussian(Spectrum(mBetaR) + forwardScatteringVariance, theta - mAlphaR);
    }

    Spectrum DualScatteringBSDF::MG_tt(Float theta, Spectrum forwardScatteringVariance) const {
        // WATCHOUT: do not square betaTT, already squared at material creation

        return Gaussian(Spectrum(mBetaTT) + forwardScatteringVariance, theta - mAlphaTT);
    }

    Spectrum DualScatteringBSDF::MG_trt(Float theta, Spectrum forwardScatteringVariance) const {
        // WATCHOUT: do not square betaTRT, already squared at material creation

        return Gaussian(Spectrum(mBetaTRT) + forwardScatteringVariance, theta - mAlphaTRT);
    }


    //! Equation 6 in dual scattering approximation

    Spectrum DualScatteringBSDF::AverageForwardScatteringAttenuation(Float thetaD) const {
        // this method is just for trying to integrate around sphere = 1
        Spectrum sum(.0);
        MyRandomSampler sampler(0.0, 1.0);

        Float cosThetaD = std::max(cos(thetaD), .0);

        const int SAMPLES = 10000;

        for (int i = 0; i < SAMPLES; ++i) {
            Vector3f wi = SampleBackHemisphere(thetaD, sampler.next());
            Vector3f wr = SampleFrontHemisphere(Point2f(sampler.next(), sampler.next()));

            Float thetaR, phiR;
            ToSphericalCoords(wr, thetaR, phiR);

            Float phiI;
            ToSphericalCoords(wi, thetaD, phiI);

            Spectrum r = 2.0 * Pi * mMarschnerBSDF->f(wr, wi) * cosThetaD;

            sum += r;
        }


        sum /= static_cast<Float> (SAMPLES);
        CHECK_GE(sum.y(), 0.0);

        return sum;
    }

    //! equation 12 in dual scattering approximation

    Spectrum DualScatteringBSDF::AverageBackwardScatteringAttenuation(Float thetaD) const {
        Spectrum integral(.0);
        Float thetaR, phiR;
        MyRandomSampler sampler(0.0, 1.0);

        const int SAMPLES = 2000;
        for (int i = 0; i < SAMPLES; ++i) {

            const Vector3f woBackward = SampleBackHemisphere(Point2f(sampler.next(), sampler.next()));
            ToSphericalCoords(woBackward, thetaR, phiR);

            // -> thetaI = thetaR - 2.0 * thetaD;
            Float thetaI = GetThetaIFromDifferenceAngle(thetaD, thetaR);
            const Vector3f wi = SampleBackHemisphere(thetaI, sampler.next());

            // remove below later
            Float phiI;
            ToSphericalCoords(wi, thetaI, phiI);
            integral += mMarschnerBSDF->f(woBackward, wi) * cos(thetaD);
        }

        return 2.0 * Pi * integral / static_cast<Float> (SAMPLES);
    }

    Spectrum DualScatteringBSDF::AverageForwardScatteringAlpha(Float thetaD) const {
        Spectrum sum(.0);
        Spectrum denominator(.0);
        Float thetaR, phiR;
        MyRandomSampler sampler(0.0, 1.0);

        const int SAMPLES = 1000;
        for (int i = 0; i < SAMPLES; ++i) {

            const Vector3f woForward = SampleFrontHemisphere(Point2f(sampler.next(), sampler.next()));
            ToSphericalCoords(woForward, thetaR, phiR);

            // -> thetaI = thetaR - 2.0 * thetaD;
            Float thetaI = GetThetaIFromDifferenceAngle(thetaD, thetaR);
            const Vector3f wi = SampleBackHemisphere(thetaI, sampler.next());

            Spectrum fR = mMarschnerBSDF->f_r(woForward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woForward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woForward, wi);

            sum += fR * mAlphaR + fTT * mAlphaTT + fTRT * mAlphaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::AverageBackwardScatteringAlpha(Float thetaD) const {

        Spectrum sum(.0), denominator(.0);
        Float thetaR, phiR;
        MyRandomSampler sampler(0.0, 1.0);

        const int SAMPLES = 1000;
        for (int i = 0; i < SAMPLES; ++i) {

            const Vector3f woBackward = SampleBackHemisphere(Point2f(sampler.next(), sampler.next()));
            ToSphericalCoords(woBackward, thetaR, phiR);

            // -> thetaI = thetaR - 2.0 * thetaD;
            Float thetaI = GetThetaIFromDifferenceAngle(thetaD, thetaR);
            const Vector3f wi = SampleBackHemisphere(thetaI, sampler.next());

            Spectrum fR = mMarschnerBSDF->f_r(woBackward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woBackward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woBackward, wi);

            sum += fR * mAlphaR + fTT * mAlphaTT + fTRT * mAlphaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::AverageForwardScatteringBeta(Float thetaD) const {

        Spectrum sum(.0), denominator(.0);
        Float thetaR, phiR;
        MyRandomSampler sampler(0.0, 1.0);

        const int SAMPLES = 1000;
        for (int i = 0; i < SAMPLES; ++i) {

            const Vector3f woForward = SampleFrontHemisphere(Point2f(sampler.next(), sampler.next()));
            ToSphericalCoords(woForward, thetaR, phiR);

            // -> thetaI = thetaR - 2.0 * thetaD;
            Float thetaI = GetThetaIFromDifferenceAngle(thetaD, thetaR);
            const Vector3f wi = SampleBackHemisphere(thetaI, sampler.next());

            Spectrum fR = mMarschnerBSDF->f_r(woForward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woForward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woForward, wi);

            sum += fR * mBetaR + fTT * mBetaTT + fTRT * mBetaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    Spectrum DualScatteringBSDF::AverageBackwardScatteringBeta(Float thetaD) const {

        Spectrum sum(.0), denominator(.0);
        Float thetaR, phiR;
        MyRandomSampler sampler(0.0, 1.0);

        const int SAMPLES = 1000;
        for (int i = 0; i < SAMPLES; ++i) {

            const Vector3f woBackward = SampleBackHemisphere(Point2f(sampler.next(), sampler.next()));
            ToSphericalCoords(woBackward, thetaR, phiR);

            // -> thetaI = thetaR - 2.0 * thetaD;
            Float thetaI = GetThetaIFromDifferenceAngle(thetaD, thetaR);
            const Vector3f wi = SampleBackHemisphere(thetaI, sampler.next());

            Spectrum fR = mMarschnerBSDF->f_r(woBackward, wi);
            Spectrum fTT = mMarschnerBSDF->f_tt(woBackward, wi);
            Spectrum fTRT = mMarschnerBSDF->f_trt(woBackward, wi);

            sum += fR * mBetaR + fTT * mBetaTT + fTRT * mBetaTRT;
            denominator += (fR + fTT + fTRT);
        }

        return sum / denominator;
    }

    //! equation 14 in dual scattering approximation

    Spectrum DualScatteringBSDF::BackscatteringAttenuation(Float thetaD) const {
        Spectrum afSquared = Sqr(mLookup->AverageForwardScatteringAttenuation(thetaD));
        Spectrum ab = mLookup->AverageBackwardScatteringAttenuation(thetaD);

        Spectrum abPow3 = Pow3(ab);
        Spectrum oneMinusAfSquared = Spectrum(1.0) - afSquared;
        Spectrum oneMinusAfSquaredPow3 = Sqr(oneMinusAfSquared) * oneMinusAfSquared;

        Spectrum A1 = ab * afSquared / oneMinusAfSquared;
        Spectrum A3 = abPow3 * afSquared / oneMinusAfSquaredPow3;

        return A1 + A3;
    }

    //! equation 16 in Dual Scattering Approximation

    Spectrum DualScatteringBSDF::BackscatteringMean(Float thetaD) const {

        // average forward and backward scattering shifts, taken from BCSDF of the hair fiber
        Spectrum alphaF = mLookup->AverageForwardScatteringAlpha(thetaD);
        Spectrum alphaB = mLookup->AverageBackwardScatteringAlpha(thetaD);

        Spectrum ab = mLookup->AverageBackwardScatteringAttenuation(thetaD);
        Spectrum af = mLookup->AverageForwardScatteringAttenuation(thetaD);

        Spectrum oneMinSquaredAf = Spectrum(1.0) - Sqr(af);
        Spectrum oneMinSquaredAfSquared = Sqr(oneMinSquaredAf);
        Spectrum oneMinSquaredAfPow3 = Pow3(oneMinSquaredAf);

        Spectrum part1 = Spectrum(1.0) - (2.0 * Sqr(ab) / oneMinSquaredAfSquared);
        Spectrum part2 = (2.0 * oneMinSquaredAfSquared + 4.0 * Sqr(af) * Sqr(ab)) / oneMinSquaredAfPow3;

        return alphaB * part1 + alphaF * part2;
    }

    //! equation 17 in Dual Scattering Approximation

    Spectrum DualScatteringBSDF::BackscatteringVariance(Float thetaD) const {

        Spectrum betaF = mLookup->AverageForwardScatteringBeta(thetaD);
        Spectrum betaB = mLookup->AverageBackwardScatteringBeta(thetaD);

        Spectrum ab = mLookup->AverageBackwardScatteringAttenuation(thetaD);
        Spectrum af = mLookup->AverageForwardScatteringAttenuation(thetaD);

        Spectrum sqrt1 = Sqrt(2.0 * Sqr(betaF) + Sqr(betaB));
        Spectrum sqrt2 = Sqrt(2.0 * Sqr(betaF) + 3.0 * Sqr(betaB));
        Spectrum abPow3 = Pow3(ab);

        Spectrum factor = Spectrum(1.0) + this->mDb * Sqr(af);
        Spectrum nom = ab * sqrt1 + abPow3 * sqrt2;
        Spectrum denom = ab + abPow3 * (2.0 * betaF + 3.0 * betaB);

        return factor * nom / denom;
    }

    Spectrum DualScatteringBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType * sampledType) const {
        if (this->mSampleUniform) {
            return UniformSample_f(wo, wi, pdf, sample);
        } else {
            Float u[6] = {sample.x, sample.y,
                (Float) distribution(generator), (Float) distribution(generator),
                (Float) distribution(generator), (Float) distribution(generator)};
            return DEonSample_f(wo, wi, pdf, u);
        }
    }

    Float DualScatteringBSDF::Pdf(const Vector3f &wo, const Vector3f & wi) const {
        if (this->mSampleUniform) {
            return UniformPdf(wo, wi);
        } else {
            return DEonPdf(wo, wi);
        }
    }

    std::string DualScatteringBSDF::ToString() const {
        return "DualScatteringBSDF";
    }

    static Float sampleSphericalGaussian(Float variance, Float u) {
        return variance * log(exp(1.0 / variance) - 2.0 * u * sinh(1.0 / variance));
    }

    static Float sampleThetaI(Float u1, Float u2, Float variance, Float thetaCone) {
        Float thetaAccent = .5 * Pi - thetaCone;

        Float cosTheta = sampleSphericalGaussian(variance, u1);
        Float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        Float a = cosTheta * cos(thetaAccent);
        Float b = sinTheta * cos(2.0 * Pi * u2) * sin(thetaAccent);

        return SafeASin(a + b);
    }

    static Float BoxMuller(Float u1, Float u2) {
        return sqrt(-2.0 * log(u1)) * cos(2.0 * Pi * u2);
    }

    static void RetrieveSpecularWeights(Float cosGammaI, Float gammaT, Float cosThetaT, const Spectrum& sigmaA, Float etaT, Float *weights) {
        //printf("\n\n1. Called with: cosGammaI: %f, cosThetaT: %f\n", cosGammaI, cosThetaT);
        Float gammaI = acos(cosGammaI);

        Float specAr = AttenuationSpecR(cosGammaI, etaT) / (2.0 * fabs(DPhiDh_R(gammaI)));
        Float specAtt = AttenuationSpecTT(cosGammaI, gammaT, cosThetaT, sigmaA, etaT) / fabs(2.0 * DPhiDh(1, gammaI, etaT));
        Float specAtrt = AttenuationSpecTRT(cosGammaI, gammaT, cosThetaT, sigmaA, etaT) / fabs(2.0 * DPhiDh(2, gammaI, etaT));
        Float specSum = specAr + specAtt + specAtrt;

        weights[0] = specAr / specSum;
        weights[1] = specAtt / specSum;
        weights[2] = specAtrt / specSum;

        CHECK_NEAR(weights[0] + weights[1] + weights[2], 1.0, 0.02);
    }

    static void RetrieveSpecularWeightsFromRelativePhi(Float dphi, Float cosThetaT, const Spectrum &sigmaA, Float etaT, Float* weights) {
        // in this case we have been given a delta phi between two directions.
        // To compute the weights

        //printf("2. Called with: cosThetaT: %f, dphi: %f\n", cosThetaT, dphi);
        Float specAr = 0.0f, specAtt = 0.0f, specAtrt = 0.0f;
        // R
        {
            Float gammaI = SolveGammaRoot_R(dphi);
            Float cosGammaI = cos(gammaI);
            specAr = AttenuationSpecR(cosGammaI, etaT)
                    / (2.0 * fabs(DPhiDh_R(gammaI)));
        }

        // TT
        {
            Float gammaI[3];
            int nRoots = SolveGammaRoots(1, dphi, etaT, gammaI);
            Float gammaT = GammaT(gammaI[0], etaT);
            Float cosGammaI = cos(gammaI[0]);
            //printf("2. TT gammaI: %f, cosGammaI: %f\n", gammaI[0], cosGammaI);
            specAtt = AttenuationSpecTT(cosGammaI, gammaT, cosThetaT, sigmaA, etaT)
                    / (fabs(2.0 * DPhiDh(1, gammaI[0], etaT)));
        }

        // TRT
        {
            Float roots[3];
            int nRoots = SolveGammaRoots(2, dphi, etaT, roots);
            specAtrt = 0.0f;
            //printf("nRoots TRT: %d\n", nRoots);

            for (int i = 0; i < nRoots; ++i) {
                Float gammaI = roots[i];
                Float gammaT = GammaT(gammaI, etaT);
                Float cosGammaI = cos(gammaI);

                //printf("2. TRT cosGammaI: %f\n\n\n", cosGammaI);
                specAtrt += AttenuationSpecTRT(cosGammaI, gammaT, cosThetaT, sigmaA, etaT)
                        / (fabs(2.0 * DPhiDh(2, gammaI, etaT)));
            }
        }

        Float specSum = specAr + specAtt + specAtrt;
        weights[0] = specAr / specSum;
        weights[1] = specAtt / specSum;
        weights[2] = specAtrt / specSum;

        CHECK_NEAR(weights[0] + weights[1] + weights[2], 1.0, 1e-5);
    }

    Spectrum DualScatteringBSDF::DEonSample_f(const Vector3f &wo, Vector3f *wi, Float *pdf, const Float u[6]) const {

        //1. An incoming direction is given
        Float thetaO, phiO;
        ToSphericalCoords(wo, thetaO, phiO);
        //Float cosThetaO = cos(thetaO);

        // 2. uniformly choose a random offset [-1, 1]
        Float h = 2.0 * u[0] - 1.0;
        if (!(h >= -1.0001 && h <= 1.0001)) {
            printf("ERROR: h (%f) is not valid, u[0] = %f\n", h, u[0]);
        }
        // TODO: check if etaT = 1.55, or should we adjust for theta angles that we dont know yet?
        Float etaT = 1.55; //this->mEta;
        Float gammaI = SafeASin(h);
        Float gammaT = SafeASin(h / etaT);
        Float cosGammaI = cos(gammaI);

        // TODO: is this correct? To use thetaI instead of thetaR
        Float sinThetaO = sin(thetaO);
        Float sinThetaT = sinThetaO / etaT;
        Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        // 3. Compute attenuations Aspec(h, p) for each lobe assuming
        // no deflection away from specular cone
        //TODO: no deflection away from specular cone (means probably beta = 0)
        //        Float specAr = AttenuationSpecR(cosGammaI, etaT);
        //        Float specAtt = AttenuationSpecTT(cosGammaI, gammaT, cosThetaT, mMarschnerBSDF->getSigmaA(), etaT);
        //        Float specAtrt = AttenuationSpecTRT(cosGammaI, gammaT, cosThetaT, mMarschnerBSDF->getSigmaA(), etaT);
        //        Float specSum = specAr + specAtt + specAtrt;
        //        Float w[3] = {specAr / specSum, specAtt / specSum, specAtrt / specSum};

        Float w[3];
        RetrieveSpecularWeights(cosGammaI, gammaT, cosThetaT, mMarschnerBSDF->getSigmaA(), etaT, w);

        // 4. We select a lobe in proportion to the specular attenuations
        Float alpha[3] = {mAlphaR, mAlphaTT, mAlphaTRT};
        Float beta[3] = {mBetaR, mBetaTT, mBetaTRT};

        Float lobeSelect = u[1];
        int p = lobeSelect < w[0] ? 0
                : lobeSelect < (w[0] + w[1]) ? 1 : 2;

        // 5. We now know which Mp function (of variance vp) to sample, giving theta thetaO and thetaD
        Float variance = beta[p] * beta[p];
        Float thetaI = sampleThetaI(u[2], u[3], variance, -thetaO + alpha[p]);

        // TODO: which is correct?
        Float thetaD = DifferenceAngle(thetaI, thetaO);
        //Float thetaD = .5 * (thetaI + thetaO);

        Float etaPerp = sqrt(mEta * mEta - sin(thetaD)) / cos(thetaD);

        // 6. We sample a random Gaussian variable g and compute the relative outgoing azimuth
        Float g = BoxMuller(u[4], u[5]);
        Float phi = Phi(p, gammaI, gammaT) + g * sqrt(variance); //sqrt(vp) == abs(beta)
        Float phiI = UnwrapPhi(phiO + phi);

        // Evaluation of the model with Ap replaced by wp
        if (std::isnan(thetaI)) {
            printf("thetaI: %f, phiI: %f, u2: %f, u3: %f, variance: %f, -thetaO: %f, alphaP: %f, p: %d, phiO: %f, phi: %f\n", thetaI, phiI, u[2], u[3], variance, thetaO, alpha[p], p, phiO, phi);
        }
        *wi = FromSphericalCoords(thetaI, phiI);
        //*pdf = this->mMarschnerBSDF->f_Weighted(wo, *wi, w[0], w[1], w[2]);
        *pdf = this->DEonPdf(wo, *wi);

        // 7. We return a sample weight
        //Float Ap = this->mMarschnerBSDF->f_p(p, wi, *wo).y();
        //Float weight = Ap / w[p];

        // TODO: what is a weight? Has it to do with pdf?
        //*pdf = 1.0 / weight; ///Ap / w[p];
        return f(wo, *wi);
    }



    Float DualScatteringBSDF::DEonPdf(const Vector3f &wo, const Vector3f &wi) const {
        Float thetaI, phiI, thetaR, phiR;
        ToSphericalCoords(wi, thetaI, phiI);
        ToSphericalCoords(wo, thetaR, phiR);

        Float theta_d = DifferenceAngle(thetaI, thetaR);
        Float dphi = RelativeAzimuth(phiI, phiR);
        Float theta_h = HalfAngle(thetaI, thetaR);
        Float phi_h = HalfAngle(phiI, phiR);

        //Float etaPerp, etaPar;
        //ToBravais(mEta, thetaR, etaPerp, etaPar);
        Float etaT = 1.55; //etaPerp;

        Float sinThetaR = sin(thetaR);
        Float sinThetaT = sinThetaR / etaT;
        Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        Float weights[3];
        if (!(dphi < 3*Pi && dphi > -3*Pi)) {
            printf("dPhi: %f -- phiI: %f -- phiR: %f, wi: %f %f %f -- wr: %f %f %f\n", dphi, phiI, phiR, wi.x, wi.y, wi.z, wo.x, wo.y, wo.z);
        }
        RetrieveSpecularWeightsFromRelativePhi(dphi, cosThetaT, mMarschnerBSDF->getSigmaA(), etaT, weights);

//        Float MR = this->mMarschnerBSDF->M_r(thetaI, thetaR);
//        Float MTT = this->mMarschnerBSDF->M_tt(thetaI, thetaR);
//        Float MTRT = this->mMarschnerBSDF->M_trt(thetaI, thetaR);
//        float R = (weights[0] * MR + weights[1] * MTT + weights[2] * MTRT)/CosineSquared(theta_d);

        Float r = this->mMarschnerBSDF->f_Weighted(wo, wi, weights[0], weights[1], weights[2]);

        return r;
    }


    //
    // Uniform sampling
    //

    Spectrum DualScatteringBSDF::UniformSample_f(const Vector3f &wo, Vector3f *wi, Float* pdf, const Point2f& u) const {
        Float theta = acos(2.0 * u[0] - 1.0);
        Float phi = 2.0 * Pi * u[1];

        Float x = sin(theta) * cos(phi);
        Float y = sin(theta) * sin(phi);
        Float z = cos(theta);

        *wi = Vector3f(x, y, z);
        *pdf = Inv4Pi;

        return f(wo, *wi);
    }

    Float DualScatteringBSDF::UniformPdf(const Vector3f& wo, const Vector3f& wi) const {
        return Inv4Pi;
    }

} // namespace pbrt

//
//Float DualScatteringBSDF::PdfMarschnerR(const Vector3f &wo, const Vector3f & wi) const {
//        MarschnerAngles angles = MarschnerAngles(wo, wi, mEta, mMarschnerBSDF->getEccentricity());
//        Float cosThetaI = std::max(1e-5, cos(angles.thetaI));
//        Float thetaS = angles.thetaH - mAlphaR;
//
//        Vector3f Rperp = Vector3f(.0, wo.y, wo.z);
//        Rperp /= Rperp.Length();
//        Vector3f Lperp = Vector3f(0.0, wi.y, wi.z);
//        Lperp /= Lperp.Length();
//
//        Float denom = -.5 / Sqr(mBetaR);
//        Float M = exp(thetaS * thetaS * denom) / (mBetaR * sqrt(2.0 * Pi));
//        Float N = sqrt(std::max(0.0, .5 * (1.0 + Dot(Rperp, Lperp))));
//
//        Float pdf = M * N / (8.0 * cosThetaI);
//        CHECK_GE(pdf, 0.0);
//        return pdf;
//    }
