
/*
 * File:   MarschnerMaterial.cpp
 * Author: jeffrey lemein
 *
 * Created on February 11, 2019, 11:34 PM
 */

#include "marschner.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "hair.h"

#include <algorithm>

#include "hairutil.h"

namespace pbrt {

    static int SolveGammaRoots(int p, Float phi, Float etaPerp, Float gammaRoots[3]) {
        CHECK_GT(etaPerp, 0.0);
        CHECK(1.0 / etaPerp >= 0.0 && 1.0 / etaPerp <= 1.0);

        Float C = asin(1.0 / etaPerp);

        Float a = -8.0 * p * C / (Pi * Pi * Pi);
        Float c = 6.0 * p * C / Pi - 2.0;
        Float d = (p % 2) * Pi - phi;

        // by wrapping like this, we have exactly one possible d value to solve for
        while (d > Pi) d -= 2 * Pi;
        while (d < -Pi) d += 2 * Pi;

        int numberRoots = SolveDepressedCubic(a, c, d, gammaRoots);
        CHECK(numberRoots == 1 || numberRoots == 3);

        // filter roots that are invalid
        int numberValidRoots = 0;
        for (int i = 0; i < numberRoots; ++i) {
            Float gamma = RangeBoundGamma(gammaRoots[i]);
            if (fabs(gamma) <= .5 * Pi) {
                gammaRoots[numberValidRoots++] = gamma;
            }
        }

        return numberValidRoots;
    }

    /**
     * Solves the root for reflection (R) mode of Marschner.
     * Gamma is returned instead of the root h (for efficiency reasons)
     *
     *  h = sin gamma, so if you want h, then do arcsin(gamma)
     *
     * @param phi
     * @return Gamma_i
     */
    static inline Float SolveGammaRoot_R(Float phi) {
        return -phi / 2.0;
    }

    static Float DPhiDh_R(Float gamma_i) {
        return -2.0 / AssurePositiveNonZero(sqrt(1.0 - SineSquared(gamma_i)));
    }

    static Float DPhiDh(int p, Float gammaI, Float etaPerp) {
        Float c = asin(1.0 / etaPerp);
        Float a = 8.0 * p * c / (Pi * Pi * Pi);
        Float b = 6.0 * p * c / Pi - 2.0;

        return (-3.0 * a * gammaI * gammaI + b) / SafeSqrt(1.0 - SineSquared(gammaI));
    }

    /**
     * Second derivative of root function
     * @param p The scattering mode (1=R, 2=TT, 3=TRT)
     * @param gammaI
     * @param etaPerp
     * @return
     */
    static Float DPhi2Dh2(int p, Float gammaI, Float etaPerp) {
        Float c = asin(1.0 / etaPerp);
        Float a = 8.0 * p * c / (Pi * Pi * Pi);
        return (-6.0 * a * gammaI) / (1.0 - Sqr(sin(gammaI)));
    }

    /**
     * @param etaPerp should be smaller than 2.0
     * @return
     */
    static Float GammaCaustic(Float etaPerp) {
        // root for caustic is (hc or -hc)
        Float hc = sqrt((4.0 - Sqr(etaPerp)) / 3.0);
        return SafeASin(-hc);
    }

    /**
     * @param etaPerp The index of refraction
     * @return the angle (in radians) where a caustic will appear
     */
    static Float PhiCaustic(Float etaPerp) {
        if (etaPerp < 2.0) {
            Float gammaCaustic = GammaCaustic(etaPerp);
            return Phi(2, gammaCaustic, GammaT(gammaCaustic, etaPerp));
        } else {
            return 0.0;
        }
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

    /*******************************
     * MarschnerMaterial
     *******************************/

    MarschnerBSDF* MarschnerMaterial::CreateMarschnerBSDF(SurfaceInteraction *si,
            MemoryArena &arena) const {

        Spectrum sigmaA = mHairRadius * mSigmaA->Evaluate(*si);

        return ARENA_ALLOC(arena, MarschnerBSDF)(*si, mAr, mAtt, mAtrt,
                mBr, mBtt, mBtrt,
                mHairRadius, mEta, sigmaA, mEccentricity, mGlintScaleFactor,
                mCausticWidth, mCausticFadeRange, mCausticIntensityLimit);
    }

    void MarschnerMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
            MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const {

        Spectrum sigmaA = mSigmaA->Evaluate(*si);

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);

        si->bsdf->Add(CreateMarschnerBSDF(si, arena));
    }

    MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp) {

        Float alphaR = mp.FindFloat("alphaR", Float(-7.5));
        Float alphaTT = mp.FindFloat("alphaTT", Float(-.5 * alphaR));
        Float alphaTRT = mp.FindFloat("alphaTRT", Float(-1.5 * alphaR));

        Float betaR = mp.FindFloat("betaR", Float(7.5));
        Float betaTT = mp.FindFloat("betaTT", Float(.5 * betaR));
        Float betaTRT = mp.FindFloat("betaTRT", Float(2.0 * betaR));

        Float alpha[3] = {Radians(alphaR), Radians(alphaTT), Radians(alphaTRT)};
        Float beta[3] = {Radians(betaR), Radians(betaTT), Radians(betaTRT)};

        Float hairRadius = mp.FindFloat("hairRadius", Float(1.0));
        Float eta = 1.55;
        Float eccentricity = mp.FindFloat("eccentricity", Float(1.0));

        Float glintScaleFactor = mp.FindFloat("glintScaleFactor", Float(2.5));
        //        if (glintScaleFactor < 0.5 || glintScaleFactor > 5.0)
        //            Warning("Glint scale factor should be between 0.5 and 5.0, but is %f", glintScaleFactor);

        Float causticWidth = mp.FindFloat("causticWidth", 15.0);
        //        if (causticWidth < 10.0 || causticWidth > 25.0)
        //            Warning("Caustic width should be between 10 and 25 degrees, but is %f degrees", causticWidth);
        Float causticWidthRadians = Radians(causticWidth);

        Float causticFade = mp.FindFloat("causticFadeRange", 0.3);
        if (causticFade < 0.2 || causticFade > 0.4) {
            Warning("caustic fade range should be between 0.2 and 0.4, but is %f", causticFade);
        }
        Float causticIntensityLimit = mp.FindFloat("causticIntensityLimit", Float(0.5));
        if (causticIntensityLimit < 0.0) {
            Error("caustic intensity limit should be positive, but is %f", causticIntensityLimit);
        } else if (causticIntensityLimit > 0.5) {
            Warning("caustic intensity limit should not be larger than 0.5, but is %f", causticIntensityLimit);
        }

        Float defaultSigmaA[3] = {0.432, 0.612, 0.98};

        std::shared_ptr<Texture < Spectrum>> sigmaA = mp.GetSpectrumTexture("sigmaA", Spectrum::FromRGB(defaultSigmaA)); //should be defined as color

        return new MarschnerMaterial(alpha, beta, hairRadius, eta, eccentricity, glintScaleFactor, causticWidthRadians, causticFade, causticIntensityLimit, sigmaA);
    }

    /*******************************
     * MarschnerBSDF
     *******************************/

    MarschnerBSDF::MarschnerBSDF(const SurfaceInteraction& si,
            Float alphaR, Float alphaTT, Float alphaTRT,
            Float betaR, Float betaTT, Float betaTRT,
            Float hairRadius,
            Float eta, Spectrum sigmaA, Float eccentricity,
            Float glintScaleFactor,
            Float causticWidth, Float causticFadeRange, Float causticIntensityLimit
            )
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
    mNs(si.shading.n), mNg(si.n), mDpdu(si.dpdu), mDpdv(si.dpdv),
    mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT),
    mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT),
    mHairRadius(hairRadius), mEta(eta), mSigmaA(sigmaA), mEccentricity(eccentricity), mGlintScaleFactor(glintScaleFactor),
    mCausticWidth(causticWidth), mCausticFadeRange(causticFadeRange), mCausticIntensityLimit(causticIntensityLimit) {
    };

    Float MarschnerBSDF::M_r(Float theta_h) const {

        return Gaussian(mBetaR, theta_h - mAlphaR);
    }

    Float MarschnerBSDF::M_tt(Float theta_h) const {

        return Gaussian(mBetaTT, theta_h - mAlphaTT);
    }

    Float MarschnerBSDF::M_trt(Float theta_h) const {

        return Gaussian(mBetaTRT, theta_h - mAlphaTRT);
    }

    Spectrum MarschnerBSDF::N_r(Float dphi, Float etaPerp) const {

        Float gammaI = SolveGammaRoot_R(dphi);
        Float fresnel = FrDielectric(cos(gammaI), 1.0, etaPerp);

        CHECK_EQ(PhiR(gammaI), dphi);

        // reflection is only determined by Fresnel (no absorption)
        return fresnel / (2.0 * fabs(DPhiDh_R(gammaI)));
    }

    Spectrum MarschnerBSDF::N_tt(Float dphi, Float etaPerp, Float cosThetaT) const {

        Float gammaI;
        int nRoots = SolveGammaRoots(1, dphi, etaPerp, &gammaI);
        if (nRoots == 0) {
            return Spectrum(.0);
        }

        CHECK_LE(nRoots, 1);
        CHECK_LE(fabs(gammaI), .5 * Pi);

        Float fresnel = FrDielectric(cos(gammaI), 1.0, etaPerp);

        return Sqr(1.0 - fresnel)
                * Transmittance(mSigmaA, GammaT(gammaI, etaPerp), cosThetaT)
                / (fabs(2.0 * DPhiDh(1, gammaI, etaPerp)));
    }

    Spectrum MarschnerBSDF::N_trt(Float phi, Float etaPerp, Float cosThetaT) const {

        Float roots[3];
        int nRoots = SolveGammaRoots(2, phi, etaPerp, roots);

        Spectrum sum(.0);

        for (int i = 0; i < nRoots; ++i) {
            Float gammaI = roots[i];
            Float gammaT = GammaT(gammaI, etaPerp);

            CHECK_LE(fabs(gammaI), .5 * Pi);
            CHECK_NEAR(Phi(2, gammaI, gammaT) - phi, 0.0, 1e-2);

            Float fresnel = FrDielectric(cos(gammaI), 1.0, etaPerp);
            Spectrum T = Transmittance(mSigmaA, gammaT, cosThetaT);
            Spectrum absorption = Sqr(1.0 - fresnel) * fresnel * T * T;
            Spectrum L = absorption / (fabs(2.0 * DPhiDh(2, gammaI, etaPerp)));

            // remove glints from response
            applySurfaceRoughness(L, absorption, phi, etaPerp);

            sum += L;
        }

        return sum;
    }

    void MarschnerBSDF::applySurfaceRoughness(Spectrum& L, const Spectrum& absorption, Float phi, Float etaPerp) const {

        Float t;
        Float gammaCaustic;
        Float causticIntensity;

        if (etaPerp < 2.0) {
            gammaCaustic = GammaCaustic(etaPerp);
            causticIntensity = std::min(mCausticIntensityLimit,
                    (Float) (2.0 * sqrt(2.0 * mCausticWidth / fabs(DPhi2Dh2(2, gammaCaustic, etaPerp)))));
            t = 1.0;
        } else {
            gammaCaustic = 0.0;
            causticIntensity = mCausticIntensityLimit;
            t = SmoothStep(2.0, 2.0 + mCausticFadeRange, etaPerp);
        }

        Float phiCaustic = PhiCaustic(etaPerp);
        Float gaussianL = Gaussian(mCausticWidth, ClampPhi(phi - phiCaustic));
        Float gaussianR = Gaussian(mCausticWidth, ClampPhi(phi + phiCaustic));
        Float gaussianCenter = Gaussian(mCausticWidth, .0);

        CHECK_GE(gaussianL, 0.0);
        CHECK_GE(gaussianR, 0.0);
        CHECK_GT(gaussianCenter, 0.0);
        CHECK(t >= 0.0 && t <= 1.0);

        L *= (1.0 - t * gaussianL / gaussianCenter);
        L *= (1.0 - t * gaussianR / gaussianCenter);
        L += t * mGlintScaleFactor * absorption * causticIntensity * (gaussianL + gaussianR);
    }

    Spectrum MarschnerBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
        // x axis goes with the fiber, from root to tip
        // y represents normal to hair fiber (major axis)
        // z axis represents normal to hair fiber (minor axis)

        Float theta_i, phi_i, theta_r, phi_r;
        ToSphericalCoords(wi, theta_i, phi_i);
        ToSphericalCoords(wo, theta_r, phi_r);

        Float theta_d = DifferenceAngle(theta_i, theta_r);
        Float phi = RelativeAzimuth(phi_i, phi_r);
        Float theta_h = HalfAngle(theta_i, theta_r);
        Float phi_h = HalfAngle(phi_i, phi_r);

        // If there is a longitudinal inclination, then the cross section is
        // elliptic, complicating our scattering model.
        // By adjusting the index of refraction (eta) using Bravais' law,
        // we have qualitatively the same effects as an elliptic cross section,
        // but we can now continue to use a circular cross section in the
        // Marschner scattering model.

        // TODO: check if bravais index is based on theta_r or theta_i or maybe theta_d ??
        Float etaPerp, etaPar;
        ToBravais(mEta, theta_r, etaPerp, etaPar);

        // take into account eccentricity (by again adjusting the eta)
        if (mEccentricity != 1.0) {
            etaPerp = EtaEccentricity(mEccentricity, etaPerp, theta_h);
        }

        Float sinThetaR = sin(theta_r);
        Float sinThetaT = sinThetaR / etaPerp;
        Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        Spectrum result = (
                M_r(2.0 * theta_h) * N_r(phi, etaPerp)
                + M_tt(2.0 * theta_h) * N_tt(phi, etaPerp, cosThetaT)
                + M_trt(2.0 * theta_h) * N_trt(phi, etaPerp, cosThetaT)
                ) / CosineSquared(theta_d);

        return result;
    }

    Spectrum MarschnerBSDF::f_p(int p, const Vector3f &wo, const Vector3f &wi) const {
        if (p == 0) {
            return f_r(wo, wi);
        } else if (p == 1) {
            return f_tt(wo, wi);
        } else {
            return f_trt(wo, wi);
        }
    }

    Spectrum MarschnerBSDF::f_r(const Vector3f &wo, const Vector3f &wi) const {
        const MarschnerAngles ma(wo, wi, mEta, mEccentricity);

        return M_r(2.0 * ma.thetaH) * N_r(ma.phi, ma.etaPerp) / CosineSquared(ma.thetaD);
    }

    Spectrum MarschnerBSDF::f_tt(const Vector3f &wo, const Vector3f &wi) const {
        const MarschnerAngles ma(wo, wi, mEta, mEccentricity);

        return M_tt(2.0 * ma.thetaH) * N_tt(ma.phi, ma.etaPerp, ma.cosThetaT) / CosineSquared(ma.thetaD);
    }

    Spectrum MarschnerBSDF::f_trt(const Vector3f &wo, const Vector3f &wi) const {
        const MarschnerAngles ma(wo, wi, mEta, mEccentricity);

        return M_trt(2.0 * ma.thetaH) * N_trt(ma.phi, ma.etaPerp, ma.cosThetaT) / CosineSquared(ma.thetaD);
    }

    Spectrum MarschnerBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const {

        Float theta = acos(2.0 * sample.x - 1.0);
        Float phi = 2.0 * Pi * sample.y;

        Float x = sin(theta) * cos(phi);
        Float y = sin(theta) * sin(phi);
        Float z = cos(theta);

        *wi = Vector3f(x, y, z);
        *pdf = this->Pdf(wo, *wi);

        return f(wo, *wi);
    }

    const Spectrum& MarschnerBSDF::getSigmaA() const {
        return mSigmaA;
    }

    Float MarschnerBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
        return PiOver4;
    }

    std::string MarschnerBSDF::ToString() const {
        std::stringstream ss;


        Float rgb[3];
        this->mSigmaA.ToRGB(rgb);

        ss << "alphaR: " << this->mAlphaR
                << " alphaTT: " << this->mAlphaTT
                << " alphaTRT: " << this->mAlphaTRT
                << " betaR: " << this->mBetaR
                << " betaTT: " << this->mBetaTT
                << " betaTRT: " << this->mBetaTRT
                << " causticFadeRange: " << this->mCausticFadeRange
                << " causticIntensityLimit: " << this->mCausticIntensityLimit
                << " causticWidth: " << this->mCausticWidth
                << " eccentricity: " << this->mEccentricity
                << " eta: " << this->mEta
                << " glintScaleFactor: " << this->mGlintScaleFactor
                << " hairRadius: " << this->mHairRadius
                << " sigmaA: " << rgb[0] << ", " << rgb[1] << ", " << rgb[2];

        return "MarschnerBSDF -> " + ss.str();
    }

    Float MarschnerBSDF::getEccentricity() const {
        return mEccentricity;
    }

    MarschnerAngles::MarschnerAngles(const Vector3f& wo, const Vector3f& wi, Float eta, Float eccentricity) {
        ToSphericalCoords(wi, this->thetaI, this->phiI);
        ToSphericalCoords(wo, this->thetaR, this->phiR);

        this->thetaD = DifferenceAngle(this->thetaI, this->thetaR);
        this->phi = RelativeAzimuth(this->phiI, this->phiR);
        this->thetaH = HalfAngle(this->thetaI, this->thetaR);
        this->phiH = HalfAngle(this->phiI, this->phiR);

        // If there is a longitudinal inclination, then the cross section is
        // elliptic, complicating our scattering model.
        // By adjusting the index of refraction (eta) using Bravais' law,
        // we have qualitatively the same effects as an elliptic cross section,
        // but we can now continue to use a circular cross section in the
        // Marschner scattering model.

        // TODO: check if bravais index is based on theta_r or theta_i or maybe theta_d ??
        ToBravais(eta, this->thetaR, this->etaPerp, this->etaPar);

        // take into account eccentricity (by again adjusting the eta)
        if (eccentricity != 1.0) {
            this->etaPerp = EtaEccentricity(eccentricity, etaPerp, this->thetaH);
        }

        this->sinThetaR = sin(this->thetaR);
        this->sinThetaT = sinThetaR / this->etaPerp;
        this->cosThetaT = SafeSqrt(1 - Sqr(this->sinThetaT));
    }

} // namespace pbrt
