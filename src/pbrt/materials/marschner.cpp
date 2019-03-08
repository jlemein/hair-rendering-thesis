
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

namespace pbrt {

    /**
     * Safe acos function, to prevent 'nan' as result when parameter x is
     * (slightly) larger than or smaller than 1.0 or -1.0 respectively
     * (due to rounding errors).
     * @param x Argument to the acosine function. Values smaller than
     * -1.0 and larger than 1.0 are mapped to -1.0 and 1.0 respectively.
     * @return The inverse cosine result for parameter x.
     */
    static inline Float SafeACos(Float x) {
        if (x <= -1.0) {
            return Pi;
        }
        if (x >= 1.0) {
            return 0.0;
        }
        return acos(x);
    }

    static void PrintSpectrum(const char* s, const Spectrum& sp) {
        Float rgb[3];
        sp.ToRGB(rgb);
        printf("%s [%f %f %f]\n", s, rgb[0], rgb[1], rgb[2]);
    }

    /**
     * Utility function to adjust a value to be positive and nonzero.
     * Reason for this function is that rounding errors can lead to values
     * of zero or negative values close to zero.
     * This value in fact reaches the limit to zero. This function limits the
     * minium to 1e-5.
     * @param value Value to assure to be nonzero and positive
     * @return The input value or 1e-5
     */
    static inline Float AssurePositiveNonZero(Float value) {
        return std::max(Float(1e-5), value);
    }

    static void ToSphericalCoords(const Vector3f& w, Float& theta, Float& phi) {
        theta = PiOver2 - SafeACos(w.x);
        phi = atan2(w.y, w.z);

        CHECK_LE(abs(theta), PiOver2);
        CHECK_LE(abs(phi), Pi);
    }

    static Float DifferenceAngle(Float theta_i, Float theta_r) {
        return 0.5 * (theta_r - theta_i);
    }

    // TODO: remove this function if DifferencePhi is almost similar

    static Float RelativeAzimuth(Float phi_i, Float phi_r) {
        Float phi = phi_r - phi_i;

        // map to range [-pi, pi] so that our cubic solver works
        while (phi > Pi) {
            phi -= 2.0 * Pi;
        }
        while (phi < -Pi) {
            phi += 2.0 * Pi;
        }

        CHECK_LE(abs(phi), Pi);
        return phi;
    }

    Float DifferencePhi(Float phi1, Float phi2) {
        Float dphi = phi1 - phi2;
        while (dphi > Pi) {
            dphi -= 2 * Pi;
        }
        while (dphi < -Pi) {
            dphi += 2 * Pi;
        }
        return dphi;
    }

    static Float HalfAngle(Float a, Float b) {
        return 0.5 * (a + b);
    }

    static inline Float CosineSquared(Float x) {
        return Sqr(cos(x));
    }

    static inline Float SineSquared(Float x) {
        return Sqr(sin(x));
    }

    /**
     * Returns the value of a normalized Gaussian function at point 'x', with a standard deviation of 'width'
     * @param width The width of the gaussian function
     * @param x Position to sample the gaussian
     * @returns Sampled position of the gaussian at position x
     */
    static Float Gaussian(Float width, Float x) {
        Float a = 1.0 / (width * sqrt(2.0 * Pi));
        Float c = width; //width of the curve is beta (might also be 0.5 * sigma)

        Float nom = Sqr(x);
        Float den = 2.0 * Sqr(width);

        return a * exp(-nom / den);
    }

    /**
     * Slightly faster variant when you need to get both Bravais indices
     * @param eta Index of refraction
     * @param theta Longitudinal angle of incidence (in radians)
     * @param bravaisPerpendicular Output value that will hold perpendicular component of bravais index
     * @param bravaisParallel Output value that will hold parallel component of bravais index
     */
    static void ToBravais(Float eta, Float theta, Float& etaPerpendicular, Float& etaParallel) {
        Float rootPart = sqrt(Sqr(eta) - SineSquared(theta));

        // sin gamma = h, where 1 < h < 1
        // Gamma represents the angle between the incident ray and the normal of a dielectric cylinder
        // This means -pi <= gamma <= pi
        // This also means cos(gamma) is 0 <= cos(gamma) <= 1

        Float cosTheta = Clamp(cos(theta), 1e-5, 1.0);
        CHECK(cosTheta > 0.0 && cosTheta <= 1.0);

        etaPerpendicular = rootPart / cosTheta;
        etaParallel = Sqr(eta) * cosTheta / rootPart;
    }

    static Float BravaisPerpendicular(Float eta, Float theta) {
        return sqrt(Sqr(eta) - Sqr(sin(theta))) / cos(theta);
    }

    static Float BravaisParallel(Float eta, Float theta) {
        return Sqr(eta) * cos(theta) / sqrt(Sqr(eta) - Sqr(sin(theta)));
    }

    /**
     * Fresnel reflection for s-polarized light (perpendicular to incidence plane)
     * @param ni Index of refraction for material of incident side
     * @param nt Index of refraction for material on transmitted side
     * @param gamma_i Angle of incidence
     * @param gamma_t Angle of refraction
     * @return
     */
    static Float FresnelReflectionS(Float ni, Float nt, Float gamma_i) {
        Float sinGammaT = ni / nt * sin(gamma_i);

        // TODO: check if I can assume total reflection when "absolute" is bigger
        // than one. I did this, because sinGammaT is occasionally larger than 1
        // or smaller than -1.

        if (abs(sinGammaT) >= 1.0) {
            return 1.0;
        } else {
            Float gamma_t = SafeASin(sinGammaT);
            Float cos_gamma_i = cos(gamma_i);
            Float cos_gamma_t = cos(gamma_t);

            return Clamp(Sqr((ni * cos_gamma_i - nt * cos_gamma_t)
                    / (ni * cos_gamma_i + nt * cos_gamma_t)), 0.0, 1.0);

        }
    }

    /**
     * Fresnel reflection for p-polarized light (parallel to incidence plane)
     * @param ni Index of refraction for material of incident side
     * @param nt Index of refraction for material on transmitted side
     * @param gamma_i Angle of incidence
     * @param gamma_t Angle of refraction
     * @return
     */
    static Float FresnelReflectionP(Float ni, Float nt, Float gamma_i) {
        Float sinGammaT = ni / nt * sin(gamma_i);

        // TODO: check if I can assume total reflection when "absolute" is bigger
        // than one. I did this, because sinGammaT is occasionally larger than 1
        // or smaller than -1.
        // I think it can be because gamma_i represents the angle with the surface normal
        // and reflection is symmetric (if I am right?)
        // if gamma is +90 or -90 degrees, it doesn't matter for fresnel equations


        if (abs(sinGammaT) >= 1.0) {
            // total reflection
            return 1.0;
        } else {
            Float gamma_t = SafeASin(sinGammaT);
            Float cos_gamma_i = cos(gamma_i);
            Float cos_gamma_t = cos(gamma_t);
            return Clamp(Sqr((nt * cos_gamma_i - ni * cos_gamma_t) / (nt * cos_gamma_i + ni * cos_gamma_t)), 0.0, 1.0);
        }
    }

    /**
     * Fresnel reflection for polarized light
     * @param etaPerp
     * @param etaPar
     * @param gamma_i
     * @return
     */
    static Float Fresnel(Float etaPerp, Float etaPar, Float gammaI) {
        CHECK_GT(etaPerp, 0.0);
        CHECK_GT(etaPar, 0.0);

        // TODO: check if gamma_i should be projected for parallel and perpendicular planes
        // Marschner indicates that by using the altered index of refractions, that you can
        // just use gamma as angle (appendix B)

        // TODO: must we use gamma_t based on eta, or calculate gamma_t for projected and perpendicular cases
        // use Snell's law to find transmitted angle

        Float fresnelS = FresnelReflectionS(1.0, etaPerp, gammaI);
        CHECK(fresnelS >= 0.0 && fresnelS <= 1.0);

        // for s-polarized light
        Float fresnelP = FresnelReflectionP(1.0, etaPar, gammaI);
        CHECK(fresnelP > -0.00001 && fresnelP < 1.00001);

        return 0.5 * (fresnelP + fresnelS);
    }

    static inline Float DiscriminantCardano(Float p, Float q) {
        return 0.25 * q * q + p * p * p / 27.0;
    }

    /**
     * To solve a depressed cubic, which is a cubic polynomial of the form
     * ax^3 + cx + d = 0
     * @param a
     * @param c
     * @param c
     * @param roots The solved
     * @returns then number of roots found
     */
    static int SolveDepressedCubic(Float a, Float c, Float d, Float roots[]) {
        Float p = c / a;
        Float q = d / a;

        Float D = DiscriminantCardano(p, q);

        if (D >= 0) {
            Float alpha = pow(-0.5 * q + sqrt(D), Float(1.0 / 3.0));
            Float beta = -p / (3.0 * alpha);

            roots[0] = alpha + beta;
            return 1;
        } else {
            Float R = 2.0 * pow(sqrt(0.25 * q * q - D), 1.0 / 3.0);
            Float tanPhi = -2.0 * sqrt(-D) / q;
            Float phi = atan(tanPhi);

            // when phi is smaller than 0, then negate result
            if (phi < 0) {
                R *= -1.0;
            }

            roots[0] = R * cos(phi / 3.0);
            roots[1] = R * cos((phi + 2.0 * Pi) / 3.0);
            roots[2] = R * cos((phi + 4.0 * Pi) / 3.0);
            return 3;
        }
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

    static Float SolveGammaRoot_TT(Float phi, Float etaPerp) {
        CHECK_GT(etaPerp, 0.0);
        CHECK(1.0 / etaPerp >= 0.0 && 1.0 / etaPerp <= 1.0);
        Float constant = asin(1.0 / etaPerp);

        Float a = -8.0 * constant / (Pi * Pi * Pi);
        Float c = 6.0 * constant / Pi - 2.0;
        Float d = Pi - phi;

        // TODO: check if this is needed?
        //        while (d > Pi) d -= 2 * Pi;
        //        while (d < -Pi) d += 2 * Pi;

        Float roots[3];
        int numberRoots = SolveDepressedCubic(a, c, d, roots);
        CHECK_EQ(numberRoots, 1);

        return roots[0];
    }

    static int SolveGammaRoots(int p, Float phi, Float etaPerp, Float gammaRoots[3]) {
        CHECK_GT(etaPerp, 0.0);
        CHECK(1.0 / etaPerp >= 0.0 && 1.0 / etaPerp <= 1.0);
        Float C = asin(1.0 / etaPerp);

        Float a = -8.0 * p * C / (Pi * Pi * Pi);
        Float c = 6.0 * p * C / Pi - 2.0;
        Float d = p * Pi - phi;

        // TODO: check if this is needed?
        //        while (d > Pi) d -= 2 * Pi;
        //        while (d < -Pi) d += 2 * Pi;

        int numberRoots = SolveDepressedCubic(a, c, d, gammaRoots);
        CHECK(numberRoots == 1 || numberRoots == 3);

        return numberRoots;
    }

    static Spectrum Transmittance(const Spectrum& sigmaA, Float gammaT, Float cosThetaT) {

        // PBRT implementation
        Float cosGammaT = cos(gammaT);
        return Exp(-sigmaA * (2 * cosGammaT / cosThetaT));

        // My implementation, which is correct?
        // Float cosGamma2T = AssurePositiveNonZero(cos(2.0 * gammaT));
        // return Exp(-2.0 * (sigmaA / cosThetaT) * (1.0 + cosGamma2T));
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
     * Receives the incoming gammaI direction with the corresponding
     * refracted direction gammaT. It computes the resulting phi when scattered
     * through the cylinder.
     *
     * @param p
     * @param gammaI
     * @param gammaT
     * @return
     */
    Float Phi(int p, Float gammaI, Float gammaT) {
        Float phi = 2.0 * p * gammaT - 2.0 * gammaI + p * Pi;
        //        while (phi > Pi) phi -= 2.0 * Pi;
        //        while (phi < -Pi) phi += 2.0 * Pi;
        return phi;
    }

    Float PhiR(Float gammaI) {
        return -2.0 * gammaI;
    }

    /** Temporarily use Logistic from PBRT source code, to see if problems are gone */
    //    inline Float Logistic(Float x, Float s) {
    //        x = std::abs(x);
    //        return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
    //    }
    //
    //    inline Float LogisticCDF(Float x, Float s) {
    //        return 1 / (1 + std::exp(-x / s));
    //    }
    //
    //    inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
    //        CHECK_LT(a, b);
    //        return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
    //    }

    /** End of temporary PBRT code */

    //    static Float GlintLobe(Float phi) {
    //        Float s = 0.117160;
    //        return TrimmedLogistic(phi, s, -Pi, Pi);
    //        //return Clamp(cos(phi), 0.0, 1.0);
    //    }

    /*******************************
     * MarschnerMaterial
     *******************************/

    void MarschnerMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
            MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const {

        Spectrum sigmaA = mSigmaA->Evaluate(*si);

        //        Float rgb[3];
        //        sigmaA.ToRGB(rgb);
        //        printf("SigmaA: %f %f %f\n", rgb[0], rgb[1], rgb[2]);

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);

        Float h = 2.0 * si->uv[1] - 1.0;

        si->bsdf->Add(ARENA_ALLOC(arena, MarschnerBSDF)(*si, h, mAr, mAtt, mAtrt, mBr, mBtt, mBtrt, mEta, sigmaA));
    }

    MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp) {
        Float Ar = mp.FindFloat("Ar", Radians(-7.5));
        Float Br = mp.FindFloat("Br", Radians(7.5));
        Float hairRadius = 1.0;
        Float eta = 1.55;
        Float eccentricity = mp.FindFloat("eccentricity", Float(1.0));
        Float glintScaleFactor = 2.5;
        Float causticWidth = Radians(15.0);
        Float causticFade = 0.3;
        Float causticLimit = 0.5;

        Float rgb[3] = {0.432, 0.612, 0.98};

        std::shared_ptr<Texture < Spectrum>> sigmaA = mp.GetSpectrumTexture("sigmaA", Spectrum::FromRGB(rgb)); //should be defined as color
        std::shared_ptr<Texture < Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));

        return new MarschnerMaterial(Ar, Br, hairRadius, eta, eccentricity, glintScaleFactor, causticWidth, causticFade, causticLimit, sigmaA, Kd);
    }

    /*******************************
     * MarschnerBSDF
     *******************************/

    MarschnerBSDF::MarschnerBSDF(const SurfaceInteraction& si,
            Float h,
            Float alphaR, Float alphaTT, Float alphaTRT,
            Float betaR, Float betaTT, Float betaTRT,
            Float eta, Spectrum sigmaA
            )
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
    mH(h), mNs(si.shading.n), mNg(si.n), mDpdu(si.dpdu), mDpdv(si.dpdv),
    mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT),
    mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT),
    mEta(eta), mSigmaA(sigmaA) {

        CHECK(abs(mH) <= 1.0);
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

    Spectrum MarschnerBSDF::N_r(Float dphi, Float etaPerp, Float etaPar) const {

        // Compute the phi offset between incoming and reflected direction
        //Float offsetPhi = PhiR(gammaI);
        //Float dphi = DifferencePhi(phi, offsetPhi);

        // we only know phi_r at the moment [-pi, pi]
        // now find the gamma that contributes to this scattering direction
        Float gammaI = SolveGammaRoot_R(dphi);
        Float fresnel = FrDielectric(cos(gammaI), 1.0, etaPerp); //Fresnel(etaPerp, etaPar, gammaI);

        CHECK_EQ(PhiR(gammaI), dphi);

        // reflection is only determined by Fresnel
        return fresnel / (2.0 * fabs(DPhiDh_R(gammaI)));
    }

    Float inline GammaT(Float gammaI, Float etaPerp) {
        const Float Pi3 = Pi * Pi*Pi;
        const Float C = SafeASin(1.0 / etaPerp);

        return gammaI * 3.0 * C / Pi - gammaI * gammaI * gammaI * 4.0 * C / Pi3;
    }

    Spectrum MarschnerBSDF::N_tt(Float dphi, Float etaPerp, Float etaPar, Float cosThetaT) const {

        //Float dphi = DifferencePhi(phi, Phi(ScatteringMode::TT, gammaI, gammaT));
        Float gammaI = SolveGammaRoot_TT(dphi, etaPerp);
        Float gammaT = GammaT(gammaI, etaPerp);

        //CHECK_LT(fabs(dphi - Phi(1, gammaI, gammaT)), 0.01);

        Float fresnel = FrDielectric(cos(gammaI), 1.0, etaPerp); //Fresnel(etaPerp, etaPar, gammaI);

        return Sqr(1.0 - fresnel)
                * Transmittance(mSigmaA, gammaT, cosThetaT)
                / (fabs(2.0 * DPhiDh(ScatteringMode::TT, gammaI, etaPerp)));
    }

    Float smoothstep(Float a, Float b, Float x) {
        CHECK_GT(b - a, 0.0);
        return Clamp((b - x) / (b - a), 0.0, 1.0);
    }

    static Float ClampPhi(Float phi, Float min = -Pi, Float max = Pi) {
        while (phi > max) phi -= 2.0 * Pi;
        while (phi < min) phi += 2.0 * Pi;
        return phi;
    }

    Spectrum MarschnerBSDF::N_trt(Float phi, Float etaPerp, Float etaPar, Float cosThetaT) const {

        // Surface roughness parameters
        Float mCausticIntensityLimit = 0.5;
        Float mFadeRangeCausticMerge = 0.4; //[.2; .4]
        Float mCausticWidth = Radians(25.0); // between 10 and 25 degrees
        Float mGlintScaleFactor = 5; // between 0.5 to 5

        Float t;
        Float causticIntensity;
        Float phiC;

        // Compute Ntrt
        //
        Float roots[3];
        int nRoots = SolveGammaRoots(2, phi, etaPerp, roots);

        Spectrum sum(.0);

        for (int i = 0; i < nRoots; ++i) {
            Float gammaI = roots[i];
            Float gammaT = GammaT(gammaI, etaPerp);

            Float fresnel = FrDielectric(cos(gammaI), 1.0, etaPerp); //Fresnel(etaPerp, etaPar, gammaI);
            Float fresnelI = Fresnel(1.0 / etaPerp, 1.0 / etaPar, gammaT);
            //printf("fresnel: %f, fresnelI: %f  ---  gammaI: %f, gammaT: %f\n", fresnel, fresnelI, gammaI, gammaT);

            Spectrum T = Transmittance(mSigmaA, gammaT, cosThetaT);
            Spectrum Absorption = Sqr(1.0 - fresnel) * fresnel * T * T;
            Spectrum L = Absorption / (fabs(2.0 * DPhiDh(ScatteringMode::TRT, gammaI, etaPerp)));

            if (etaPerp < 2.0) {
                // root for caustic is (hc or -hc)
                Float hc = sqrt((4.0 - Sqr(etaPerp)) / 3.0);
                Float gammaC = SafeASin(hc);
                Float gammaTC = GammaT(gammaC, etaPerp);
                phiC = Phi(2, gammaI, gammaTC);

                //TODO: Check if this is correct
                Float squaredDPhiDH = Sqr(DPhiDh(2, gammaC, etaPerp));
                //Float squaredDPhiDH = Sqr(DPhiDh(2, gammaC*gammaC));

                causticIntensity = std::min(mCausticIntensityLimit, (Float) (2.0 * sqrt(2.0 * mCausticWidth / fabs(squaredDPhiDH))));
                t = 1.0;
            } else {
                phiC = 0.0;
                causticIntensity = mCausticIntensityLimit;
                t = smoothstep(2.0, 2.0 + mFadeRangeCausticMerge, etaPerp);
            }
            //printf("t = %f\n", t);

            // compute roughness
            Float gaussianCenter = Gaussian(mCausticWidth, .0);
            Float gaussianL = Gaussian(mCausticWidth, ClampPhi(phi - phiC));
            Float gaussianR = Gaussian(mCausticWidth, ClampPhi(phi + phiC));

            CHECK_GE(gaussianL, 0.0);
            CHECK_GE(gaussianR, 0.0);
            CHECK_GT(gaussianCenter, 0.0);
            CHECK(t >= 0.0 && t <= 1.0);

            //printf("phiC: %f\n", phiC);
            //PrintSpectrum("Np", L);
            L *= (1.0 - t * gaussianL / gaussianCenter);
            L *= (1.0 - t * gaussianR / gaussianCenter);
            L += t * mGlintScaleFactor * Absorption * causticIntensity * (gaussianL + gaussianR);
            //PrintSpectrum("L", L);
            sum += L;
        }

        //PrintSpectrum("sum", sum);
        return sum;
        //        return Spectrum(0.0);

    }


    //    Spectrum MarschnerBSDF::N_trt(Float phi, Float etaPerp, Float etaPar, Float gammaI, Float gammaT, Float sinGammaT, Float cosThetaT) const {
    //
    //        // find roots
    //        //Float dphi = DifferencePhi(phi, Phi(2.0, gammaI, gammaT));
    //
    //
    //        Float rGammaI[3];
    //        int nRoots = SolveGammaRoots(2, phi, etaPerp, rGammaI);
    //        Float fresnel = Fresnel(etaPerp, etaPar, gammaI);
    //        Float T = Transmittance(mSigmaA, gammaT, cosThetaT);
    //
    //        Spectrum sum(.0);
    //
    //        //Float rGammaT = SafeASin(sin(rGammaI[0])/etaPerp);
    //
    //        Spectrum A = Sqr(1.0 - fresnel) * fresnel * Sqr(T);
    //
    //        return A / (fabs(2.0 * DPhiDh(ScatteringMode::TRT, rGammaI[0], etaPerp)));
    //
    //    }

    Spectrum MarschnerBSDF::N_p(int p, Float relativePhi) const {

        return Spectrum(0.01);
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

        Float etaPerp, etaPar;

        // TODO: check if bravais index is based on theta_r or theta_i or maybe theta_d ??
        ToBravais(mEta, theta_r, etaPerp, etaPar);

        // Compute all useful properties
        Float gammaI = SafeASin(mH);
        Float sinGammaT = mH / etaPerp;
        Float cosGammaT = SafeSqrt(1.0 - Sqr(sinGammaT));
        Float gammaT = SafeASin(sinGammaT);

        Float sinThetaR = sin(theta_r);
        Float sinThetaT = sinThetaR / mEta;
        Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        Spectrum result = (
                M_r(theta_h) * N_r(phi, etaPerp, etaPar)
                + M_tt(theta_h) * N_tt(phi, etaPerp, etaPar, cosThetaT)
                + M_trt(theta_h) * N_trt(phi, etaPerp, etaPar, cosThetaT)
                ) / CosineSquared(theta_d);

        return result;
    }

    Spectrum MarschnerBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const {

        Float theta = Pi * sample.x;
        Float phi = 2.0 * Pi * sample.y;

        Float x = sin(theta) * cos(phi);
        Float y = sin(theta) * sin(phi);
        Float z = cos(theta);

        //TODO: sampling is not performed uniform, so pdf should be adjusted
        *pdf = this->Pdf(wo, *wi);
        *wi = Vector3f(x, y, z);

        //        return Spectrum(.0);
        return f(wo, *wi);
    }

    Float MarschnerBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {

        return PiOver4;
    }

    std::string MarschnerBSDF::ToString() const {
        return "MarschnerBSDF";
    }

} // namespace pbrt
