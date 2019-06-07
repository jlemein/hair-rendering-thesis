/*
 * File:   hairutil.cpp
 * Author: jeffrey lemein
 *
 * Created on March 22, 2019, 2:46 PM
 */

#include "hair.h"
#include "sampling.h"
#include "hairutil.h"

namespace pbrt {

    /**
     * @param min minimum value for the range to sample from
     * @param max maximum value for the range to sample from
     */
    MyRandomSampler::MyRandomSampler(Float min, Float max) {
        std::random_device rd;
        mDistribution = new std::uniform_real_distribution<Float>(min, max);
        mGenerationStrategy = new std::mt19937(rd());
    }

    /**
     * @returns next sample within the bounds defined in constructor
     */
    Float MyRandomSampler::next() const {
        return (*mDistribution)(*mGenerationStrategy);
    }

    bool CHECK_SPECTRUM_GE(const Spectrum& Ab, Float value) {
        static Float rgb[3];
        Ab.ToRGB(rgb);
        CHECK_GE(rgb[0], value);
        CHECK_GE(rgb[1], value);
        CHECK_GE(rgb[2], value);
    }

    /**
     * Safe acos function, to prevent 'nan' as result when parameter x is
     * (slightly) larger than or smaller than 1.0 or -1.0 respectively
     * (due to rounding errors).
     * @param x Argument to the acosine function. Values smaller than
     * -1.0 and larger than 1.0 are mapped to -1.0 and 1.0 respectively.
     * @return The inverse cosine result for parameter x.
     */
    Float SafeACos(Float x) {
        if (x <= -1.0) {
            return Pi;
        }
        if (x >= 1.0) {
            return 0.0;
        }
        return acos(x);
    }

    /**
     * Outputs the RGB values of the spectrum prefixed by a message
     * @param s Reference to a PBRT Spectrum
     * @param str Optional descriptive text to use as a prefix
     */
    void PrintSpectrum(const char* str, const Spectrum& sp) {
        Float rgb[3];
        sp.ToRGB(rgb);
        printf("%s [%f %f %f]\n", str, rgb[0], rgb[1], rgb[2]);
    }

    Spectrum Sqr(const Spectrum& s) {
        return s*s;
    }

    Spectrum Pow3(const Spectrum& s) {
        return s * s*s;
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
    Float AssurePositiveNonZero(Float value) {
        return std::max(Float(1e-5), value);
    }

    void ToSphericalCoords(const Vector3f& w, Float& theta, Float& phi) {
        theta = PiOver2 - SafeACos(w.x);
        phi = atan2(w.y, w.z);

        CHECK_LE(abs(theta), PiOver2);
        CHECK_LE(abs(phi), Pi);
    }

    Vector3f FromSphericalCoords(Float theta, Float phi) {
        Float x = cos(PiOver2 - theta);
        Float y = sin(phi) * cos(theta);
        Float z = cos(phi) * cos(theta);

        return Vector3f(x, y, z);
    }

    Float DifferenceAngle(Float theta_i, Float theta_r) {
        return 0.5 * (theta_r - theta_i);
    }

    Float GetThetaRFromDifferenceAngle(Float thetaD, Float thetaI) {
        return 2.0 * thetaD + thetaI;
    }

    Float GetThetaIFromDifferenceAngle(Float thetaD, Float thetaR) {
        return thetaR - 2.0 * thetaD;
    }

    // TODO: remove this function if DifferencePhi is almost similar

    Float RelativeAzimuth(Float phi_i, Float phi_r) {
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

    Float HalfAngle(Float a, Float b) {
        return 0.5 * (a + b);
    }

    Float CosineSquared(Float x) {
        return Sqr(cos(x));
    }

    Float SineSquared(Float x) {
        return Sqr(sin(x));
    }

    Float ClampPhi(Float phi, Float min, Float max) {
        while (phi > max) phi -= 2.0 * Pi;
        while (phi < min) phi += 2.0 * Pi;
        return phi;
    }

    /**
     * Returns the value of a normalized Gaussian function at point 'x', with a standard deviation of 'width'
     * @param width The width of the gaussian function
     * @param x Position to sample the gaussian
     * @returns Sampled position of the gaussian at position x
     */
    Float Gaussian(Float width, Float x) {
        Float a = 1.0 / (width * sqrt(2.0 * Pi));
        Float c = width; //width of the curve is beta (might also be 0.5 * sigma)

        Float nom = Sqr(x);
        Float den = 2.0 * Sqr(width);

        return a * exp(-nom / den);
    }

    Spectrum Gaussian(const Spectrum& width, const Spectrum& x) {
        Spectrum a = Spectrum(1.0) / (width * sqrt(2.0 * Pi));
        Spectrum c = width; //width of the curve is beta (might also be 0.5 * sigma)

        Spectrum nom = Sqr(x);
        Spectrum den = 2.0 * Sqr(width);

        return a * Exp(-nom / den);
    }

    Float SmoothStep(Float a, Float b, Float x) {
        CHECK_GT(b - a, 0.0);
        return Clamp((b - x) / (b - a), 0.0, 1.0);
    }

    /**
     * Slightly faster variant when you need to get both Bravais indices
     * @param eta Index of refraction
     * @param theta Longitudinal angle of incidence (in radians)
     * @param bravaisPerpendicular Output value that will hold perpendicular component of bravais index
     * @param bravaisParallel Output value that will hold parallel component of bravais index
     */
    void ToBravais(Float eta, Float theta, Float& etaPerpendicular, Float& etaParallel) {
        Float rootPart = sqrt(Sqr(eta) - SineSquared(theta));
        Float cosTheta = Clamp(cos(theta), 1e-5, 1.0);
        CHECK(cosTheta > 0.0 && cosTheta <= 1.0);

        etaPerpendicular = rootPart / cosTheta;
        etaParallel = Sqr(eta) * cosTheta / rootPart;
    }

    Float BravaisPerpendicular(Float eta, Float theta) {
        return sqrt(Sqr(eta) - Sqr(sin(theta))) / cos(theta);
    }

    Float BravaisParallel(Float eta, Float theta) {
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
    Float FresnelReflectionS(Float ni, Float nt, Float gamma_i) {
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
    Float FresnelReflectionP(Float ni, Float nt, Float gamma_i) {
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
    Float Fresnel(Float etaPerp, Float etaPar, Float gammaI) {
        CHECK_GT(etaPerp, 0.0);
        CHECK_GT(etaPar, 0.0);

        Float fresnelS = FresnelReflectionS(1.0, etaPerp, gammaI);
        CHECK(fresnelS >= 0.0 && fresnelS <= 1.0);

        // for s-polarized light
        Float fresnelP = FresnelReflectionP(1.0, etaPar, gammaI);
        CHECK(fresnelP > -0.00001 && fresnelP < 1.00001);

        return 0.5 * (fresnelP + fresnelS);
    }

    inline Float DiscriminantCardano(Float p, Float q) {
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
    int SolveDepressedCubic(Float a, Float c, Float d, Float roots[]) {
        Float p = c / a;
        Float q = d / a;

        Float D = DiscriminantCardano(p, q);

        if (D >= 0) {
            Float alpha = cbrt(-0.5 * q + sqrt(D));
            Float beta = cbrt(-0.5 * q - sqrt(D));

            roots[0] = alpha + beta;

            return 1;
        } else {
            CHECK_GE(0.25 * q * q - D, 0.0);
            Float R = 2.0 * cbrt(sqrt(0.25 * q * q - D));
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
     * Bounds gamma to it's valid range between [-pi/2; pi/2]
     * If gamma is higher or lower than this range, then gamma is
     * adjusted to fall back in this range.
     * The result is actually asin(sin(gamma)), but this leads to
     * precision errors.
     * @param gamma
     * @return
     */
    Float RangeBoundGammaInversed(Float gamma) {
        // wrap gamma back to range [-1.5pi ; pi/2]
        while (gamma > .5 * Pi) gamma -= 2 * Pi;
        while (gamma < -1.5 * Pi) gamma += 2 * Pi;

        // gamma >= pi/2, then gamma is valid and correct
        if (gamma >= -.5 * Pi) {
            return gamma;
        } else {
            // value below -pi/2 should be inversely mapped to [-pi/2 ; pi/2]
            return -Pi - gamma;
        }
    }

    Float RangeBoundGamma(Float gamma) {
        // wrap gamma back to range [-1.5pi ; pi/2]
        while (gamma > Pi) gamma -= 2.0 * Pi;
        while (gamma < -Pi) gamma += 2.0 * Pi;
        return gamma;
    }

    /**
     * Samples a vector in the front hemisphere.
     * This means theta can vary from -Pi/2 to Pi/2
     * Phi is bounded between [-Pi, -Pi/2] U [Pi/2, Pi]
     * @param uv
     * @return
     */
    Vector3f SampleBackHemisphere(const Point2f & uv) {
        Vector3f w = UniformSampleSphere(uv);

        if (w.z < 0.0) {
            w.z *= -1.0;
        }
        return w;
    }

    /**
     * Samples in the front hemisphere for a given theta
     * @param theta Fixed theta to use
     * @param u random number between 0 and 1 for phi generation
     * @return Unit vector in the front hemisphere
     */
    Vector3f SampleBackHemisphere(Float theta, Float u) {
        Float phi = (-1.0 + 2.0 * u) * Pi;
        Vector3f w = FromSphericalCoords(theta, phi);

        if (w.z < 0.0) {
            w.z *= -1.0;
        }
        return w;
    }

    /**
     * Samples in the front hemisphere for a given theta
     * @param uv 2D sample to be used to sample a vector on unit sphere
     * @return Unit vector in the back hemisphere
     */
    Vector3f SampleFrontHemisphere(const Point2f & uv) {
        Vector3f w = UniformSampleSphere(uv);
        if (w.z > 0.0) {

            w.z *= -1.0;
        }
        return w;
    }

    /**
     * Samples in the back hemisphere for a given theta
     * @param theta Fixed theta to use
     * @param u random number between 0 and 1 for phi generation
     * @return Unit vector in the back hemisphere
     */
    Vector3f SampleFrontHemisphere(Float theta, Float u) {
        Float phi = (-1.0 + 2.0 * u) * Pi;
        //Float phi = (-.5 + u) * Pi;
        Vector3f w = FromSphericalCoords(theta, phi);
        if (w.z > 0.0) {

            w.z *= -1.0;
        }
        return w;
    }

}