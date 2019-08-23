/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   hairutil.h
 * Author: jeffrey
 *
 * Created on March 22, 2019, 2:46 PM
 */

#ifndef HAIRUTIL_H
#define HAIRUTIL_H

#include "pbrt.h"
#include <random>

namespace pbrt {
    
    class MyRandomSampler {
    public:
        /**
         * @param min minimum value for the range to sample from
         * @param max maximum value for the range to sample from
         */
        MyRandomSampler(Float min, Float max);
        
        /**
         * @returns next sample within the bounds defined in constructor
         */
        Float next() const;
        
    private:
        std::mt19937 *mGenerationStrategy;
        std::uniform_real_distribution<Float> *mDistribution;
    };

    bool CHECK_SPECTRUM_GE(const Spectrum& Ab, Float value);
    
/**
     * Safe acos function, to prevent 'nan' as result when parameter x is
     * (slightly) larger than or smaller than 1.0 or -1.0 respectively
     * (due to rounding errors).
     * @param x Argument to the acosine function. Values smaller than
     * -1.0 and larger than 1.0 are mapped to -1.0 and 1.0 respectively.
     * @return The inverse cosine result for parameter x.
     */
    Float SafeACos(Float x) ;

    void PrintSpectrum(const char* str, const Spectrum& sp);
    
    Spectrum Sqr(const Spectrum& s);

    Spectrum Pow3(const Spectrum& s);
    
    Float Sign(Float x);

    /**
     * Utility function to adjust a value to be positive and nonzero.
     * Reason for this function is that rounding errors can lead to values
     * of zero or negative values close to zero.
     * This value in fact reaches the limit to zero. This function limits the
     * minium to 1e-5.
     * @param value Value to assure to be nonzero and positive
     * @return The input value or 1e-5
     */
    Float AssurePositiveNonZero(Float value) ;

    void ToSphericalCoords(const Vector3f& w, Float& theta, Float& phi);
    
    Vector3f FromSphericalCoords(Float theta, Float phi);

    Float DifferenceAngle(Float theta_i, Float theta_r);
    
    Float GetThetaIFromDifferenceAngle(Float thetaD, Float thetaR);
    Float GetThetaRFromDifferenceAngle(Float thetaD, Float thetaI);

    // TODO: remove this function if DifferencePhi is almost similar

    Float RelativeAzimuth(Float phi_i, Float phi_r);

    Float DifferencePhi(Float phi1, Float phi2) ;
    

    Float HalfAngle(Float a, Float b);

    Float CosineSquared(Float x) ;

    Float SineSquared(Float x) ;

    Float ClampPhi(Float phi, Float min = -Pi, Float max = Pi);

    /**
     * Returns the value of a normalized Gaussian function at point 'x', with a standard deviation of 'width'
     * @param width The width of the gaussian function
     * @param x Position to sample the gaussian
     * @returns Sampled position of the gaussian at position x
     */
    Float Gaussian(Float width, Float x);
    Spectrum Gaussian(const Spectrum& width, const Spectrum& x);
    
    Float SmoothStep(Float a, Float b, Float x);
    
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
    Float Phi(int p, Float gammaI, Float gammaT);
    Float PhiApprox(int p, Float gammaI, Float etaPerp);
    Float PhiR(Float gammaI);
    
    Float UnwrapPhi(Float phi);

    Float DPhiDh_R(Float gamma_i);
    Float DPhiDh(int p, Float gammaI, Float etaPerp);
    
    /**
     * Returns the specular amount
     * @param sigmaA
     * @param gammaT
     * @param cosThetaT
     * @return 
     */
    Float AttenuationSpec(int p, Float cosGammaI, Float gammaT, Float cosThetaT, const Spectrum& sigmaA, Float etaT);
    Float AttenuationSpecR(Float cosGammaI, Float etaT);
    Float AttenuationSpecTT(Float cosGammaI, Float gammaT, Float cosThetaT, const Spectrum& sigmaA, Float etaT);
    Float AttenuationSpecTRT(Float cosGammaI, Float gammaT, Float cosThetaT, const Spectrum& sigmaA, Float etaT);
    Float Attenuation(int p, Float cosGammaI);
    
    Spectrum Transmittance(const Spectrum& sigmaA, Float gammaT, Float cosThetaT);
    Spectrum TransmittanceR(Float cosTheta, Float etaT = 1.55, Float etaI = 1.0);
    Spectrum TransmittanceTT(Float cosTheta, Float etaT = 1.55, Float etaI = 1.0);
    Spectrum TransmittanceTRT(Float cosTheta, Float etaT = 1.55, Float etaI = 1.0);
    
    Float GammaT(Float gammaI, Float etaPerp);
    
    Float SolveGammaRoot_R(Float phi);
    int SolveGammaRoots(int p, Float phi, Float etaPerp, Float gammaRoots[3]);

    /**
     * Slightly faster variant when you need to get both Bravais indices
     * @param eta Index of refraction
     * @param theta Longitudinal angle of incidence (in radians)
     * @param bravaisPerpendicular Output value that will hold perpendicular component of bravais index
     * @param bravaisParallel Output value that will hold parallel component of bravais index
     */
    void ToBravais(Float eta, Float theta, Float& etaPerpendicular, Float& etaParallel) ;

    Float BravaisPerpendicular(Float eta, Float theta);

    Float BravaisParallel(Float eta, Float theta) ;

    /**
     * Fresnel reflection for s-polarized light (perpendicular to incidence plane)
     * @param ni Index of refraction for material of incident side
     * @param nt Index of refraction for material on transmitted side
     * @param gamma_i Angle of incidence
     * @param gamma_t Angle of refraction
     * @return
     */
    Float FresnelReflectionS(Float ni, Float nt, Float gamma_i);

    /**
     * Fresnel reflection for p-polarized light (parallel to incidence plane)
     * @param ni Index of refraction for material of incident side
     * @param nt Index of refraction for material on transmitted side
     * @param gamma_i Angle of incidence
     * @param gamma_t Angle of refraction
     * @return
     */
    Float FresnelReflectionP(Float ni, Float nt, Float gamma_i) ;

    /**
     * Fresnel reflection for polarized light
     * @param etaPerp
     * @param etaPar
     * @param gamma_i
     * @return
     */
    Float Fresnel(Float etaPerp, Float etaPar, Float gammaI) ;

    Float DiscriminantCardano(Float p, Float q) ;

    /**
     * To solve a depressed cubic, which is a cubic polynomial of the form
     * ax^3 + cx + d = 0
     * @param a
     * @param c
     * @param c
     * @param roots The solved
     * @returns then number of roots found
     */
    int SolveDepressedCubic(Float a, Float c, Float d, Float roots[]);

    /**
     * Bounds gamma to it's valid range between [-pi/2; pi/2]
     * If gamma is higher or lower than this range, then gamma is
     * adjusted to fall back in this range.
     * The result is actually asin(sin(gamma)), but this leads to
     * precision errors.
     * @param gamma
     * @return
     */
    Float RangeBoundGammaInversed(Float gamma);

    Float RangeBoundGamma(Float gamma);
    
    /**
     * Samples a vector in the front hemisphere.
     * This means theta can vary from -Pi/2 to Pi/2
     * Phi is bounded between [-Pi, -Pi/2] U [Pi/2, Pi]
     * @param uv
     * @return
     */
    Vector3f SampleBackHemisphere(const Point2f & uv);

    /**
     * Samples in the front hemisphere for a given theta
     * @param theta Fixed theta to use
     * @param u random number between 0 and 1 for phi generation
     * @return Unit vector in the front hemisphere
     */
    Vector3f SampleBackHemisphere(Float theta, Float u);

    /**
     * Samples in the front hemisphere for a given theta
     * @param uv 2D sample to be used to sample a vector on unit sphere
     * @return Unit vector in the back hemisphere
     */
    Vector3f SampleFrontHemisphere(const Point2f & uv);

    /**
     * Samples in the back hemisphere for a given theta
     * @param theta Fixed theta to use
     * @param u random number between 0 and 1 for phi generation
     * @return Unit vector in the back hemisphere
     */
    Vector3f SampleFrontHemisphere(Float theta, Float u);
    
}
    
#endif /* HAIRUTIL_H */

