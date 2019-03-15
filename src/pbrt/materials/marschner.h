
/*
 Dual hair scattering algorithm for PBRT
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef MARSCHNERMATERIAL_H
#define MARSCHNERMATERIAL_H

// materials/marschner.h*
#include "pbrt.h"
#include "material.h"
#include "reflection.h"

namespace pbrt {

    enum ScatteringMode { R=0, TT=1, TRT=2 };
    
    // PlasticMaterial Declarations
    class MarschnerMaterial : public Material {
      public:
        MarschnerMaterial(const Float alpha[], 
                        const Float beta[], 
                        Float hairRadius,
                        Float eta, 
                        Float eccentricity, 
                        Float glintScaleFactor, 
                        Float causticWidth, 
                        Float causticFade, 
                        Float causticLimit,
                        const std::shared_ptr<Texture<Spectrum>> &sigmaA) 
            : mAr(alpha[0]), mAtt(alpha[1]), mAtrt(alpha[2]),                
            mBr(beta[0]), mBtt(beta[1]), mBtrt(beta[2]),
            mHairRadius(hairRadius), 
            mEta(eta), 
            mEccentricity(eccentricity), 
            mGlintScaleFactor(glintScaleFactor), 
            mCausticWidth(causticWidth), 
            mCausticFadeRange(causticFade), 
            mCausticIntensityLimit(causticLimit),
            mSigmaA(sigmaA)
            {}

        void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                        TransportMode mode,
                                        bool allowMultipleLobes) const;

      private:
        Float mAr, mAtt, mAtrt;
        Float mBr, mBtt, mBtrt;
        Float mHairRadius, mEta, mEccentricity, mGlintScaleFactor, 
        mCausticWidth, mCausticFadeRange, mCausticIntensityLimit;
        std::shared_ptr<Texture<Spectrum>> mSigmaA, mKd;
    };
    
    // MarschnerBSDF Declarations
    class MarschnerBSDF : public BxDF {        
      public:
        MarschnerBSDF(const SurfaceInteraction& si, Float h, 
                Float alphaR, Float alphaTT, Float alphaTRT, 
                Float betaR, Float betaTT, Float betaTRT, 
                Float hairRadius, Float eta, Spectrum sigmaA,
                Float eccentricity, Float glintScaleFactor, Float causticWidth, Float causticFadeRange, Float causticIntensityLimit);
        
        /**
         * (Required) Returns the value of the distribution function for the given pair of directions
         *
         * @param wo outgoing direction
         * @param wi incoming direction
         * @return
         */
        virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
        
        /**
         * Returns the value of the distribution function for the given pair of directions.
         *
         * Used for handling scattering that is described by delta distributions
         * as well as for randomly sampling directions from BxDFs that scatter
         * light along multiple directions.
         *
         * @param wo outgoing direction
         * @param wi incoming direction
         * @param sample
         * @param pdf
         * @param sampledType
         * @return
         */
        virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType = nullptr) const;
        
        /**
         * Returns the probability density function between the incoming and outgoing direction.
         * In other words, how likely is it that incoming radiance is reflected
         * (or refracted) towards the outgoing direction.
         * @param wo
         * @param wi
         * @return 
         */       
        virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
        
        /**
         * 
         * @return 
         */
        virtual std::string ToString() const;
        
        // Not implemented for now
        //virtual Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples) const;
        //virtual Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2) const;

    private:
        Float mH;
        const Normal3f mNs, mNg;
        const Vector3f mDpdu, mDpdv;

        // Marschner params
        const Float mAlphaR, mAlphaTT, mAlphaTRT, mBetaR, mBetaTT, mBetaTRT;
        const Float mHairRadius = 1.0;
        const Float mEta, mEccentricity, mGlintScaleFactor, 
        mCausticWidth, mCausticFadeRange, mCausticIntensityLimit;
        Spectrum mSigmaA;
        

        Float M_r(Float theta_h) const;
        Float M_tt(Float theta_h) const;
        Float M_trt(Float theta_h) const;
        
        Spectrum N_r(Float relativePhi, Float etaPerp, Float etaPar) const;
        Spectrum N_tt(Float relativePhi, Float etaPerp, Float etaPar, Float cosThetaT) const;
        Spectrum N_trt(Float relativePhi, Float etaPerp, Float etaPar, Float cosThetaT) const;
        
        Spectrum N_p(int p, Float relativePhi) const;
    };

    MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // MARSCHNERMATERIAL_H
