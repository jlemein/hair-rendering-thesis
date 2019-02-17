
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

    // PlasticMaterial Declarations
    class MarschnerMaterial : public Material {
      public:
        MarschnerMaterial(Float Ar, 
                        Float Br, 
                        Float hairRadius,
                        Float eta, 
                        Float eccentricity, 
                        Float glintScaleFactor, 
                        Float causticWidth, 
                        Float causticFade, 
                        Float causticLimit,
                        const std::shared_ptr<Texture<Spectrum>> &sigmaA, 
                        const std::shared_ptr<Texture<Spectrum>> &Kd) 
            : mAr(Ar), mAtt(-.5 * Ar), mAtrt(-1.5*Ar),                
            mBr(Br), mBtt(0.5*Br), mBtrt(2.0*Br),
            mHairRadius(hairRadius), 
            mEta(eta), 
            mEccentricity(eccentricity), 
            mGlintScaleFactor(glintScaleFactor), 
            mCausticWidth(causticWidth), 
            mCausticFade(causticFade), 
            mCausticLimit(causticLimit),
            mSigmaA(sigmaA), 
            mKd(Kd) {}

        void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                        TransportMode mode,
                                        bool allowMultipleLobes) const;

      private:
        Float mAr, mAtt, mAtrt;
        Float mBr, mBtt, mBtrt;
        Float mHairRadius, mEta, mEccentricity, mGlintScaleFactor, mCausticWidth, mCausticFade, mCausticLimit;
        std::shared_ptr<Texture<Spectrum>> mSigmaA, mKd;
    };
    
    // MarschnerBSDF Declarations
    class MarschnerBSDF : public BxDF {
      public:
        MarschnerBSDF(const SurfaceInteraction& si);
        
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
        //virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType = nullptr) const;
        
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
          const Normal3f mNs, mNg;
          const Vector3f mDpdu, mDpdv;
    };

    MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // MARSCHNERMATERIAL_H
