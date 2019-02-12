
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
//    std::shared_ptr<Texture<Spectrum>> Kd, Ks;
//    std::shared_ptr<Texture<Float>> roughness, bumpMap;
//    const bool remapRoughness;
      
    // Input parameters to Marschner hair scattering
    Float mAr, mAtt, mAtrt;
    Float mBr, mBtt, mBtrt;
    Float mHairRadius, mEta, mEccentricity, mGlintScaleFactor, mCausticWidth, mCausticFade, mCausticLimit;
    std::shared_ptr<Texture<Spectrum>> mSigmaA, mKd;
};

MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // MARSCHNERMATERIAL_H
