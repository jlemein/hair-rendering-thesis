
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
    // PlasticMaterial Public Methods
    MarschnerMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
                    const std::shared_ptr<Texture<Spectrum>> &Ks,
                    const std::shared_ptr<Texture<Float>> &roughness,
                    const std::shared_ptr<Texture<Float>> &bumpMap,
                    bool remapRoughness)
        : Kd(Kd),
          Ks(Ks),
          roughness(roughness),
          bumpMap(bumpMap),
          remapRoughness(remapRoughness) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // DualscatteringMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kd, Ks;
    std::shared_ptr<Texture<Float>> roughness, bumpMap;
    const bool remapRoughness;
};

MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // MARSCHNERMATERIAL_H
