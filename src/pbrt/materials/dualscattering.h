
/*
 Dual hair scattering algorithm for PBRT
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_DUALSCATTERING_H
#define PBRT_MATERIALS_DUALSCATTERING_H

// materials/dualscattering.h*
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// PlasticMaterial Declarations
class DualscatteringMaterial : public Material {
  public:
    // PlasticMaterial Public Methods
    DualscatteringMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
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

DualscatteringMaterial *CreateDualscatteringMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_DUALSCATTERING_H
