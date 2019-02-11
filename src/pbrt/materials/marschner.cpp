/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   MarschnerMaterial.cpp
 * Author: jeffrey
 *
 * Created on February 11, 2019, 11:34 PM
 */

// materials/plastic.cpp*
#include "dualscattering.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"

namespace pbrt {

    // MarschnerMaterial Method Definitions

    void MarschnerMaterial::ComputeScatteringFunctions(
            SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
            bool allowMultipleLobes) const {
        // Perform bump mapping with _bumpMap_, if present
        if (bumpMap) Bump(bumpMap, si);
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
        // Initialize diffuse component of plastic material
        Spectrum kd = Kd->Evaluate(*si).Clamp();
        if (!kd.IsBlack())
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));

        // Initialize specular component of plastic material
        Spectrum ks = Ks->Evaluate(*si).Clamp();
        if (!ks.IsBlack()) {
            Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
            // Create microfacet distribution _distrib_ for plastic material
            Float rough = roughness->Evaluate(*si);
            if (remapRoughness)
                rough = TrowbridgeReitzDistribution::RoughnessToAlpha(rough);
            MicrofacetDistribution *distrib =
                    ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(rough, rough);
            BxDF *spec =
                    ARENA_ALLOC(arena, MicrofacetReflection)(ks, distrib, fresnel);
            si->bsdf->Add(spec);
        }
    }

    MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp) {
        Error("Marschner material is working!!");
        std::shared_ptr<Texture < Spectrum>> Kd =
                mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
        std::shared_ptr<Texture < Spectrum>> Ks =
                mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
        std::shared_ptr<Texture < Float>> roughness =
                mp.GetFloatTexture("roughness", .1f);
        std::shared_ptr<Texture < Float>> bumpMap =
                mp.GetFloatTextureOrNull("bumpmap");
        bool remapRoughness = mp.FindBool("remaproughness", true);
        return new DualscatteringMaterial(Kd, Ks, roughness, bumpMap, remapRoughness);
    }

} // namespace pbrt
