
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// materials/plastic.cpp*
#include "dualscattering.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "materials/marschner.h"

namespace pbrt {

    // DualscatteringMaterial Method Definitions

    void DualscatteringMaterial::ComputeScatteringFunctions(
            SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
            bool allowMultipleLobes) const {

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);

        MarschnerBSDF* marschnerBSDF = mMarschnerMaterial->CreateMarschnerBSDF(si, arena);

        si->bsdf->Add(ARENA_ALLOC(arena, DualScatteringBSDF)(*si, mEta, marschnerBSDF));
    }

    DualscatteringMaterial *CreateDualscatteringMaterial(const TextureParams &mp) {
        Error("Dualscattering material working!!");

        MarschnerMaterial* marschnerMaterial = CreateMarschnerMaterial(mp);

        Float eta = 1.55;

        // collect dual scattering properties

        return new DualscatteringMaterial(eta, marschnerMaterial);
    }

    /*******************************
     * DualScatteringBSDF
     *******************************/

    DualScatteringBSDF::DualScatteringBSDF(const SurfaceInteraction& si, Float eta,
            MarschnerBSDF* marschnerBSDF)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)), mEta(eta), mMarschnerBSDF(marschnerBSDF) {
    };

    Spectrum DualScatteringBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
        return mMarschnerBSDF->f(wo, wi);
    }

    Spectrum DualScatteringBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const {

        return mMarschnerBSDF->Sample_f(wo, wi, sample, pdf, sampledType);
    }

    Float DualScatteringBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {

        return PiOver4;
    }

    std::string DualScatteringBSDF::ToString() const {
        return "DualScatteringBSDF";
    }

} // namespace pbrt
