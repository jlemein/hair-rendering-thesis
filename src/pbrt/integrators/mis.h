
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef INTEGRATORS_MIS_H
#define INTEGRATORS_MIS_H

// integrators/mis.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

namespace pbrt {

// MISIntegrator Declarations
class MISIntegrator : public SamplerIntegrator {
public:
    MISIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                   std::shared_ptr<Sampler> sampler,
                   const Bounds2i &pixelBounds, Float rrThreshold = 1,
                   const std::string &lightSampleStrategy = "spatial");
    
    //virtual void Preprocess(const Scene &scene, Sampler &sampler) {}
    //void Render(const Scene &scene);
    virtual Spectrum Li(const RayDifferential &ray, const Scene &scene,
                        Sampler &sampler, MemoryArena &arena,
                        int depth = 0) const;
    
private:
    // MISIntegrator Private Data
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;
    
};

MISIntegrator *CreateMISIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera);

} // namespace pbrt

#endif  // INTEGRATORS_PATH_H
