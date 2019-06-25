#include "integrators/mis.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

namespace pbrt {

    MISIntegrator::MISIntegrator(int maxDepth,
            std::shared_ptr<const Camera> camera,
            std::shared_ptr<Sampler> sampler,
            const Bounds2i &pixelBounds, Float rrThreshold,
            const std::string &lightSampleStrategy)
    : SamplerIntegrator(camera, sampler, pixelBounds),
    maxDepth(maxDepth),
    rrThreshold(rrThreshold),
    lightSampleStrategy(lightSampleStrategy) {
    }

    Spectrum MISIntegrator::Li(const RayDifferential &ray, const Scene &scene,
            Sampler &sampler, MemoryArena &arena,
            int depth) const {

        return Spectrum(.5);
    }

    MISIntegrator* CreateMISIntegrator(const ParamSet &params,
            std::shared_ptr<Sampler> sampler,
            std::shared_ptr<const Camera> camera) {
        int maxDepth = params.FindOneInt("maxdepth", 5);

        int np;
        const int *pb = params.FindInt("pixelbounds", &np);
        Bounds2i pixelBounds = camera->film->GetSampleBounds();
        if (pb) {
            if (np != 4)
                Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                    np);
            else {
                pixelBounds = Intersect(pixelBounds,
                        Bounds2i{
                    {pb[0], pb[2]},
                    {pb[1], pb[3]}
                });
                if (pixelBounds.Area() == 0)
                    Error("Degenerate \"pixelbounds\" specified.");
            }
        }
        Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
        std::string lightStrategy =
                params.FindOneString("lightsamplestrategy", "spatial");


        return new MISIntegrator(maxDepth, camera, sampler, pixelBounds,
                rrThreshold, lightStrategy);
    }

}