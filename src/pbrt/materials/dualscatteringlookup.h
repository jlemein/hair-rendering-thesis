
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_DUALSCATTERING_LOOKUP_H
#define PBRT_MATERIALS_DUALSCATTERING_LOOKUP_H

// materials/dualscattering.h*
#include "pbrt.h"
#include "material.h"
#include "marschner.h"
#include <vector>
#include <string>


namespace pbrt {
    
    class DualScatteringBSDF;
    class OpenVdbReader;
    
    class DualScatteringLookup {
    public:
        /**
         * Gets 
         * @param bsdf
         * @return 
         */
        static const DualScatteringLookup* Get(const DualScatteringBSDF* bsdf);
        
        /**
         * Resets the lookup table, so that next time lookup table is recalculated.
         * Necessary for plotting different hair models, so that each hair model will
         * have it's own lookup data.
         * TODO: not optimal yet
         * @return 
         */
        static void Reset();


        Spectrum AverageForwardScatteringAttenuation(Float thetaD) const;
        Spectrum AverageBackwardScatteringAttenuation(Float thetaD) const;

        Spectrum AverageBackwardScatteringAlpha(Float thetaD) const;
        Spectrum AverageBackwardScatteringBeta(Float thetaD) const;
        Spectrum AverageForwardScatteringAlpha(Float thetaD) const;
        Spectrum AverageForwardScatteringBeta(Float thetaD) const;

        const OpenVdbReader* getVdbReader() const;

    private:
        static DualScatteringLookup* instance;



        DualScatteringLookup(const DualScatteringBSDF* dualScatteringBSDF) : mDualScatteringBSDF(dualScatteringBSDF) {}
        void Init();

        void ReadVoxelGrid();
        void PrecomputeAverageScatteringAttenuation();

        const DualScatteringBSDF* mDualScatteringBSDF;

        // Lookup data
        OpenVdbReader* mVdbReader;
        std::vector<Spectrum> mAverageForwardScatteringAttenuation;
        std::vector<Spectrum> mAverageBackwardScatteringAttenuation;
        std::vector<Spectrum> mAverageForwardScatteringAlpha, mAverageBackwardScatteringAlpha;
        std::vector<Spectrum> mAverageForwardScatteringBeta, mAverageBackwardScatteringBeta;

        const int LOOKUP_TABLE_SIZE = 400;

    };
}

#endif // PBRT_MATERIALS_DUALSCATTERING_LOOKUP_H