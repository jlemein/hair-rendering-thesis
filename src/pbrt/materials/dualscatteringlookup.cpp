#include "dualscatteringlookup.h"
#include "dualscattering.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "materials/marschner.h"
#include "hairutil.h"
#include "hair.h"
#include "light.h"
#include "samplers/random.h"
#include "sampling.h"
#include <random>
#include <mutex>
#include "materials/dualscattering.h"
#include "scene.h"
#include "materials/hairutil.h"
#include "materials/OpenVdbReader.h"

namespace pbrt {

    DualScatteringLookup * DualScatteringLookup::instance = 0;

    std::mutex myMutex;

    const DualScatteringLookup * DualScatteringLookup::Get(const DualScatteringBSDF * bsdf) {

        // We are locking here, because the singleton requires initialization
        std::lock_guard<std::mutex> myLock(myMutex);

        if (DualScatteringLookup::instance == 0) {

            instance = new DualScatteringLookup(bsdf);
            instance->Init();
        }

        return instance;
    }

    void DualScatteringLookup::Init() {

        this->ReadVoxelGrid();
        this->PrecomputeAverageScatteringAttenuation();
    }

    void DualScatteringLookup::ReadVoxelGrid() {
        printf("\rReading VoxelGrid FileName: %s\n", this->mDualScatteringBSDF->mVoxelGridFileName.c_str());
        try {
            mVdbReader = new OpenVdbReader(this->mDualScatteringBSDF->mVoxelGridFileName);
            mVdbReader->initialize();
        } catch (const std::exception& e) {

            printf("\r[ERROR]: Reading VoxelGrid failed: %s\n", e.what());
            throw e;
        }
        printf("[SUCCESS]: Read voxel grid\n");
    }

    const OpenVdbReader * DualScatteringLookup::getVdbReader() const {

        return mVdbReader;
    }

    void DualScatteringLookup::PrecomputeAverageScatteringAttenuation() {
        printf("\rPRECOMPUTING FOR LOOKUP TABLE SIZE: %d\n", LOOKUP_TABLE_SIZE);
        for (int i = 0; i < LOOKUP_TABLE_SIZE; ++i) {

            Float ratio = i / static_cast<Float> (LOOKUP_TABLE_SIZE - 1);
            Float thetaD = (-.5 + ratio) * Pi;

            mAverageBackwardScatteringAttenuation.push_back(mDualScatteringBSDF->AverageBackwardScatteringAttenuation(thetaD));
            mAverageForwardScatteringAttenuation.push_back(mDualScatteringBSDF->AverageForwardScatteringAttenuation(thetaD));

            mAverageForwardScatteringAlpha.push_back(mDualScatteringBSDF->AverageForwardScatteringAlpha(thetaD));
            mAverageBackwardScatteringAlpha.push_back(mDualScatteringBSDF->AverageBackwardScatteringAlpha(thetaD));

            mAverageForwardScatteringBeta.push_back(mDualScatteringBSDF->AverageForwardScatteringBeta(thetaD));
            mAverageBackwardScatteringBeta.push_back(mDualScatteringBSDF->AverageBackwardScatteringBeta(thetaD));
        }

        for (int i = 0; i < LOOKUP_TABLE_SIZE; ++i) {
            //                    PrintSpectrum("AvgBackScatterAttenuation:", mAverageBackwardScatteringAttenuation[i]);
            //                    PrintSpectrum("AvgBackScatterAlpha:", mAverageBackwardScatteringAlpha[i]);
            //                    PrintSpectrum("AvgBackScatterBeta:", mAverageBackwardScatteringBeta[i]);
            Float ratio = i / static_cast<Float> (LOOKUP_TABLE_SIZE - 1);
            Float thetaD = (-.5 + ratio) * Pi;

            printf("thetaD: %f -- forw: %f -- backw: %f\n", thetaD, mAverageForwardScatteringAttenuation[i].y(), mAverageBackwardScatteringAttenuation[i].y());
            //                    PrintSpectrum("AvgForwardScatterAttenuation:", mAverageForwardScatteringAttenuation[i]);
            //                    PrintSpectrum("AvgForwardScatterAlpha:", mAverageForwardScatteringAlpha[i]);
            //                    PrintSpectrum("AvgForwardScatterBeta:", mAverageForwardScatteringBeta[i]);
        }

        printf(" [DONE]\n");
    }

    Spectrum Lookup(Float value, Float min, Float max, const std::vector<Spectrum>& lookupTable) {
        Float range = max - min;
        int lookupIndex = round((lookupTable.size() - 1) * (value - min) / range);

        CHECK_GE(lookupIndex, 0);
        CHECK_LT(lookupIndex, lookupTable.size());

        return lookupTable[lookupIndex];
    }

    Spectrum DualScatteringLookup::AverageBackwardScatteringAttenuation(Float thetaD) const {

        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageBackwardScatteringAttenuation);
    }

    Spectrum DualScatteringLookup::AverageForwardScatteringAttenuation(Float thetaD) const {

        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageForwardScatteringAttenuation);
    }

    Spectrum DualScatteringLookup::AverageBackwardScatteringAlpha(Float thetaD) const {

        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageBackwardScatteringAlpha);
    }

    Spectrum DualScatteringLookup::AverageBackwardScatteringBeta(Float thetaD) const {

        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageBackwardScatteringBeta);
    }

    Spectrum DualScatteringLookup::AverageForwardScatteringAlpha(Float thetaD) const {

        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageForwardScatteringAlpha);
    }

    Spectrum DualScatteringLookup::AverageForwardScatteringBeta(Float thetaD) const {

        return Lookup(thetaD, -.5 * Pi, .5 * Pi, mAverageForwardScatteringBeta);
    }
}