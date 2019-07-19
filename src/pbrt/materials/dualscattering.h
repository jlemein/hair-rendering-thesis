
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
#include "marschner.h"
#include <vector>
#include <string>

namespace pbrt {

    class VisibilityTester;
    class DualScatteringLookup;
    
    class DualscatteringMaterial : public Material {
      public:
        // PlasticMaterial Public Methods
        DualscatteringMaterial(Float eta, MarschnerMaterial* marschnerMaterial, 
                Float alphaR, Float alphaTT, Float alphaTRT, 
                Float betaR, Float betaTT, Float betaTRT,
                Float df, Float db, Float scatterCount,
                std::string voxelGridFileName) 
        : mEta(eta), mMarschnerMaterial(marschnerMaterial), 
                mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT), 
                mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT), 
                mDf(df), mDb(db), mScatterCount(scatterCount), 
                mVoxelGridFileName(voxelGridFileName) {
        }

        void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                        TransportMode mode,
                                        bool allowMultipleLobes) const;

      private:
          MarschnerMaterial* mMarschnerMaterial = nullptr;
          Float mEta;
          Float mAlphaR, mAlphaTT, mAlphaTRT, mBetaR, mBetaTT, mBetaTRT;
          Float mDf, mDb;
          Float mScatterCount;
          std::string mVoxelGridFileName;
    };

    struct GlobalScatteringInformation {
        Spectrum transmittance;
        bool isDirectIlluminated;
        Spectrum variance;
    };

    class DualScatteringBSDF : public BxDF {
    public:
        DualScatteringBSDF(const SurfaceInteraction& si, const Scene* scene,
                Float eta,
                MarschnerBSDF* marschnerBSDF,
                Float alphaR, Float alphaTT, Float alphaTRT,
                Float betaR, Float betaTT, Float betaTRT,
                Float df, Float db, Float scatterCount,
                std::string voxelGridFileName);

        Vector3f WorldToLocal(const Vector3f& v) const;
        Vector3f LocalToWorld(const Vector3f& v) const;

        /**
            * (Required) Returns the value of the distribution function for the given pair of directions
            *
            * @param wo outgoing direction
            * @param wi incoming direction
            * @return
            */
        virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
        virtual Spectrum f(const Vector3f &wo, const Vector3f &wi, const Scene* scene, const VisibilityTester* visibilityTester) const;

        /**
         * @override
         * Here occlusion is handled. It is in the material description, because some materials have approximations to the amount
         * of global scattering coming at a point (as the dual scattering approximation has).
         */
        virtual Spectrum Li(const Vector3f& wi, const Spectrum& Li, const Scene& scene, const VisibilityTester& visibilityTester) const;

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
        virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType = nullptr) const;
        
        /**
         * Returns the probability density function between the incoming and outgoing direction.
         * In other words, how likely is it that incoming radiance is reflected
         * (or refracted) towards the outgoing direction.
         * @param wo
         * @param wi
         * @return 
         */       
        virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
        
        // uniform sampling
        Spectrum UniformSample_f(const Vector3f &wo, Vector3f *wi, Float *pdf, const Point2f& u) const;
        Float UniformPdf(const Vector3f &wo, const Vector3f &wi) const;
        
        // sampling according to DEon et al.
        Spectrum DEonSample_f(const Vector3f &wi, Vector3f *wo, Float *pdf, const Float u[6]) const;
        Float DEonPdf(const Vector3f &wo, const Vector3f &wi) const;

        virtual std::string ToString() const;

            // expensive functions
        Spectrum AverageForwardScatteringAttenuation(Float thetaD) const;
        Spectrum AverageBackwardScatteringAttenuation(Float thetaD) const;

        Spectrum AverageForwardScatteringAlpha(Float thetaD) const;
        Spectrum AverageBackwardScatteringAlpha(Float thetaD) const;
        Spectrum AverageForwardScatteringBeta(Float thetaD) const;
        Spectrum AverageBackwardScatteringBeta(Float thetaD) const;

    private:
        const Scene* mScene = 0;
        MarschnerBSDF* mMarschnerBSDF;
        Float mEta;
        Point3f mPosition;
        const Shape* mShape;
        const Transform mWorldToObject, mObjectToWorld;
        const Bounds3f mObjectBound, mWorldBound;
        Float mDf, mDb;
        Float mScatterCount;
        Float mBetaR, mBetaTT, mBetaTRT;
        Float mAlphaR, mAlphaTT, mAlphaTRT;
        const DualScatteringLookup* mLookup;
        std::string mVoxelGridFileName;
        const Normal3f ns, ng;
        const Vector3f ss, ts;
        const bool mSampleUniform;

        GlobalScatteringInformation GatherGlobalScatteringInformation(const Scene* scene, const VisibilityTester* visiblityTester, const Vector3f& wd, Float thetaD) const;

        Spectrum ForwardScatteringTransmittance(Float n, Float thetaD) const;
        Spectrum ForwardScatteringVariance(Float n, Float thetaD) const;

        Spectrum BackscatteringAttenuation(Float theta) const;
        Spectrum BackscatteringMean(Float theta) const;
        Spectrum BackscatteringVariance(Float theta) const;  

        Spectrum MG_r(Float theta, Spectrum forwardScatteringVariance) const;
        Spectrum MG_tt(Float theta, Spectrum forwardScatteringVariance) const;
        Spectrum MG_trt(Float theta, Spectrum forwardScatteringVariance) const;

        Spectrum EvaluateForwardScatteredMarschner(Float theta_r, Float theta_h,
            Float theta_d, Float phi, Spectrum forwardScatteredVariance) const;

        Spectrum f_R(const Vector3f& wo, const Vector3f& wi) const;
        Spectrum f_TT(const Vector3f& wo, const Vector3f& wi) const;
        Spectrum f_TRT(const Vector3f& wo, const Vector3f& wi) const;

        friend DualScatteringLookup;
    };

    DualscatteringMaterial *CreateDualscatteringMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_DUALSCATTERING_H
