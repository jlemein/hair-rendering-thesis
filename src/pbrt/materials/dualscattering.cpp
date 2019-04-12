
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
#include "hairutil.h"
#include "hair.h"

namespace pbrt {

    // DualscatteringMaterial Method Definitions

    void DualscatteringMaterial::ComputeScatteringFunctions(
            SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
            bool allowMultipleLobes) const {

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);

        MarschnerBSDF* marschnerBSDF = mMarschnerMaterial->CreateMarschnerBSDF(si, arena);

        si->bsdf->Add(ARENA_ALLOC(arena, DualScatteringBSDF)(*si, mEta, marschnerBSDF,
                mAlphaR, mAlphaTT, mAlphaTRT, mBetaR, mBetaTT, mBetaTRT));
    }

    DualscatteringMaterial *CreateDualscatteringMaterial(const TextureParams &mp) {
        MarschnerMaterial* marschnerMaterial = CreateMarschnerMaterial(mp);

        // TODO: Check if dual scattering really wants to square the widths
        marschnerMaterial->mBr = Sqr(marschnerMaterial->mBr);
        marschnerMaterial->mBtt = Sqr(marschnerMaterial->mBtrt);
        marschnerMaterial->mBtrt = Sqr(marschnerMaterial->mBtrt);

        Float eta = 1.55;

        // collect dual scattering properties

        return new DualscatteringMaterial(eta, marschnerMaterial,
                marschnerMaterial->mAr, marschnerMaterial->mAtt, marschnerMaterial->mAtrt,
                marschnerMaterial->mBr, marschnerMaterial->mBtt, marschnerMaterial->mBtrt);
    }

    static Float Transmittance() {

    }

    static Float Variance() {

    }

    /*******************************
     * DualScatteringBSDF
     *******************************/

    DualScatteringBSDF::DualScatteringBSDF(const SurfaceInteraction& si, Float eta,
            MarschnerBSDF* marschnerBSDF,
            Float alphaR, Float alphaTT, Float alphaTRT,
            Float betaR, Float betaTT, Float betaTRT)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
    mEta(eta),
    mMarschnerBSDF(marschnerBSDF),
    mPosition(si.p),
    mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT),
    mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT) {
    };

    /**
     *
     * @param eccentricity Eccentricity value (1 means circular)
     * @param etaPerp The index of refraction of the fiber
     * @param thetaH The longitudinal half angle
     * @return the adjusted index of refraction taking into account eccentricity
     */
    static Float EtaEccentricity(Float eccentricity, Float etaPerp, Float thetaH) {
        Float eta1 = 2.0 * (etaPerp - 1.0) * Sqr(eccentricity) - etaPerp + 2.0;
        Float eta2 = 2.0 * (etaPerp - 1.0) / Sqr(eccentricity) - etaPerp + 2.0;

        return 0.5 * ((eta1 + eta2) + cos(2.0 * thetaH)*(eta1 - eta2));
    }

    Spectrum DualScatteringBSDF::f(const Vector3f &wo, const Vector3f &wi) const {

        Float thetaI, phiI, thetaR, phiR;
        ToSphericalCoords(wi, thetaI, phiI);
        ToSphericalCoords(wo, thetaR, phiR);

        Float thetaD = DifferenceAngle(thetaI, thetaR);
        Float phi = RelativeAzimuth(phiI, phiR);
        Float thetaH = HalfAngle(thetaI, thetaR);
        Float phiH = HalfAngle(phiI, phiR);

        // Global multiple scattering
        // TODO: Check if wd = wi?
        const Vector3f& wd = wi;
        static GlobalScatteringInformation gsi;
        GatherGlobalScatteringInformation(wd, gsi);

        // Shading

        // compute local multiple scattering contribution
        // TODO: Is this correct?
        Float theta = thetaH; //0.0; Is this correct?

        Spectrum Ab = BackscatteringAttenuation(theta);
        Float deltaB = BackscatteringMean(theta);

        Float backScatterVariance = BackscatteringVariance(theta);
        Float forwardScatterVariance = ForwardScatteringVariance(theta);
        Spectrum fBack = Spectrum(.0); //2.0 * Ab * Gaussian(backScatterVariance + forwardScatterVariance,
        //thetaD + thetaR - deltaB) / (Pi * Sqr(cos(theta)));

        // compute BCSDF of the fiber due to direct illumination
        //TODO: find a way to use beta squared
        Spectrum fDirectS = mMarschnerBSDF->f(wo, wi);
        Spectrum FDirect = gsi.directIlluminationFraction * (fDirectS + mDb * fBack);

        // Compute BCSDF of the fiber due to forward scattered illumination similarly
        Spectrum fScatterS = EvaluateForwardScatteredMarschner(thetaR, thetaH, thetaD, phi, forwardScatterVariance);
        //TODO: Does this Pi belong here, or is it a typo in paper??
        Spectrum FScatter = (gsi.transmittance - gsi.directIlluminationFraction) * mDf * (fScatterS + Pi * mDb * fBack);

        // combine direct and forward scattered components
        return (FDirect + 0.0 * FScatter) * cos(thetaI);
    }

    void DualScatteringBSDF::GatherGlobalScatteringInformation(const Vector3f& wd, GlobalScatteringInformation& gsi) const {
        //if (this->FindScatteringCount(wd) == 0) {
        gsi.directIlluminationFraction = 1.0;
        gsi.transmittance = 1.0;
        gsi.variance = 0.0;
        //        } else {
        //            gsi.directIlluminationFraction = 0.0;
        //            gsi.transmittance = ForwardTransmittance();
        //            gsi.variance = Variance();
        //        }
    }

    int DualScatteringBSDF::FindScatteringCount(const Vector3f& wd) const {
        //return VoxelGrid.interpolateP(mPosition, wd);
        return 2;
    }

    Spectrum DualScatteringBSDF::EvaluateForwardScatteredMarschner(Float thetaR, Float thetaH, Float thetaD,
            Float phi, Float forwardScatteredVariance) const {

        // TODO: check if we can remove eta perp calculation and let marschner
        // be responsible for it
        Float etaPerp, etaPar;
        ToBravais(mEta, thetaR, etaPerp, etaPar);

        // take into account eccentricity (by again adjusting the eta)
        Float eccentricity = mMarschnerBSDF->getEccentricity();
        if (eccentricity != 1.0) {
            etaPerp = EtaEccentricity(eccentricity, etaPerp, thetaH);
        }

        // Very similar to Marschner, but MG_r, MG_tt and MG_trt are used
        // to incorporate the forward scattering variance
        Float sinThetaR = sin(thetaR);
        Float sinThetaT = sinThetaR / etaPerp;
        Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        Spectrum result = (
                MG_r(2.0 * thetaH, forwardScatteredVariance) * mMarschnerBSDF->N_r(phi, etaPerp)
                + MG_tt(2.0 * thetaH, forwardScatteredVariance) * mMarschnerBSDF->N_tt(phi, etaPerp, cosThetaT)
                + MG_trt(2.0 * thetaH, forwardScatteredVariance) * mMarschnerBSDF->N_trt(phi, etaPerp, cosThetaT)
                ) / CosineSquared(thetaD);

        return result;
    }

    Spectrum DualScatteringBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const {

        return mMarschnerBSDF->Sample_f(wo, wi, sample, pdf, sampledType);
    }

    Float DualScatteringBSDF::MG_r(Float theta, Float forwardScatteringVariance) const {
        // WATCHOUT: do not square betaR, already squared at material creation
        return Gaussian(mBetaR + forwardScatteringVariance, theta - mAlphaR);
    }

    Float DualScatteringBSDF::MG_tt(Float theta, Float forwardScatteringVariance) const {
        // WATCHOUT: do not square betaTT, already squared at material creation
        return Gaussian(mBetaTT + forwardScatteringVariance, theta - mAlphaTT);
    }

    Float DualScatteringBSDF::MG_trt(Float theta, Float forwardScatteringVariance) const {
        // WATCHOUT: do not square betaTRT, already squared at material creation
        return Gaussian(mBetaTRT + forwardScatteringVariance, theta - mAlphaTRT);
    }

    /**
     * Equation 5 in dual scattering approximation
     */
    Spectrum DualScatteringBSDF::ForwardScatteringTransmittance() const {
        //        Float af = mMarschnerBSDF->f(wo, wi) * cos(thetaD);
        //        return df * pow(af)
        return Spectrum(.0);
    }

    Float DualScatteringBSDF::ForwardScatteringVariance(Float theta) const {
        return 1.0;
    }

    Spectrum DualScatteringBSDF::BackscatteringAttenuation(Float theta) const {
        //        Float response = 0.0;
        //
        //        for (Float phi = -Pi; phi < Pi; phi += 0.001) {
        //            for (Float theta = -.5 * Pi; theta < .5 * Pi; theta += 0.001) {
        //                Vector3f wi;
        //                FromSphericalCoords(theta, phi, wi);
        //
        //                response += mMarschnerBSDF->f(wd, wi);
        //            }
        //        }
        //        return 0.5 * Pi * response * cos(thetaD);
        return Spectrum(.0);
    }

    Float DualScatteringBSDF::BackscatteringMean(Float theta) const {
        return 1.0;
    }

    Float DualScatteringBSDF::BackscatteringVariance(Float theta) const {
        return 1.0;
    }

    Float DualScatteringBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {

        return PiOver4;
    }

    std::string DualScatteringBSDF::ToString() const {
        return "DualScatteringBSDF";
    }

} // namespace pbrt
