
/*
 * File:   MarschnerMaterial.cpp
 * Author: jeffrey lemein
 *
 * Created on February 11, 2019, 11:34 PM
 */

#include "marschner.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "hair.h"

namespace pbrt {

    /**
     * Safe acos function, to prevent 'nan' as result when parameter x is
     * (slightly) larger than or smaller than 1.0 or -1.0 respectively
     * (due to rounding errors).
     * @param x Argument to the acosine function. Values smaller than
     * -1.0 and larger than 1.0 are mapped to -1.0 and 1.0 respectively.
     * @return The inverse cosine result for parameter x.
     */
    static inline Float SafeACos(Float x) {
        if (x <= -1.0) {
            return Pi;
        }
        if (x >= 1.0) {
            return 0.0;
        }
        return acos(x);
    }

    static inline Float SafeZero(Float x) {
        return std::max(Float(1e-5), x);
    }

    static void ToSphericalCoords(const Vector3f& w, Float& theta, Float& phi) {
        theta = PiOver2 - SafeACos(w.x);
        phi = atan2(w.y, w.z);

        CHECK_LE(abs(theta), PiOver2);
        CHECK_LE(abs(phi), Pi);
    }

    static Float DifferenceAngle(Float theta_i, Float theta_r) {
        return 0.5 * (theta_r - theta_i);
    }

    static Float RelativeAzimuth(Float phi_i, Float phi_r) {
        // TODO: check if needs to be wrapped around [0, 2pi]
        //mod(abs(phi_r - phi_i), 2*PI);

        // I think phi should be mapped between [-pi, +pi] to make
        // finding roots work that is accurate for gamma in range [-pi/2, pi/2]

        //return phi_r - phi_i;

        CHECK_LE(abs(phi_i), Pi);
        CHECK_LE(abs(phi_r), Pi);

        if (phi_i > Pi) {
            phi_i -= 2.0 * Pi;
        }
        if (phi_i < -Pi) {
            phi_i += 2.0 * Pi;
        }
        if (phi_r > Pi) {
            phi_r -= 2.0 * Pi;
        }
        if (phi_r < -Pi) {
            phi_r += 2.0 * Pi;
        }

        Float phi = phi_r - phi_i;

        if (phi > Pi) {
            phi -= 2.0 * Pi;
        }
        if (phi < -Pi) {
            phi += 2.0 * Pi;
        }
        if (abs(phi) > Pi) {
            printf("phi_i = %f, phi_r = %f, phi = %f\n", phi_i, phi_r, phi);
        }
        CHECK_LE(abs(phi), Pi);
        return phi;
    }

    static Float HalfAngle(Float a, Float b) {
        return 0.5 * (a + b);
    }

    static inline Float CosineSquared(Float x) {
        return Sqr(cos(x));
    }

    static inline Float SineSquared(Float x) {
        return Sqr(sin(x));
    }

    /**
     * Returns the value of a normalized Gaussian function at point 'x', with a standard deviation of 'width'
     * @param width The width of the gaussian function
     * @param x Position to sample the gaussian
     * @returns Sampled position of the gaussian at position x
     */
    static Float Gaussian(Float width, Float x) {
        Float a = 1.0 / (width * sqrt(2.0 * Pi));
        float c = width; //width of the curve is beta (might also be 0.5 * sigma)

        float nom = Sqr(x);
        float den = 2.0 * Sqr(width);

        return a * exp(-nom / den);
    }

    /**
     * Slightly faster variant when you need to get both Bravais indices
     * @param eta Index of refraction
     * @param gamma Angle of incidence (in radians) or sine angle_of_incidence ??
     * @param bravaisPerpendicular Output value that will hold perpendicular component of bravais index
     * @param bravaisParallel Output value that will hold parallel component of bravais index
     */
    static void ToBravais(Float eta, Float gamma, Float& etaPerpendicular, Float& etaParallel) {
        Float rootPart = sqrt(Sqr(eta) - SineSquared(gamma));

        //TODO: sin gamma = h, where 1 < h < 1
        // Gamma represents the angle between the incident ray and the normal of a dielectric cylinder
        // This means -pi < gamma < pi
        // This also means cos(gamma) is 0 < cos(gamma) <= 1

        Float cosGamma = cos(gamma);

        if (cosGamma < 0.0) {
            printf("cosgamma: %f -- gamma is: %f\n", cosGamma, gamma);
            CHECK_GE(cosGamma, 0.0);
        }

        etaPerpendicular = rootPart / cosGamma;
        etaParallel = Sqr(eta) * cosGamma / rootPart;
    }

    // TODO: gamma is angle of incidence or sine of angle of incidence??

    static Float BravaisPerpendicular(Float eta, Float gamma) {
        return sqrt(Sqr(eta) - Sqr(sin(gamma))) / cos(gamma);
    }

    // TODO: gamma is angle of incidence or sine of angle of incidence??

    static Float BravaisParallel(Float eta, Float gamma) {
        return Sqr(eta) * cos(gamma) / sqrt(Sqr(eta) - Sqr(sin(gamma)));
    }

    static Float FresnelSplitted(Float ni, Float nt, Float gamma_i, Float gamma_t) {
        Float cos_gamma_i = cos(gamma_i);
        Float cos_gamma_t = cos(gamma_t);

        // for p-polarized light
        Float rParallel = (nt * cos_gamma_i - ni * cos_gamma_t) / (nt * cos_gamma_i + ni * cos_gamma_t);

        // for s-polarized light
        Float rPerpendicular = (ni * cos_gamma_i - nt * cos_gamma_t) / (ni * cos_gamma_i + nt * cos_gamma_t);

        return 0.5 * (rParallel + rPerpendicular);
    }

    /**
     * Fresnel reflection for s-polarized light (perpendicular to incidence plane)
     * @param ni Index of refraction for material of incident side
     * @param nt Index of refraction for material on transmitted side
     * @param gamma_i Angle of incidence
     * @param gamma_t Angle of refraction
     * @return
     */
    static Float FresnelReflectionS(Float ni, Float nt, Float gamma_i) {
        Float sinGammaByN = ni * sin(gamma_i) / nt;

        if (sinGammaByN >= 1.0) {
            return 1.0;
        } else {
            Float gamma_t = SafeASin(sinGammaByN);
            Float cos_gamma_i = cos(gamma_i);
            Float cos_gamma_t = cos(gamma_t);
            return (ni * cos_gamma_i - nt * cos_gamma_t) / (ni * cos_gamma_i + nt * cos_gamma_t);
        }
    }

    /**
     * Fresnel reflection for p-polarized light (parallel to incidence plane)
     * @param ni Index of refraction for material of incident side
     * @param nt Index of refraction for material on transmitted side
     * @param gamma_i Angle of incidence
     * @param gamma_t Angle of refraction
     * @return
     */
    static Float FresnelReflectionP(Float ni, Float nt, Float gamma_i) {
        Float sinGammaByN = ni * sin(gamma_i) / nt;
        if (sinGammaByN >= 1.0) {
            return 1.0;
        } else {
            Float gamma_t = SafeASin(sinGammaByN);
            Float cos_gamma_i = cos(gamma_i);
            Float cos_gamma_t = cos(gamma_t);
            return (nt * cos_gamma_i - ni * cos_gamma_t) / (nt * cos_gamma_i + ni * cos_gamma_t);
        }
    }

    /**
     * Fresnel reflection for polarized light
     * @param etaPerp
     * @param etaPar
     * @param gamma_i
     * @return
     */
    static Float Fresnel(Float etaPerp, Float etaPar, Float gamma_i) {
        CHECK_GT(etaPerp, 0.0);
        CHECK_GT(etaPar, 0.0);

        // TODO: check if gamma_i should be projected for parallel and perpendicular planes
        // Marschner indicates that by using the altered index of refractions, that you can
        // just use gamma as angle (appendix B)

        // TODO: must we use gamma_t based on eta, or calculate gamma_t for projected and perpendicular cases
        // use Snell's law to find transmitted angle
        //Float gamma_t = SafeASin(sin(gamma_i) / eta);

        // Something goes wrong here, because sin(gamma) divided
        Float x = sin(gamma_i) / etaPerp;
        Float y = sin(gamma_i) / etaPar;
        //if (x < -1.0001 || x > 1.0001) {
        //printf("safe asin failed, because y = %f ;; etaPar = %f\n", y, etaPar);

        // }
        CHECK(x >= -1.0001 && x <= 1.0001);
        //        Float gamma_t_perpendicular = SafeASin(sin(gamma_i) / etaPerp);
        Float fresnelS = FresnelReflectionS(1.0, etaPerp, gamma_i);
        CHECK_LE(fresnelS, 1.0);

        // for s-polarized light
        //        Float gamma_t_parallel = SafeASin(sin(gamma_i) / etaPar);
        Float fresnelP = FresnelReflectionP(1.0, etaPar, gamma_i);
        CHECK_LE(fresnelP, 1.0);

        return 0.5 * (fresnelP + fresnelS);
    }

    /**
     * Solves the root for reflection (R) mode of Marschner.
     * Gamma is returned instead of the root h (for efficiency reasons)
     *
     *  h = sin gamma, so if you want h, then do arcsin(gamma)
     *
     * @param phi
     * @return Gamma_i
     */
    static inline Float SolveGammaRoot_R(Float phi) {
        return -phi / 2.0;
    }

    static Spectrum Transmittance(const Spectrum& sigmaA, Float gamma_t) {
        return Exp(-2.0 * sigmaA * (1.0 + cos(2.0 * gamma_t)));
    }

    static Float DPhiDh_R(Float gamma_i) {
        return 1.0 / SafeZero(sqrt(1.0 - SineSquared(gamma_i)));
    }

    /*******************************
     * MarschnerMaterial
     *******************************/

    void MarschnerMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
            MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const {

        // Allocate a bsdf that contains the collection of BRDFs and BTDFs
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, this->mEta);
        si->bsdf->Add(ARENA_ALLOC(arena, MarschnerBSDF)(*si, mAr, mAtt, mAtrt, mBr, mBtt, mBtrt, mEta, mSigmaA));
    }

    MarschnerMaterial *CreateMarschnerMaterial(const TextureParams &mp) {
        Float Ar = mp.FindFloat("Ar", Radians(-7.5));
        Float Br = mp.FindFloat("Br", Radians(7.5));
        Float hairRadius = 1.0;
        Float eta = 1.55;
        Float eccentricity = mp.FindFloat("eccentricity", Float(1.0));
        Float glintScaleFactor = 2.5;
        Float causticWidth = Radians(15.0);
        Float causticFade = 0.3;
        Float causticLimit = 0.5;

        Float rgb[3] = {0.4, 0.2, 0.1};
        std::shared_ptr<Texture < Spectrum>> sigmaA = mp.GetSpectrumTexture("sigmaA", Spectrum::FromRGB(rgb)); //should be defined as color
        std::shared_ptr<Texture < Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));

        return new MarschnerMaterial(Ar, Br, hairRadius, eta, eccentricity, glintScaleFactor, causticWidth, causticFade, causticLimit, sigmaA, Kd);
    }

    /*******************************
     * MarschnerBSDF
     *******************************/

    MarschnerBSDF::MarschnerBSDF(const SurfaceInteraction& si,
            Float alphaR, Float alphaTT, Float alphaTRT,
            Float betaR, Float betaTT, Float betaTRT,
            Float eta, std::shared_ptr<Texture < Spectrum>> sigmaA
            )
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
    mNs(si.shading.n), mNg(si.n), mDpdu(si.dpdu), mDpdv(si.dpdv),
    mAlphaR(alphaR), mAlphaTT(alphaTT), mAlphaTRT(alphaTRT),
    mBetaR(betaR), mBetaTT(betaTT), mBetaTRT(betaTRT),
    mEta(eta), mSigmaA(sigmaA) {
    };

    Float MarschnerBSDF::M_r(Float theta_h) const {
        return Gaussian(mBetaR, theta_h - mAlphaR);
    }

    Float MarschnerBSDF::M_tt(Float theta_h) const {
        return Gaussian(mBetaTT, theta_h - mAlphaTT);
    }

    Float MarschnerBSDF::M_trt(Float theta_h) const {
        return Gaussian(mBetaTRT, theta_h - mAlphaTRT);
    }

    Spectrum MarschnerBSDF::N_r(Float relativePhi) const {
        // Reflection has only 1 root
        // TODO: let relativePhi be in specific bounds.
        // Marschner says that rootsolving works for relativePhi between [-pi, +pi]
        //
        Float gamma_i = SolveGammaRoot_R(relativePhi);

        Float etaPerp, etaPar;
        ToBravais(mEta, gamma_i, etaPerp, etaPar);

        // reflection is only determined by Fresnel
        return Fresnel(etaPerp, etaPar, gamma_i) / (2.0 * DPhiDh_R(gamma_i));
    }

    Spectrum MarschnerBSDF::N_tt(Float relativePhi) const {
        return N_p(1, relativePhi);
    }

    Spectrum MarschnerBSDF::N_trt(Float relativePhi) const {
        return N_p(2, relativePhi);
    }

    Spectrum MarschnerBSDF::N_p(int p, Float relativePhi) const {
        return Spectrum(0.01);
    }

    Spectrum MarschnerBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
        // x axis goes with the fiber, from root to tip
        // y represents normal to hair fiber (major axis)
        // z axis represents normal to hair fiber (minor axis)

        Float theta_i, phi_i, theta_r, phi_r;
        ToSphericalCoords(wi, theta_i, phi_i);
        ToSphericalCoords(wo, theta_r, phi_r);

        Float theta_d = DifferenceAngle(theta_i, theta_r);
        Float phi = RelativeAzimuth(phi_i, phi_r);
        Float theta_h = HalfAngle(theta_i, theta_r);
        Float phi_h = HalfAngle(phi_i, phi_r);



        Spectrum result = (
                M_r(theta_h) * N_r(phi)
                + M_tt(theta_h) * N_tt(phi)
                + M_trt(theta_h) * N_trt(phi)) / CosineSquared(theta_d);

        return result;
    }

    //    Spectrum MarschnerBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const {
    //
    //    }

    Float MarschnerBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
        return PiOver4;
    }

    std::string MarschnerBSDF::ToString() const {
        return "MarschnerBSDF";
    }



    //    class marschner
    //(
    //	//Alpha: shift of cuticle scales
    //	uniform float Ar = radians(-7.5); //, Att, Atrt; //alpha values (scattering deviation to account for tilted cuticle scales)
    //	uniform float Att = -.5 * Ar;
    //	uniform float Atrt = -1.5 * Ar;
    //
    //	//Beta: width of the unit-integral zero-mean lobe function
    //	uniform float Br = radians(7.5);
    //	uniform float Btt = Br / 2.0;
    //	uniform float Btrt = 2.0 * Br;
    //
    //	float hairRadius = 1.0;
    //	float eta = 1.55;			//Index of refraction
    //	color sigma_a = color(0.2, 0.3, 0.8);	//Absorption coefficient
    //	float eccentricity = 1.0; 		//hair eccentricity (1.0 for circular cross section)
    //
    //	//Parameters used for TRT calculation
    //	float glintScaleFactor = 2.5;	 // Glint scale factor
    //	float causticWidth = radians(15);// Azimuthal width of caustic
    //	float causticFade = 0.3;	 // Fade range for caustic merge
    //	float causticLimit = 0.5;		// Caustic intensity limit
    //
    //	//Additional diffuse component (suggested by Marschner to improve rendering result)
    //	float Kd = 0.0;
    //
    //	//for debugging purposes to block out or intensify the R, TT or TRT components resp.
    //	color enableComponents = color(1,1,1);
    //)
    //{
    //	/**
    //	 * Returns the cubic root of a value v
    //	 */
    //	public float cubicRoot(float v)
    //	{
    //		float sign = -1.0;
    //		if( v >= 0.0 )
    //			sign = 1.0;
    //
    //		return sign * pow(abs(v), 1.0/3.0);
    //	}
    //
    //	/**
    //	 * First derivative of the equation from Marschner's paper
    //	 */
    //	float dPhi(float p, eta, h; )
    //	{
    //		float gamma = asin(h);
    //		float pi3 = PI * PI * PI;
    //
    //		float c = asin(1.0/eta);
    //		float dGamma_i = (6.0*p*c/PI - 2.0) - 24.0 * (p*c/pi3) * gamma*gamma;
    //
    //		return dGamma_i / max(1e-5, sqrt(1.0 - h*h));
    //	}
    //
    //	float d2Phi( uniform float p; varying float eta; varying float h )
    //	{
    //	        // computes the second derivative of the polynomial relative to h.
    //	        float gamma = asin( h );
    //	        uniform float pi3 = PI*PI*PI;
    //	        float c = asin(1/eta);
    //	        float dGamma = -2*3*8*(p*c/pi3)*gamma;
    //	        float denom = pow( 1-h*h, 3/2 );
    //	        return (dGamma*h)/max(1e-5,denom);
    //	}
    //
    //
    //
    //	/**
    //	 * Returns the value of a normalized Gaussian function at point 'x', with a standard deviation of 'width'
    //	 */
    //	float Gaussian( float width, x; )
    //	{
    //		float a = 1.0 / (width * sqrt(2.0 * PI)); //unit integral gaussian function
    ////		float b = 0.0; //zero mean function
    //		float c = width; //width of the curve is beta (might also be 0.5 * sigma)
    //
    ////		float nom = (x - b) * (x - b);
    //		float nom = x*x;
    //		float den = 2.0 * c * c;
    //
    //		return a * exp( - nom / den );
    //	}
    //
    //	/**
    //	 * Returns the perpendicular component of the Bravais index
    //	 */
    //	float EtaPerp(float eta, theta)
    //	{
    //		float sin_theta = sin(theta);
    //
    //		return sqrt(eta*eta - sin_theta*sin_theta) / cos(theta);
    //	}
    //
    //	/**
    //	 * Returns the parallel component of the Bravais index
    //	 */
    //	float EtaParr(float eta, gamma)
    //	{
    //		float sin_gamma = sin(gamma);
    //
    //		return eta*eta*cos(gamma) / sqrt(eta*eta - sin_gamma*sin_gamma);
    //	}
    //
    //	/**
    //	 *	eta:   index of refraction (usually 1.55 for hair fibers)
    //	 *	a:     eccentricity of cross section (1.0 for circular cross section)
    //	 *  phi_h:
    //	 */
    //	float EtaAsterisk(float a; float phi_h)
    //	{
    //		//TO DO: check if implementation is correct
    //		float eta1 = 2.0 * (eta - 1) * a*a - eta + 2.0;
    //		float eta2 = 2.0 * (eta - 1) / (a*a) - eta + 2.0;
    //
    //		return 0.5 * ((eta1 + eta2) + cos(2.0 * phi_h)*(eta1 - eta2));
    //	}
    //
    //	/**
    //	 * Returns the Fresnel coefficient for the given angle
    //	 * eta_perp: The perpendicular index of refraction
    //	 * eta_parr: The parallel index of refraction
    //	 * incidenceAngle: incidence angle in radians [-pi, pi]
    //	 */
    //	float Fresnel(float eta_perp, eta_parr, incidenceAngle; )
    //	{
    //		float ni, nt;
    //		float r_parr = 1.0;
    //		float r_perp = 1.0;
    //
    //		float angle = abs(incidenceAngle);
    //		if( angle > 0.5*PI )
    //			angle = PI - angle; //[0, pi]
    //
    //		//parallel light reflectance
    //		ni = 1.0;
    //		nt = eta_parr;
    //		float s = (ni/nt)*sin(angle);
    //		s *= s;
    //
    //		if( s <= 1.0 )
    //		{
    //			float a = ni * sqrt(1.0 - s);
    //			float b = nt * cos(angle);
    //			r_parr = (a - b) / (a + b);
    //			r_parr = min(1.0, r_parr*r_parr);
    //		}
    //
    //		//perpendicular light reflectance
    //		ni = 1.0;
    //		nt = eta_perp;
    //
    //		s = (ni/nt)*sin(angle);
    //		s *= s;
    //
    //		if( s <= 1 )
    //		{
    //			float a = ni * cos(angle);
    //			float b = nt * sqrt(1 - s);
    //			r_perp = (a - b) / (a + b);
    //			r_perp = min(1, r_perp*r_perp);
    //		}
    //
    //		//gamma is angle between the incident vector and its projection on the normal plane
    //		return 0.5 * (r_perp + r_parr);
    //	}
    //
    //
    //	float LinearSolver(float a; float b; output float root )
    //	{
    //		float rootCount = -1;
    //        	if (a != 0)
    //        	{
    //                	root = -b / a;
    //                	rootCount = 1;
    //        	}
    //        	else if (b != 0)
    //        	{
    //               		rootCount = 0;
    //        	}
    //        	return rootCount;
    //	}
    //
    //	/*
    //	 * Solves a quadratic equation: ax^2 + bx + c = 0
    //	 * The roots of this equation are stored in 'roots'
    //	 * Returns the number of roots found
    //	 */
    //	float QuadraticSolver(float a, b, c; output float roots[3] )
    //	{
    //		uniform float epsilon = 1e-16;
    //		float rootCount = 0;
    //
    //		if( abs(a) < epsilon )
    //		{
    //			//roots[0] = -c / b;
    //			rootCount = LinearSolver( b, c, roots[0] );
    //		}
    //		else
    //		{
    //
    //			float D = b*b - 4.0*a*c;
    //
    //			if( abs(D) < epsilon )
    //			{
    //				roots[0] = -b / (2.0*a);
    //				rootCount = 1;
    //			}
    //			else if( D > 0 )
    //			{
    //				float s = sqrt(D);
    //				roots[0] = (-b + s) / (2 * a);
    //				roots[1] = (-b - s) / (2 * a);
    //				rootCount = 2;
    //			}
    //		}
    //
    //		return rootCount;
    //	}
    //
    //	float NormalizedCubicSolver( float A; float B; float C; output float roots[] )
    //	{
    //		uniform float epsilon = 1e-16;
    //		float rootCount = 0;
    //
    //		if( abs(C) < epsilon )
    //		{
    //			rootCount = QuadraticSolver(1, A, B, roots);
    //			roots[rootCount] = 0;
    //			rootCount = rootCount + 1;
    //		}
    //		else
    //		{
    //			float Q = (3*B - A*A)/9;
    //			float R = (9*A*B - 27*C - 2*A*A*A)/54;
    //			float D = Q*Q*Q + R*R;
    //
    //			if( D > 0 )
    //			{
    //				float sqrtD = sqrt(D);
    //				float S = cubicRoot(R + sqrtD);
    //				float T = cubicRoot(R - sqrtD);
    //				roots[0] = (-A/3 + (S+T));
    //				rootCount = 1;
    //			}
    //			else
    //			{
    //				float th = acos( R/sqrt(-(Q*Q*Q)) );
    //				float sqrtQ = sqrt(-Q);
    //				roots[0] = (2*sqrtQ*cos(th/3) - A/3);
    //				roots[1] = (2*sqrtQ*cos((th + 2*PI)/3) - A/3);
    //				roots[2] = (2*sqrtQ*cos((th + 4*PI)/3) - A/3);
    //				rootCount = 3;
    //			}
    //		}
    //		return rootCount;
    //	}
    //
    //	/**
    //	 * Solves the following cubic equation: ax^3 + bx^2 + cx + d = 0
    //	 * Roots are stored in array roots[3] and number of roots is returned.
    //	 */
    //	float CubicSolver(float a, b, c, d; output float roots[3])
    //	{
    //		float epsilon = 1e-16;
    //		float rootCount;
    //
    //		if( abs(a) < epsilon )
    //		{
    //			rootCount = QuadraticSolver(b, c, d, roots);
    //		}
    //		else
    //		{
    //			float A = b/a;
    //			float B = c/a;
    //			float C = d/a;
    //			rootCount = NormalizedCubicSolver(A, B, C, roots);
    //		}
    //		return rootCount;
    //		//float A = b/a;
    //		//float B = c/a;
    //		//float C = d/a;
    //
    //		//float Q = (3.0*B - A*A)/9.0;
    //		//float R = (9.0*A*B - 27.0*C - 2.0*A*A*A)/54.0;
    //		//float D = Q*Q*Q + R*R;
    //
    //		//if( D > 0.0 )
    //		//{
    //		//	float sqrtD = sqrt(D);
    //		//	float S = cubicRoot(R + sqrtD);
    //		//	float T = cubicRoot(R - sqrtD);
    //		//	roots[0] = (-A/3.0 + (S+T));
    //		//	return 1;
    //		//}
    //		//else
    //		//{
    //		//	float th = acos( R/sqrt(-(Q*Q*Q)) );
    //		//	float sqrtQ = sqrt(-Q);
    //		//	roots[0] = (2.0*sqrtQ*cos(th/3) - A/3);
    //		//	roots[1] = (2.0*sqrtQ*cos((th + 2*PI)/3) - A/3);
    //		//	roots[2] = (2.0*sqrtQ*cos((th + 4*PI)/3) - A/3);
    //		//	return 3;
    //		//}
    //	}
    //
    //	/**
    //	 * Finds the root for a given reflection component (p), for a particular IOR and angle.
    //	 * Roots are stored in array roots[3] and number of roots is returned.
    //	 */
    //	float h(float p, eta_bravais, phi; output float roots[3])
    //	{
    //		float pi3 = PI * PI * PI;
    //		float c = asin(1.0 / eta_bravais);
    //
    //		//Quadratic terms
    //		float c3 = -8.0 * (p*c/pi3);	//gamma^3
    //		float c2 = 0.0;					//gamma^2
    //		float c1 = 6.0*p*c/PI - 2.0;	//gamma^1
    //		float c0 = p*PI - phi;				//gamma^0
    //
    //		float nRoots = CubicSolver(c3, c2, c1, c0, roots);
    //
    //		return nRoots;
    //	}
    //
    //	/**
    //	 * Absorption function returning the amount of light after the light leaves the hair fiber
    //	 * For the R-case, the Fresnel coefficient will be returned
    //	 * For the TT and TRT-cases the reflected value will also incorporate the absorption inside the hair fiber (this gives the color to TT and TRT)
    //	 */
    //	color A2(float p, gamma_i, eta_perp, eta_parr, theta_i;)
    //	{
    //		float fresnel = Fresnel(eta_perp, eta_parr, gamma_i);
    //
    //		//For the R component (p = 0), return the Fresnel coefficient
    //		if( p == 0 )
    //		{
    //			return color(fresnel, fresnel, fresnel);
    //		}
    //
    //		float h = sin( gamma_i );
    //		float gamma_t = asin( h/eta_perp );
    //		float theta_t = -sign(theta_i) * acos( (eta_perp/eta) * cos(theta_i));
    //		float cosTheta_t = cos(theta_t);
    //
    //		float length = (2.0 * hairRadius * cos(gamma_t)) / max(1e-5, cosTheta_t);
    //
    //
    //		float invF = Fresnel(1.0/eta_perp, 1.0/eta_parr, gamma_t);
    //		float oneMinF = (1.0 - fresnel);
    //
    //		vector T = vector( 	exp(-sigma_a[0] * length),
    //					exp(-sigma_a[1] * length),
    //					exp(-sigma_a[2] * length) );
    //
    //		if( p == 1 )
    //			return color( (1.0 - fresnel) * (1.0 - invF) * T);
    //		else
    //			return color( (1.0 - fresnel) * (1.0 - invF) * invF * T*T );
    //	}
    //
    //	/**
    //	 * Absorption function returning the amount of light after the light leaves the hair fiber
    //	 * For the R-case, the Fresnel coefficient will be returned
    //	 * For the TT and TRT-cases the reflected value will also incorporate the absorption inside the hair fiber (this gives the color to TT and TRT)
    //	 */
    //	color A(uniform float p, eta_perp, eta_parr, h;)
    //	{
    //		float gamma_i = asin(h);
    //		float gamma_t = asin( h/eta_perp );
    //
    //		float F = Fresnel(eta_perp, eta_parr, gamma_i);
    //		float oneMinF = (1.0 - F);
    //
    //		//For the R component (p = 0), return the Fresnel coefficient
    //		if( p == 0 )
    //		{
    //			return color(F, F, F);
    //		}
    //
    //		//float theta_t = -sgn(theta_i) * acos( (eta_perp / eta)
    //		//float length = p * (2 * radius * cos(gamma_t)) / cos(theta_t);
    //
    //
    //		//compute transmittance
    //		//TO DO: account for eccentricity values other than a = 1.0
    //		color sigma_a2 = sigma_a;// / cos(theta_t);
    //		vector T = vector(   exp(-2.0 * sigma_a2[0] * (1.0 + cos(2.0 * gamma_t))),
    //				     exp(-2.0 * sigma_a2[1] * (1.0 + cos(2.0 * gamma_t))),
    //				     exp(-2.0 * sigma_a2[2] * (1.0 + cos(2.0 * gamma_t))) );
    //
    //		//Vector Tunit = Vector( exp(-2.0 * sigma_a2[0] * (1.0 + cos(2.0*gamma_i))),
    //		//					   exp(-2.0 * sigma_a2[1] * (1.0 + cos(2.0*gamma_i))),
    //		//					   exp(-2.0 * sigma_a2[2] * (1.0 + cos(2.0*gamma_i))) );
    //
    //		//float hairRadius = 2.0;
    //
    //		//double theta_t = sign(theta_i) * acos(cos(theta_i) * eta/x);
    //		//double l = (2.0 * hairRadius * cos(gamma_t)) / cos(theta_t);
    //
    //		//Vector T = Vector( pow(Tunit.x, hairRadius),
    //		//				   pow(Tunit.y, hairRadius),
    //		//				   pow(Tunit.z, hairRadius) );
    //
    //		//double T = exp(-2.0 * sigma_a2 * (1.0 + cos(2.0 * gamma_t)));
    //
    //		if( p == 1 )
    //			return color((oneMinF*oneMinF) * T);
    //		else
    //		{
    //			float invF = Fresnel(1.0/eta_perp, 1.0/eta_parr, gamma_t);
    //			vector T2 = T*T;//Vector(T.x*T.x, T.y*T.y, T.z*T.z);
    //
    //			return color((oneMinF*oneMinF) * invF * T2);
    //		}
    //	}
    //
    //
    //	/**
    //	 * Normal plane scattering function, accounting for the light reflection in the azimuthal direction
    //	 * [Same as RSL code]
    //	 */
    //	color Np(float p, eta_perp, eta_parr, targetAngle, theta_i; )
    //	{
    //		//retrieves the roots of the function.
    //		float roots[3] = { 0, 0, 0 };
    //		float numRoots = h(p, eta_perp, targetAngle, roots);
    //
    //		color sum = 0;
    //
    //		float r;
    //		for(  r=0; r<numRoots; r += 1 )
    //		{
    //			float gamma_i = roots[r];
    //
    //			//if angle is smaller than 90 degrees
    //			if( abs(gamma_i) < 0.5*PI )
    //			{
    //				float h = sin(gamma_i);
    //
    //				color absorption = A2(p, gamma_i, eta_perp, eta_parr, theta_i);
    //				float dexitAngle = dPhi(p, eta_perp, h);
    //				float denom = max( 1e-5, 2 * abs( dexitAngle ) );
    //				sum += (absorption / denom);
    //
    //			}
    //		}
    //
    //
    //		return sum;
    //	}
    //
    //	float targetAngle( float p, phi )
    //	{
    //		float targetAngle = abs(phi);
    //
    //		if( p != 1 )
    //		{
    //			if( targetAngle > PI )
    //				targetAngle -= 2*PI;
    //
    //			targetAngle += p*PI;
    //		}
    //
    //		return targetAngle;
    //	}
    //
    //	float exitAngle( uniform float p; float eta; float h )
    //	{
    //		// use polynomial that approximates original equation.
    //		uniform float pi3 = PI*PI*PI;
    //	        float gamma = asin(h);
    //	        float c = asin(1/eta);
    //	        return (6*p*c/PI - 2)*gamma-8*(p*c/pi3)*gamma*gamma*gamma + p*PI;
    //	}
    //
    //
    //	public color Ntrt(float eta_perp; float eta_parr; float targetAngle; float theta_i)
    //	{
    //		float hc, dh, t;
    //
    //		if( eta_perp < 2 )
    //		{
    //			float c = asin(1/eta_perp);
    //			uniform float pi3 = PI*PI*PI;
    //
    //			float gammac = sqrt( (12*c/PI - 2) / (24*(2*c/pi3)));
    //			hc = abs(sin(gammac));
    //			float ddexitAngle = d2Phi(2, eta_perp, hc);
    //			dh = min(causticLimit, 2*sqrt(2*causticWidth / abs(ddexitAngle)));
    //			t = 1;
    //		}
    //		else
    //		{
    //			hc = 0;
    //			dh = causticLimit;
    //			t = 1 - smoothstep(2, 2 + causticFade, eta_perp);
    //		}
    //
    //		float oc1 = exitAngle(2, eta_perp, hc);
    //		float oc2 = exitAngle(2, eta_perp, -hc);
    //
    //		//float phi_c = exitAngle(2, eta_perp, hc);
    //		float causticLeft = Gaussian(causticWidth, targetAngle - oc1);
    //		float causticRight = Gaussian(causticWidth, targetAngle - oc2);
    //		float causticCenter = Gaussian(causticWidth, 0);
    //
    //		//not sure about this one
    //		color glintAbsorption = A2(2, asin(hc), eta_perp, eta_parr, theta_i);
    //
    //		color L = Np(2, eta_perp, eta_parr, targetAngle, theta_i);
    //		L *= 1 - t*causticLeft/causticCenter;
    //		L *= 1 - t*causticRight/causticCenter;
    //		L += glintAbsorption * t * glintScaleFactor * dh * (causticLeft + causticRight);
    //
    //		return L;
    //	}
    //
    //	float GaussianVar( float x; float variance; )
    //	{
    //		return exp(-x*x/(2*variance)) / sqrt(2*PI*variance);
    //	}
    //
    //	/**
    //	 * Marschner reflection function, giving the amount of energy reflection according to Marschner's model (R, TT and TRT)
    //	 * 	wi, wo: in local coordinates, meaning that z-coordinate is in the normal plane,
    //	 * Return: a spectrum/color variable determining the amount of reflected colors
    //	 */
    //	public color f( float theta_i, phi_i, theta_r, phi_r; )
    //	{
    //		//compute derived angles
    //		float theta_h = 0.5 * (theta_i + theta_r);
    //		float theta_d = 0.5 * (theta_r - theta_i);
    //		float phi_h = 0.5 * (phi_i + phi_r);
    //		float phi = mod(abs(phi_r - phi_i), 2*PI);
    //		//float phi = phi_r - phi_i;
    //
    //		float eta_perp = EtaPerp(eta, theta_d);
    //		float eta_parr = EtaParr(eta, theta_d);
    //
    //		float eta_ast = EtaAsterisk(eccentricity, phi_h);
    //		float eta_perp2 = EtaPerp(eta_ast, theta_d);
    //		float eta_parr2 = EtaParr(eta_ast, theta_d);
    //
    //		color sum = 0;
    //
    //		float mtt, ntt;
    //
    //		if( enableComponents[0] > 0.0 ) {
    //                  float Mr = Gaussian(Br, theta_h - Ar);
    //                  color Nr = Np(0, eta_perp, eta_parr, targetAngle(0, phi), theta_i);
    //		  //printf("Mr * Nr: %f %f %f\n", Mr * Nr[0], Mr * Nr[1], Mr * Nr[2]);
    //                  sum += enableComponents[0] * Mr * Nr;
    //		}
    //		if( enableComponents[1] > 0.0 ) {
    //                  float Mtt = Gaussian(Btt, theta_h - Att);
    //                  color Ntt = Np(1, eta_perp, eta_parr, targetAngle(1, phi), theta_i);
    //  		  mtt = Mtt;
    //		  ntt = Ntt[0];
    //                  sum += enableComponents[1] * Mtt * Ntt;
    //                }
    //
    //		if( enableComponents[2] > 0.0 ) {
    //                  float Mtrt = Gaussian(Btrt, theta_h - Atrt);
    //                  color Ntrt = Ntrt(eta_perp2, eta_parr2, targetAngle(2, phi), theta_i);
    //                  sum += enableComponents[2] * Mtrt * Ntrt;
    //		}
    //
    //		float cos_theta_d = cos(theta_d);
    //		sum = sum / (cos_theta_d * cos_theta_d);
    //
    //		float res = mtt * ntt;
    //		//if( res > 1.0 )
    //		//	printf("mTT: %f -- nTT: %f -- mTT*nTT: %f ERROR\n", mtt, ntt, res);
    //
    //		return sum;
    //	}
    //
    //	public void physicallyCorrect()
    //	{
    //		printf("Testing physical correctness\n");
    //		uniform float i,j,k,l;
    //		uniform float samples = 400;
    //
    //		float dtheta = PI/samples;
    //		float dphi = 2*PI/samples;
    //		float dw = dtheta*dphi;
    //		color minn, maxx, avg;
    //		avg = color 0;
    //
    //		for( i = 50; i < samples; i += 1 )
    //		for( j = 0; j < samples; j += 1 )
    //		{
    //			float theta_i = -.5*PI + (i/samples) * PI;
    //			float phi_i = -PI + (j/samples) * 2*PI;
    //
    //			color sum = color 0;
    //
    //			for( k = 0; k < samples; k += 1 )
    //			for( l = 0; l < samples; l += 1 )
    //			{
    //				float theta_r = -.5*PI + (k/samples) * PI;
    //				float phi_r = -PI + (l/samples) * 2*PI;
    //				//printf("theta_r: %f -- phi_r: %f\n", theta_r, phi_r);
    //
    //				//sum += f(theta_i, phi_i, theta_r, phi_r) * cos(theta_r);
    //				sum += f( theta_i, phi_i, theta_r, phi_r ) * cos(theta_r);
    //			}
    //			sum *= cos(theta_i) * dw;
    //			if( j==0 )
    //			{
    //				minn = sum;
    //				maxx = sum;
    //			}
    //			else
    //			{
    //				minn = color( min(minn[0], sum[0]), min(minn[1], sum[1]), min(minn[2], sum[2]) );
    //				maxx = color( max(maxx[0], sum[0]), max(maxx[1], sum[1]), max(maxx[2], sum[2]) );
    //			}
    //			avg += sum;
    //
    //			//printf("integral (theta: %f,  phi: %f): %f %f %f\n", theta_i, phi_i, sum[0], sum[1], sum[2]);
    //			//printf("min: %f %f %f\nmax: %f %f %f\navg: %f %f %f\n\n", minn[0], minn[1], minn[2], maxx[0], maxx[1], maxx[2], avg[0], avg[1], avg[2]);
    //		}
    //
    //		//printf("min: %f %f %f\nmax: %f %f %f\navg: %f %f %f\n\n", minn[0], minn[1], minn[2], maxx[0], maxx[1], maxx[2], avg[0]/(samples*samples), avg[1]/(samples*samples), avg[2]/(samples*samples));
    //
    //	}
    //
    //	/**
    //	 * Marschner surface shader returning the fraction of light that is reflected to the camera/eye
    //	*/
    //	public void surface(output color Ci, Oi)
    //	{
    //		normal nn = normalize(N);   //normalized normal
    //		vector T = normalize(dPdv); //tangent
    //		vector B = normalize(nn^T); //binormal
    //		nn = normal( normalize(T^B) );
    //
    //		//compute reflected theta/phi
    //		vector wr = -normalize(I); //viewing ray
    //		vector local_wr = vector(T.wr, B.wr, nn.wr);
    //
    //		//compute spherical coords (for viewing ray)
    //		float theta_r = 0.5*PI - acos(local_wr[0]);
    //		float phi_r = atan(local_wr[1], local_wr[2]);
    //
    //		Ci = 0;
    //
    //		//sample from all lights
    //		illuminance( P )
    //		{
    //			vector wi = normalize(L); //light ray
    //			vector local_wi = vector(T.wi, B.wi, N.wi);
    //
    //			//computing spherical coords (for ray towards light source)
    //			float theta_i = 0.5*PI - acos(local_wi[0]);
    //			float phi_i = atan(local_wi[1], local_wi[2]);
    //
    //			//float opacity
    //
    //			Ci += f(theta_i, phi_i, theta_r, phi_r) * cos(theta_i) * Cl
    //				+ Kd * (1.0 - sigma_a) * abs(nn.wi) * Cl;
    //
    //		}
    //
    //		//incorporating opacity
    //		Oi = 1;
    //		Ci *= Oi;
    //	}
    //}

} // namespace pbrt
