//#include "pbrt.h"
//#include "sampler.h"
//#include <atomic>
//#include <vector>
//#include <functional>
//#include <random>
//
//using namespace pbrt;
//
//namespace jlemein {
//
//    double _integrateMonteCarloFront(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {
//
//        std::default_random_engine generator;
//        std::uniform_real_distribution<double> distribution(0.0, 1.0);
//
//        double V = 2.0 * Pi;
//        double sum = .0;
//
//        for (int i = 0; i < nSamples; ++i) {
//            //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
//            Vector3f wr = SampleFrontHemisphere(Point2f(distribution(generator), distribution(generator)));
//            sum += f(wr, wi);
//        }
//
//        return V / static_cast<double> (nSamples) * sum;
//    }
//
//    double _integrateMonteCarloBack(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {
//
//        std::default_random_engine generator;
//        std::uniform_real_distribution<double> distribution(0.0, 1.0);
//
//        double V = 2.0 * Pi;
//        double sum = .0;
//
//        for (int i = 0; i < nSamples; ++i) {
//            //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
//            Vector3f wr = SampleBackHemisphere(Point2f(distribution(generator), distribution(generator)));
//            sum += f(wr, wi);
//        }
//
//        return V / static_cast<double> (nSamples) * sum;
//    }
//
//    Vector3f UniformSampleSphere(const Point2f &u) {
//        Float z = 1 - 2 * u[0];
//        Float r = std::sqrt(std::max((Float) 0, (Float) 1 - z * z));
//        Float phi = 2 * Pi * u[1];
//        return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
//    }
//
//    double _integrateMonteCarlo(const Vector3f& wi, std::function<double(const Vector3f&, const Vector3f&) > f, int nSamples) {
//
//        std::default_random_engine generator;
//        std::uniform_real_distribution<double> distribution(0.0, 1.0);
//
//        double V = 4.0 * Pi;
//        double sum = .0;
//
//        for (int i = 0; i < nSamples; ++i) {
//            //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
//            Vector3f wr = UniformSampleSphere(Point2f(distribution(generator), distribution(generator)));
//            sum += f(wr, wi);
//        }
//
//        return V / static_cast<double> (nSamples) * sum;
//    }
//
//    double _integrateMonteCarlo(std::function<double(const Float) > f, int nSamples) {
//
//        std::default_random_engine generator;
//        std::uniform_real_distribution<double> distribution(.0, 1.0);
//
//        double V = Pi;
//        double sum = .0;
//
//        for (int i = 0; i < nSamples; ++i) {
//            //Vector3f wi = Vector3f(-1.0, 0.0, 0.0);
//            Float x = -.5 * Pi + distribution(generator) * Pi;
//            sum += f(x);
//        }
//
//        return V / static_cast<double> (nSamples) * sum;
//    }
//}