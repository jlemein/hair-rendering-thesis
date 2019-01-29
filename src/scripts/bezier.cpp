/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "bezier.h"

#include <algorithm>
#include <iostream>

//void quadB(double t, Point3 p0, Point3 p1, Point3 p2) {
//    return p1 + (1.0 - t)*(1.0 - t)*(p0 - p1) + t * t * (p2 - p1);
//}

Point3 cubicBezier(double t, Point3 p0, Point3 p1, Point3 p2, Point3 p3) {
    Point3 a = p0 * (1.0 - t)*(1.0 - t)*(1.0 - t);
    Point3 b = p1 * 3.0 * (1.0 - t)*(1.0 - t) * t;
    Point3 c = p2 * 3.0 * (1.0 - t)*(1.0 - t) * t * t;
    Point3 d = p3 * t * t * t;

    return a + b + c + d;
}

BezierSpline::BezierSpline(Point3 p0, Point3 p1, Point3 p2, Point3 p3) {
    this->mControlPoints.push_back(p0);
    this->mControlPoints.push_back(p1);
    this->mControlPoints.push_back(p2);
    this->mControlPoints.push_back(p3);
};

BezierSpline::BezierSpline(const double controlPoints[], unsigned int size) {
    for (int i = 0; i + 2 < size; i += 3) {
        this->mControlPoints.push_back(Point3(controlPoints[i], controlPoints[i + 1], controlPoints[i + 2]));
    }
}

Point3 BezierSpline::sample(double t) const {
    if (t < 0.0 || t > 1.0) {
        std::cout << "WARNING: t = " << t << ": invalid t value. Value will be clamped between 0 and 1" << std::endl;
        t = std::max(0.0, std::min(1.0, t));
    }
    return cubicBezier(t, this->mControlPoints[0], this->mControlPoints[1], this->mControlPoints[2], this->mControlPoints[3]);
}