/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "bezier.h"

#include <algorithm>
#include <iostream>

static Point3 cubicBezier(double t, Point3 p0, Point3 p1, Point3 p2, Point3 p3) {

    Point3 a = p0 * (1.0 - t)*(1.0 - t)*(1.0 - t);
    Point3 b = p1 * 3.0 * t * (1.0 - t)*(1.0 - t);
    Point3 c = p2 * 3.0 * t * t * (1.0 - t)*(1.0 - t);
    Point3 d = p3 * t * t*t;

    return a + b + c + d;
}

static Point3 cubicBezier2(double t, Point3 p0, Point3 p1, Point3 p2, Point3 p3) {

    Point3 a = p0 * (1.0 - t)*(1.0 - t)*(1.0 - t);
    Point3 b = (p1 - p0) * 3.0 * t * (1.0 - t)*(1.0 - t);
    Point3 c = (p2 - p3) * 3.0 * t * t * (1.0 - t)*(1.0 - t);
    Point3 d = p3 * t * t*t;

    return a + b + c + d;
}

inline Point3 Lerp(double t, Point3 p1, Point3 p2) {
    return p1 * (1 - t) + p2 * t;
}

static Point3 EvalBezier(const Point3 cp[4], double u) {
    Point3 cp1[3] = {Lerp(u, cp[0], cp[1]), Lerp(u, cp[1], cp[2]), Lerp(u, cp[2], cp[3])};
    Point3 cp2[2] = {Lerp(u, cp1[0], cp1[1]), Lerp(u, cp1[1], cp1[2])};

    return Lerp(u, cp2[0], cp2[1]);
}

std::ostream& operator<<(std::ostream& out, const Point3& p) {
    out << "[" << p.x << " " << p.y << " " << p.z << "]";
}

double Point3::DistanceBetween(const Point3& p1, const Point3& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

BezierSpline::BezierSpline() : mShareControlPoints(false) {
};

BezierSpline::BezierSpline(const BezierSpline& spline) {
    this->mShareControlPoints = spline.mShareControlPoints;
    this->mControlPoints.insert(this->mControlPoints.begin(), spline.mControlPoints.begin(), spline.mControlPoints.end());
}

BezierSpline::BezierSpline(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3) : BezierSpline() {
    this->mControlPoints.push_back(p0);
    this->mControlPoints.push_back(p1);
    this->mControlPoints.push_back(p2);
    this->mControlPoints.push_back(p3);
};

BezierSpline::BezierSpline(const double controlPoints[], unsigned int size) : BezierSpline() {
    for (int i = 0; i + 2 < size; i += 3) {
        this->mControlPoints.push_back(Point3(controlPoints[i], controlPoints[i + 1], controlPoints[i + 2]));
    }
}

BezierSpline::BezierSpline(const Point3 controlPoints[], unsigned int size) : BezierSpline() {
    for (int i = 0; i < size; ++i) {
        this->mControlPoints.push_back(controlPoints[i]);
    }
}

void BezierSpline::addControlPoint(const Point3& p) {
    this->mControlPoints.push_back(p);
}

void BezierSpline::addControlPoints(const double controlPoints[], unsigned int size) {
    for (int i = 0; i + 2 < size; i += 3) {
        this->mControlPoints.push_back(Point3(controlPoints[i], controlPoints[i + 1], controlPoints[i + 2]));
    }
}

void BezierSpline::addControlPoints(const Point3 controlPoints[], unsigned int size) {
    for (int i = 0; i < size; ++i) {
        this->mControlPoints.push_back(controlPoints[i]);
    }
}

void BezierSpline::addControlPoints(const std::vector<Point3>& controlPoints) {
    //TODO: make an option to remove first control point or not. This assumption is not clear from function signature
    this->mControlPoints.insert(this->mControlPoints.end(),
            this->mShareControlPoints ? controlPoints.begin() + 1 : controlPoints.begin(),
            controlPoints.end());
}

bool BezierSpline::setUseSharedControlPoints(bool useSharedControlPoints) {
    this->mShareControlPoints = useSharedControlPoints;
}

bool BezierSpline::isUsingSharedControlPoints() const {
    return this->mShareControlPoints;
}

unsigned int BezierSpline::getSegmentCount() const {
    if (this->mShareControlPoints) {
        return (mControlPoints.size() - 1) / 3;
    } else {
        return mControlPoints.size() / 4;
    }
}

std::vector<Point3>& BezierSpline::getControlPoints() {
    return this->mControlPoints;
}

const std::vector<Point3>& BezierSpline::getControlPoints() const {
    return this->mControlPoints;
}

unsigned int BezierSpline::getControlPointCount() const {
    return mControlPoints.size();
}

Point3 BezierSpline::sampleCurve(double t) const {
    if (t < 0.0 || t > 1.0) {
        std::cout << "[WARNING]: t value invalid. Value will be clamped between [0, 1]" << std::endl;
        t = std::max(0.0, std::min(t, 1.0));
    }


    double integerPart;
    double u = std::modf(t * getSegmentCount(), &integerPart);

    unsigned int segmentIndex;
    if (integerPart >= getSegmentCount()) {
        segmentIndex = getSegmentCount() - 1;
        u = 1.0;
    } else {
        segmentIndex = static_cast<unsigned int> (integerPart);
    }

    return this->sampleSegment(segmentIndex, u);
}

unsigned int BezierSpline::getControlPointOffsetForSegment(unsigned int segment) const {
    return this->mShareControlPoints ? segment * 3 : segment * 4;
}

Point3 BezierSpline::sampleSegment(unsigned int segment, double t) const {
    if (segment < 0 || segment >= this->getSegmentCount()) {
        std::cout << "[ERROR]: Invalid segment index " << segment << ", segment count is " << getSegmentCount() << std::endl;
        exit(1);
    }

    if (t < 0.0 || t > 1.0) {
        std::cout << "[WARNING]: Sampling curve with invalid t value. t will be clamped between [0, 1]" << std::endl;
        t = std::max(0.0, std::min(t, 1.0));
    }

    unsigned int controlPointStart = getControlPointOffsetForSegment(segment);

    const Point3 cp[4] = {mControlPoints[controlPointStart], mControlPoints[controlPointStart + 1], mControlPoints[controlPointStart + 2], mControlPoints[controlPointStart + 3]};
    return EvalBezier(cp, t);
    //return cubicBezier(t, mControlPoints[controlPointStart], mControlPoints[controlPointStart + 1], mControlPoints[controlPointStart + 2], mControlPoints[controlPointStart + 3]);
}