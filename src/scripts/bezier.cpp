/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "bezier.h"

#include <algorithm>
#include <iostream>

Point3 cubicBezier2(double t, Point3 p0, Point3 p1, Point3 p2, Point3 p3) {

    Point3 a = p0 * (1.0 - t)*(1.0 - t)*(1.0 - t);
    Point3 b = p1 * 3.0 * t * (1.0 - t)*(1.0 - t);
    Point3 c = p2 * 3.0 * t * t * (1.0 - t)*(1.0 - t);
    Point3 d = p3 * t * t*t;

    return a + b + c + d;
}

Point3 cubicBezier(double t, Point3 p0, Point3 p1, Point3 p2, Point3 p3) {

    Point3 a = p0 * (1.0 - t)*(1.0 - t)*(1.0 - t);
    Point3 b = p0 + (p1 - p0) * 3.0 * t * (1.0 - t)*(1.0 - t);
    Point3 c = p3 + (p2 - p3) * 3.0 * t * t * (1.0 - t)*(1.0 - t);
    Point3 d = p3 * t * t*t;

    return a + b + c + d;
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
    return cubicBezier(t, mControlPoints[controlPointStart], mControlPoints[controlPointStart + 1], mControlPoints[controlPointStart + 2], mControlPoints[controlPointStart + 3]);
}