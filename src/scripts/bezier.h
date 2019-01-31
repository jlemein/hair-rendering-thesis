/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bezier.h
 * Author: jeffrey
 *
 * Created on January 29, 2019, 8:40 PM
 */


#ifndef BEZIER_H
#define BEZIER_H

#include <vector>

struct Point3 {
public:
    double x, y, z;
    Point3 operator*(double m) {
        return Point3(this->x * m, this->y*m, this->z * m);
    }
    Point3 operator+(const Point3& p) {
        return Point3(this->x + p.x, this->y + p.y, this->z + p.z);
    }
    Point3(double x, double y, double z) : x(x), y(y), z(z) {
    }
};


class BezierSpline {
private:
    std::vector<Point3> mControlPoints;
    
public:
    BezierSpline(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3);
    BezierSpline(const double controlPoints[], unsigned int size);
    BezierSpline(const Point3 controlPoints[], unsigned int size);
    
    void addControlPoints(const double controlPoints[], unsigned int size);
    void addControlPoints(const Point3 controlPoints[], unsigned int size);
    void addControlPoint(const Point3& p);
    
    unsigned int getSegmentCount() const;
    const Point3* getControlPoints() const;
    unsigned int getControlPointCount() const;
    
    Point3 sampleCurve(double t) const;
    
    /**
     * samples a specific segment
     * @param segment
     * @param t the point on the curve to be sampled between start (t=0) and end (t=1)
     * @return 
     */
    Point3 sampleSegment(unsigned int segment, double t) const;
};

#endif /* BEZIER_H */

