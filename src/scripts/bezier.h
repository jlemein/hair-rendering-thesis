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
#include <iostream>
#include <stdlib.h>

struct Point3 {
public:
    Point3(double x = 0.0, double y = 0.0, double z = 0.0f) : x(x), y(y), z(z) {
    }
    
    double x, y, z;
    Point3 operator*(double m) {
        return Point3(this->x * m, this->y*m, this->z * m);
    }
    Point3 operator+(const Point3& p) {
        return Point3(this->x + p.x, this->y + p.y, this->z + p.z);
    }
    
    Point3 operator-(const Point3& p) {
        return Point3(this->x - p.x, this->y - p.y, this->z - p.z);
    }
    
    bool operator==(const Point3& p) const {
        // compare by delta instead of == to prevent floating point precision errors
        return abs(p.x - x) + abs(p.y - y) + abs(p.z - z) < 0.0001;    
    }
    
    // static methods
    static double DistanceBetween(const Point3& p1, const Point3& p2);
    
    
    friend std::ostream& operator<<(std::ostream& out, const Point3& p);
};


class BezierSpline {
private:
    bool mShareControlPoints = false;
    std::vector<Point3> mControlPoints;
    
public:
    BezierSpline();
    BezierSpline(const BezierSpline& spline);
    BezierSpline(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3);
    BezierSpline(const double controlPoints[], unsigned int size);
    BezierSpline(const Point3 controlPoints[], unsigned int size);
    
    void addControlPoint(const Point3& p);
    void addControlPoints(const double controlPoints[], unsigned int size);
    void addControlPoints(const Point3 controlPoints[], unsigned int size);
    void addControlPoints(const std::vector<Point3>& controlPoints);
    
    bool setUseSharedControlPoints(bool useSharedControlPoints);
    bool isUsingSharedControlPoints() const;
    
    unsigned int getSegmentCount() const;
    unsigned int getControlPointOffsetForSegment(unsigned int segment) const;
    
    std::vector<Point3>& getControlPoints();
    const std::vector<Point3>& getControlPoints() const;
    
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

