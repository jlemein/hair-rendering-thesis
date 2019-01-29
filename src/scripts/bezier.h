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
    BezierSpline(Point3 p0, Point3 p1, Point3 p2, Point3 p3);
    BezierSpline(const double controlPoints[], unsigned int size);
    
    /**
     * Samples the curve
     * @param t
     * @return 
     */
    Point3 sample(double t) const;
};

#endif /* BEZIER_H */

