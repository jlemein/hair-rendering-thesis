#ifndef HAIRSTRUCT_H
#define HAIRSTRUCT_H

/**
 * @author Jeffrey Lemein
 */
#include <vector>
#include <iostream>
#include <string>

#include "bezier.h"

///**
// * @brief The Point struct represents a point in 3D space
// */
//struct Point {
//    /** @brief Position of the point */
//    float x, y, z;
//
//    friend std::ostream& operator<<(std::ostream& out, const Point& p);
//};

enum CurveType { FLAT, CYLINDER, RIBBON };

/**
 * @brief The Curve struct representing a curve in 3D space
 */
struct HairFiber {
    /** @brief type of curve */
    std::string type;

    /** @brief list of points representing the curve */
    //std::vector<Point> points;
    BezierSpline curve;

    /**
     * @brief widths for the cross section of the curve (to allow oval curves)
     */
    float width0, width1;
};

/**
 * @brief The Hair struct represents a hair model file read from pbrt file.
 * Actually its just reading a scene file (*.pbrt), but assumes only curve
 * information is present corresponding to a hair model.
 */
class Hair
{
public:
    Hair();
    
    /** @brief list of curves representing the hair strands of the model */
    std::vector<HairFiber> fibers;
    
    CurveType curveType;

    /** @brief sceneExtentMin orthogonal boundary for xyz dimensions */
    Point3 sceneExtentMin, sceneExtentMax;

    /** @brief orthogonal size for each dimension in xyz */
    Point3 size;

    /** @brief center Center point of bounding box of hair model */
    Point3 center;

    float minCurveWidth0, minCurveWidth1, maxCurveWidth0, maxCurveWidth1;
    float avgCurveWidth0, avgCurveWidth1;

    void writeFile(const std::string& outputFileName);

    friend std::istream& operator>>(std::istream& istream, Hair& hair);
    friend std::ostream& operator<<(std::ostream& ostream, Hair& hair);
    
    std::string getCurveType() const;
    
    void optimizeCurves();
};

#endif // HAIRSTRUCT_H
