/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OpenVdbReader.h
 * Author: jeffrey
 *
 * Created on February 1, 2019, 4:25 PM
 */

#ifndef OPENVDBREADER_H
#define OPENVDBREADER_H

#include <string>
#include <vector>
#include "../bezier.h"
#include <openvdb/openvdb.h>

struct InterpolationResult {
    float scatterCount;
    float averageThetaD;
};
class OpenVdbReader {
public:
    OpenVdbReader(std::string fileName);
    OpenVdbReader(const OpenVdbReader& orig);
    virtual ~OpenVdbReader();
    
    void initialize();
    
    /**
     * Interpolates the voxel grid from a 3D point in space to another point 
     * in 3D space
     * @param from
     * @param to
     * @return 
     */
    float interpolate(const Point3& from, const Point3& to, unsigned int sampleCount = 100);
    InterpolationResult interpolate(const Point3& from, const Point3& to, unsigned int sampleCount = 100) const;
    float interpolateToInfinity(const Point3& from, const Point3& direction);
    
    void printMetaDataForAllGrids() const;
    void printMetaDataForHairDensityGrid() const;
   
    void printVdbInfo() const;
    void getBoundingBox(Point3* from, Point3* p2) const;
    double getVoxelSize() const;
    openvdb::FloatGrid::Ptr getHairDensityGrid() const;
    
private:
    std::string mInputFileName;
    std::vector<openvdb::GridBase::Ptr> mGrids;
    openvdb::FloatGrid::Ptr mHairDensityGrid;
    float mInfinity = 1000.0f;
    
    Point3 mBoundingBoxMin, mBoundingBoxMax;
    
     void printMetaDataForGrid(openvdb::GridBase::Ptr) const;


};

#endif /* OPENVDBREADER_H */

