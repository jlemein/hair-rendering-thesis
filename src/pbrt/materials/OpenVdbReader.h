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

#include "pbrt.h"
#include "material.h"
#include "marschner.h"
#include <string>
#include <vector>
#include <openvdb/openvdb.h>

namespace pbrt {
    
    struct InterpolationResult {
        Float scatterCount;
        Float averageThetaD;
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
        Float interpolate(const Point3f& from, const Point3f& to, unsigned int sampleCount = 100);
        InterpolationResult interpolate(const Point3f& from, const Point3f& to, unsigned int sampleCount = 100) const;
        
        InterpolationResult interpolateToInfinity(const Point3f& from, const Vector3f& direction) const;

        void printMetaDataForAllGrids() const;
        void printMetaDataForHairDensityGrid() const;

        void printVdbInfo() const;
        const Bounds3<Float>& getBounds() const;
        Float getVoxelSize() const;
        openvdb::FloatGrid::Ptr getHairDensityGrid() const;

    private:
        std::string mInputFileName;
        std::vector<openvdb::GridBase::Ptr> mGrids;
        openvdb::FloatGrid::Ptr mHairDensityGrid;
        Float mInfinity = 1000.0f;

        //Vector3f mBoundingBoxMin, mBoundingBoxMax;
        Bounds3f mBounds;

         void printMetaDataForGrid(openvdb::GridBase::Ptr) const;


    };
}

#endif /* OPENVDBREADER_H */

