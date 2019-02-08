/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   OpenVdbReader.cpp
 * Author: jeffrey
 *
 * Created on February 1, 2019, 4:25 PM
 */

#include <exception>
#include <sstream>
#include <string>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include "OpenVdbReader.h"
#include "../hairstruct.h"

using namespace openvdb;

#ifdef __APPLE__
using namespace openvdb::v6_0::tools;
#else
using namespace openvdb::v3_1::tools;
#endif

OpenVdbReader::OpenVdbReader(std::string fileName)
: mInputFileName(fileName), mBoundingBoxMin(0.0, 0.0, 0.0), mBoundingBoxMax(0.0, 0.0, 0.0) {
}

OpenVdbReader::OpenVdbReader(const OpenVdbReader& orig) {
}

OpenVdbReader::~OpenVdbReader() {
}

void OpenVdbReader::initialize() {
    openvdb::initialize();
    openvdb::io::File file(mInputFileName);
    file.open();

    for (openvdb::io::File::NameIterator nameIter = file.beginName();
            nameIter != file.endName(); ++nameIter) {
        std::cout << "Found grid with name: " << nameIter.gridName() << std::endl;
        openvdb::GridBase::Ptr grid = file.readGrid(nameIter.gridName());
        printMetaDataForGrid(grid);

        mGrids.push_back(grid);

        if (nameIter.gridName() == "HairDensityGrid") {
            mHairDensityGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(grid);

            try {
                mVoxelSize = grid->metaValue<float>("VoxelSize");
            } catch (std::exception& e) {
                std::cout << "[ERROR]: No VoxelSize value could be found" << std::endl;
            }

            try {
                std::string bbMinStr = grid->metaValue<std::string>("BoundingBoxMin");
                std::string bbMaxStr = grid->metaValue<std::string>("BoundingBoxMax");


                std::stringstream bbMin(bbMinStr);
                std::stringstream bbMax(bbMaxStr);
                bbMin >> this->mBoundingBoxMin.x >> this->mBoundingBoxMin.y >> this->mBoundingBoxMin.z;
                bbMax >> this->mBoundingBoxMax.x >> this->mBoundingBoxMax.y >> this->mBoundingBoxMax.z;
            } catch (std::exception& e) {
                std::cout << "[ERROR]: No bounding box value is stored in hair density grid" << std::endl;
            }
        }
    }

    file.close();
}

openvdb::FloatGrid::Ptr OpenVdbReader::getHairDensityGrid() const {
    return mHairDensityGrid;
}

void OpenVdbReader::getBoundingBox(Point3* bbMin, Point3* bbMax) const {
    bbMin->x = this->mBoundingBoxMin.x;
    bbMin->y = this->mBoundingBoxMin.y;
    bbMin->z = this->mBoundingBoxMin.z;

    bbMax->x = this->mBoundingBoxMax.x;
    bbMax->y = this->mBoundingBoxMax.y;
    bbMax->z = this->mBoundingBoxMax.z;
}

float OpenVdbReader::getVoxelSize() const {
    return mVoxelSize;
}

/**
 * Interpolates the voxel grid from a 3D point in space to another point
 * in 3D space
 * @param from
 * @param to
 * @return
 */
float OpenVdbReader::interpolate(const Point3& from, const Point3& to) {
    // there is a choice of different interpolators, mainly PointSampler, BoxSampler and QuadraticSampler
    // in addition to StaggeredPointSampler, StaggeredBoxSampler and StaggeredQuadraticSampler for staggered velocity grids.

    const FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(mHairDensityGrid);
    Vec3d wsFrom(from.x, from.y, from.z);
    Vec3d wsTo(to.x, to.y, to.z);

    wsFrom /= mVoxelSize;
    wsTo /= mVoxelSize;

    // Request a value accessor for accelerated access.
    // (Because value accessors employ a cache, it is important to declare
    // one accessor per thread.)
    FloatGrid::ConstAccessor accessor = grid->getConstAccessor();

    // Instantiate the GridSampler template on the accessor type and on
    // a box sampler for accelerated trilinear interpolation.
    GridSampler<FloatGrid::ConstAccessor, PointSampler> fastSampler(accessor, grid->transform());

    int nSamples = 100;
    Vec3d sampleIncrement = (wsTo - wsFrom) / nSamples;
    openvdb::FloatGrid::ValueType value = 0.0f;

    std::cout << wsFrom << " to " << wsTo << std::endl;

    for (int i = 0; i < 100; ++i) {
        Vec3d wsSamplePosition = wsFrom + i * sampleIncrement;
        value += fastSampler.wsSample(wsSamplePosition); // / mVoxelSize);
    }

    // normalize the sampling
    std::cout << "mVoxelSize" << mVoxelSize << std::endl;
    double numVoxelCellsCrossed = Point3::DistanceBetween(from, to) / mVoxelSize;
    double integratedValue = numVoxelCellsCrossed * (value / static_cast<double> (nSamples));

    return integratedValue;
}

float OpenVdbReader::interpolateToInfinity(const Point3& from, const Point3& direction) {
    std::cout << "[WARNING]: interpolateToInfinity(Point3&, Point$) is not implemented\n";

    //return interpolate(from, direction.normalize() * mInfinity);
    return -1.0f;
}

void OpenVdbReader::printMetaDataForGrid(openvdb::GridBase::Ptr grid) const {
    std::cout << "Metadata for grid with name '" << grid->getName() << "':\n";
    for (openvdb::MetaMap::MetaIterator iter = grid->beginMeta();
            iter != grid->endMeta(); ++iter) {
        const std::string& name = iter->first;
        openvdb::Metadata::Ptr value = iter->second;
        std::string valueAsString = value->str();
        std::cout << "\t" << name << " = " << valueAsString << std::endl;
    }
}

void OpenVdbReader::printMetaDataForAllGrids() const {
    for (auto& grid : mGrids) {
        printMetaDataForGrid(grid);
        std::cout << std::endl;
    }
}

void OpenVdbReader::printMetaDataForHairDensityGrid() const {
    if (mHairDensityGrid != 0) {
        printMetaDataForGrid(mHairDensityGrid);
    } else {
        std::cout << "Hair density grid does not exist\n";
    }
}

void OpenVdbReader::printVdbInfo() const {
    for (auto& grid : mGrids) {
        std::cout << "Voxel info for grid with name '" << grid->getName() << "'\n\n";
        printMetaDataForGrid(grid);
    }
}
