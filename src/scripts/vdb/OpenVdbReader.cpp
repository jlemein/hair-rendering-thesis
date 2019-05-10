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

double OpenVdbReader::getVoxelSize() const {
    return mHairDensityGrid->transform().voxelSize().x();
}

/**
 * Interpolates the voxel grid from a 3D point in space to another point
 * in 3D space.
 * The point of origin is not taken into account when interpolating. The samples are divided
 * equally over the line, from origin to destination, including the destination point, but not including the origin point.
 * @param from
 * @param to
 * @param sampleCount The amount of samples to take along the ray from origin to destination
 * @return
 */
float OpenVdbReader::interpolate(const Point3& from, const Point3& to, unsigned int sampleCount) {
    // there is a choice of different interpolators, mainly PointSampler, BoxSampler and QuadraticSampler
    // in addition to StaggeredPointSampler, StaggeredBoxSampler and StaggeredQuadraticSampler for staggered velocity grids.

    // convert to openvdb world space coordinates
    Vec3d wsFrom(from.x, from.y, from.z);
    Vec3d wsTo(to.x, to.y, to.z);
    Vec3d isFrom = mHairDensityGrid->transform().worldToIndex(wsFrom);
    Vec3d isTo = mHairDensityGrid->transform().worldToIndex(wsTo);
    std::cout << "Interpolating from world space " << wsFrom << " -> " << wsTo << std::endl;
    std::cout << "Interpolating from index space " << isFrom << " -> " << isTo << std::endl;

    const FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(mHairDensityGrid);

    // Request a value accessor for accelerated access.
    // (Because value accessors employ a cache, it is important to declare one accessor per thread.)
    FloatGrid::ConstAccessor accessor = grid->getConstAccessor();
    GridSampler<FloatGrid::ConstAccessor, BoxSampler> fastSampler(accessor, grid->transform());

    Vec3d stepIncrement = (wsTo - wsFrom) / sampleCount;
    openvdb::FloatGrid::ValueType value = 0.0f;

    for (int i = 1; i <= sampleCount; ++i) {
        Vec3d wsSamplePosition = wsFrom + i * stepIncrement;
        //FloatGrid::ValueType vv = fastSampler.wsSample(wsSamplePosition);
        //std::cout << "Sampling at index: [" << grid->transform().worldToIndex(wsSamplePosition) << "] = " << vv << std::endl;
        FloatGrid::ValueType sampledValue = fastSampler.wsSample(wsSamplePosition);
        if (sampledValue < 0.0) {
            sampledValue = 0.0;
        }
        value += sampledValue;
    }

    // normalize the sampling
    double numVoxelCellsCrossed = Point3::DistanceBetween(from, to) / this->getVoxelSize();
    std::cout << "Voxels crossed: " << numVoxelCellsCrossed << std::endl;
    std::cout << "Value: " << value << std::endl;
    double integratedValue = numVoxelCellsCrossed * (value / static_cast<double> (sampleCount));

    return integratedValue;
}

/**
 * Interpolates the voxel grid from a 3D point in space to another point
 * in 3D space.
 * The point of origin is not taken into account when interpolating. The samples are divided
 * equally over the line, from origin to destination, including the destination point, but not including the origin point.
 * @param from
 * @param to
 * @param sampleCount The amount of samples to take along the ray from origin to destination
 * @return
 */
InterpolationResult OpenVdbReader::interpolate(const Point3& from, const Point3& to, unsigned int sampleCount) const {
    // there is a choice of different interpolators, mainly PointSampler, BoxSampler and QuadraticSampler
    // in addition to StaggeredPointSampler, StaggeredBoxSampler and StaggeredQuadraticSampler for staggered velocity grids.

    // convert to openvdb world space coordinates
    Vec3d wsFrom(from.x, from.y, from.z);
    Vec3d wsTo(to.x, to.y, to.z);
    Vec3d isFrom = mHairDensityGrid->transform().worldToIndex(wsFrom);
    Vec3d isTo = mHairDensityGrid->transform().worldToIndex(wsTo);
    std::cout << "Interpolating from world space " << wsFrom << " -> " << wsTo << std::endl;
    std::cout << "Interpolating from index space " << isFrom << " -> " << isTo << std::endl;

    const FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(mHairDensityGrid);

    // Request a value accessor for accelerated access.
    // (Because value accessors employ a cache, it is important to declare one accessor per thread.)
    FloatGrid::ConstAccessor accessor = grid->getConstAccessor();
    GridSampler<FloatGrid::ConstAccessor, BoxSampler> fastSampler(accessor, grid->transform());

    Vec3d stepIncrement = (wsTo - wsFrom) / sampleCount;
    openvdb::FloatGrid::ValueType value = 0.0f;

    for (int i = 1; i <= sampleCount; ++i) {
        Vec3d wsSamplePosition = wsFrom + i * stepIncrement;
        //FloatGrid::ValueType vv = fastSampler.wsSample(wsSamplePosition);
        //std::cout << "Sampling at index: [" << grid->transform().worldToIndex(wsSamplePosition) << "] = " << vv << std::endl;
        FloatGrid::ValueType sampledValue = fastSampler.wsSample(wsSamplePosition);
        if (sampledValue < 0.0) {
            sampledValue = 0.0;
        }
        value += sampledValue;
    }

    // normalize the sampling
    double numVoxelCellsCrossed = Point3::DistanceBetween(from, to) / this->getVoxelSize();
    std::cout << "Voxels crossed: " << numVoxelCellsCrossed << std::endl;
    std::cout << "Value: " << value << std::endl;

    InterpolationResult result;
    result.scatterCount = numVoxelCellsCrossed * (value / static_cast<double> (sampleCount));
    result.averageThetaD = 0.0;

    return result;
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
