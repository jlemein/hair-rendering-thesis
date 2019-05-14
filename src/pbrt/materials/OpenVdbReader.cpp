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
#include "geometry.h"


using namespace openvdb;

#ifdef __APPLE__
using namespace openvdb::v6_0::tools;
#else
using namespace openvdb::v3_1::tools;
#endif

namespace pbrt {

    OpenVdbReader::OpenVdbReader(std::string fileName)
    : mInputFileName(fileName) {
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

                    Point3<Float> p1, p2;
                    bbMin >> p1.x >> p1.y >> p1.z;
                    bbMax >> p2.x >> p2.y >> p2.z;
                    mBounds = Bounds3<Float>(p1, p2);

                } catch (const std::exception& e) {
                    std::cout << "[ERROR]: No bounding box value is stored in hair density grid" << std::endl;
                }
            }
        }

        file.close();
    }

    openvdb::FloatGrid::Ptr OpenVdbReader::getHairDensityGrid() const {
        return mHairDensityGrid;
    }

    const Bounds3<Float>& OpenVdbReader::getBounds() const {
        return mBounds;
    }

    Float OpenVdbReader::getVoxelSize() const {
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
    Float OpenVdbReader::interpolate(const Vector3f& from, const Vector3f& to, unsigned int sampleCount) {
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
        double numVoxelCellsCrossed = Vector3f(from - to).Length() / this->getVoxelSize();
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
    InterpolationResult OpenVdbReader::interpolate(const Vector3f& from, const Vector3f& to, unsigned int sampleCount) const {
        // there is a choice of different interpolators, mainly PointSampler, BoxSampler and QuadraticSampler
        // in addition to StaggeredPointSampler, StaggeredBoxSampler and StaggeredQuadraticSampler for staggered velocity grids.

        // convert to openvdb world space coordinates
        Vec3d wsFrom(from.x, from.y, from.z);
        Vec3d wsTo(to.x, to.y, to.z);
        Vec3d isFrom = mHairDensityGrid->transform().worldToIndex(wsFrom);
        Vec3d isTo = mHairDensityGrid->transform().worldToIndex(wsTo);
        //        std::cout << "Interpolating from world space " << wsFrom << " -> " << wsTo << std::endl;
        //        std::cout << "Interpolating from index space " << isFrom << " -> " << isTo << std::endl;

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
        double numVoxelCellsCrossed = Vector3f(from - to).Length() / this->getVoxelSize();
        //        std::cout << "Voxels crossed: " << numVoxelCellsCrossed << std::endl;
        //        std::cout << "Value: " << value << std::endl;

        InterpolationResult result;
        result.scatterCount = numVoxelCellsCrossed * (value / static_cast<double> (sampleCount));
        result.averageThetaD = 1.2;

        return result;
    }

    InterpolationResult OpenVdbReader::interpolateToInfinity(const Vector3f& from, const Vector3f& direction) const {
        Vector3f normalizedDirection = direction / direction.Length();
        Vector3f to = from + normalizedDirection * mBounds.Diagonal().Length();
        return this->interpolate(from, to);
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
}