/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <iostream>
#include <fstream>
#include <string>
#include "../hairstruct.h"
#include "OpenVdbReader.h"
using namespace std;

openvdb::FloatGrid::Ptr grid;

void listCube() {

}

void printCube(openvdb::FloatGrid::Accessor& accessor) {
    // Define a coordinate with large signed indices.
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (double k = 0; k < 3; ++k) {
                openvdb::Coord ijk(i, j, k);
                cout << "[ " << i << "; " << j << "; " << k << " ] = " << accessor.getValue(ijk) << endl;
            }
        }
    }
}

void printCubeWorld(openvdb::FloatGrid::Accessor& accessor) {
    // Define a coordinate with large signed indices.
    for (double i = 0.5; i < 3.0; i += 1.0) {
        for (double j = 0.5; j < 3.0; j += 1.0) {
            for (double k = 0.5; k < 3.0; k += 1.0) {
                openvdb::Coord ijk(i, j, k);
                cout << "[" << i << "; " << j << "; " << k << "] = " << accessor.getValue(ijk) << endl;
            }
        }
    }
}

void indexSampler() {
    const openvdb::Vec3R ijk(0.0, 1.0, 1.0);
    // Compute the value of the grid at ijk via nearest-neighbor (zero-order)
    // interpolation.
    openvdb::FloatGrid::ValueType v0 = openvdb::v3_1::tools::PointSampler::sample(grid.get()->tree(), openvdb::Vec3R(0.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType v1 = openvdb::v3_1::tools::PointSampler::sample(grid.get()->tree(), openvdb::Vec3R(1.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType v2 = openvdb::v3_1::tools::PointSampler::sample(grid.get()->tree(), openvdb::Vec3R(2.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType v3 = openvdb::v3_1::tools::PointSampler::sample(grid.get()->tree(), openvdb::Vec3R(0.0, 0.0, 1.0));
    cout << "Index Samples at 0, 1, 1 = " << v0 << endl;
    cout << "Index Samples at 1, 1, 1 = " << v1 << endl;
    cout << "Index Samples at 2, 1, 1 = " << v2 << endl;
    cout << "Index Samples at 0, 0, 0 = " << v3 << endl;
}

void worldSampler() {
    // Instantiate the GridSampler template on the grid type and on a box sampler
    // for thread-safe but uncached trilinear interpolation.
    openvdb::v3_1::tools::GridSampler<openvdb::FloatGrid, openvdb::v3_1::tools::PointSampler> sampler(*grid.get());
    // Compute the value of the grid at fractional coordinates in index space.
    openvdb::FloatGrid::ValueType i0 = sampler.isSample(openvdb::Vec3R(0.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType i1 = sampler.isSample(openvdb::Vec3R(1.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType i2 = sampler.isSample(openvdb::Vec3R(2.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType i3 = sampler.isSample(openvdb::Vec3R(0.0, 0.0, 0.0));
    cout << "Index Samples at 0, 1, 1 = " << i0 << endl;
    cout << "Index Samples at 1, 1, 1 = " << i1 << endl;
    cout << "Index Samples at 2, 1, 1 = " << i2 << endl;
    cout << "Index Samples at 0, 0, 0 = " << i3 << endl;

    // Compute the value of the grid at a location in world space.
    openvdb::FloatGrid::ValueType worldValue0 = sampler.wsSample(openvdb::Vec3R(0.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType worldValue1 = sampler.wsSample(openvdb::Vec3R(1.0, 1.0, 1.0));
    openvdb::FloatGrid::ValueType worldValue2 = sampler.wsSample(openvdb::Vec3R(2.0, 2.0, 2.0));
    openvdb::FloatGrid::ValueType worldValue3 = sampler.wsSample(openvdb::Vec3R(0.25, 0.25, 0.25));
    cout << "World Samples at 0, 1, 1 = " << worldValue0 << endl;
    cout << "World Samples at 1, 1, 1 = " << worldValue1 << endl;
    cout << "World Samples at 2, 1, 1 = " << worldValue2 << endl;
    cout << "World Samples at 0, 0, 0 = " << worldValue3 << endl;
}

void simpleTest() {
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
    grid = openvdb::FloatGrid::create(0.0);
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

    // Define a coordinate with large signed indices.
    for (int i = 0; i < 3; ++i) {
        for (double j = 0; j < 3; ++j) {
            for (double k = 0; k < 3; ++k) {
                openvdb::Coord ijk(i + .5, j + .5, k + .5);
                bool isMiddle = j == 1 && k == 1;
                accessor.setValue(ijk, isMiddle ? 1.0 : 0.0);
            }
        }
    }

    //printCube(accessor);
    indexSampler();
    cout << "\n\n";
    grid.get()->setTransform(openvdb::v3_1::math::Transform::createLinearTransform(0.25));
    worldSampler();

    // Choose fractional coordinates in index space.


    //    // Verify that the voxel value at (1000, 200000000, -30000000) is
    //    // the background value, 0.
    //    std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
    //    // Set the voxel value at (1000, 200000000, -30000000) to 2.
    //    accessor.setValue(xyz, 2.0);
    //    // Set the voxels at the two extremes of the available coordinate space.
    //    // For 32-bit signed coordinates these are (-2147483648, -2147483648, -2147483648)
    //    // and (2147483647, 2147483647, 2147483647).
    //    accessor.setValue(openvdb::Coord::min(), 3.0f);
    //    accessor.setValue(openvdb::Coord::max(), 4.0f);
    //    std::cout << "Testing sequential access:" << std::endl;
    //    // Print all active ("on") voxels by means of an iterator.
    //    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
    //        std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
    //    }
}

int main(int argc, const char** argv) {
    //simpleTest();

    string inputFileVdb;
    if (argc >= 2) {
        inputFileVdb = argv[1];
        cout << "Input file: " << inputFileVdb << endl;
    }

    OpenVdbReader vdbReader(inputFileVdb);
    vdbReader.initialize();
    auto grid = vdbReader.getHairDensityGrid();

    Point3 bbMin, bbMax;
    vdbReader.getBoundingBox(&bbMin, &bbMax);
    double voxelSize = vdbReader.getVoxelSize();

    double width = bbMax.x - bbMin.x;
    double height = bbMax.y - bbMin.y;
    double depth = bbMax.z - bbMin.z;
    cout << "width: " << width << ", height: " << height << ", depth: " << depth << endl;

    Point3 leftFrom(bbMin.x - voxelSize, bbMin.y + height / 2.0, bbMin.z + depth / 2.0);
    Point3 leftTo(leftFrom.x + width + 2.0 * voxelSize, leftFrom.y, leftFrom.z);

    Point3 bottomFrom(bbMin.x + width / 2.0, bbMin.y - voxelSize, bbMin.z + depth / 2.0);
    Point3 bottomTo(bottomFrom.x, bottomFrom.y + height + 2.0 * voxelSize, bottomFrom.z);

    Point3 nearFrom(bbMin.x + width / 2.0, bbMin.y + height / 2.0, bbMin.z - voxelSize);
    Point3 nearTo(nearFrom.x, nearFrom.y, nearFrom.z + depth + 2.0 * voxelSize);

    double interpolateX = vdbReader.interpolate(leftFrom, leftTo, 10000);
    cout << "Interpolate -x ---> +x: " << interpolateX << " hair strands" << endl << endl;

    double interpolateY = vdbReader.interpolate(bottomFrom, bottomTo, 10000);
    cout << "Interpolate -y ---> +y: " << interpolateY << " hair strands" << endl << endl;

    double interpolateZ = vdbReader.interpolate(nearFrom, nearTo, 10000);
    cout << "Interpolate -z ---> +z: " << interpolateZ << " hair strands" << endl << endl;


    return 0;
}