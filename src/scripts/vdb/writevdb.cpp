/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <openvdb/openvdb.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <exception>

#include "InputOutputUtil.h"
#include "../hairstruct.h"
#include "../bezier.h"
using namespace std;

//FIXME: cmake for one version of openvdb
#ifdef __APPLE__
using namespace openvdb::v6_0::math;
#else
using namespace openvdb::v3_1::math;
#endif

const int DEFAULT_SAMPLES_PER_SEGMENT = 10;
const double DEFAULT_VOXEL_SIZE = 1.0;

// ----------------------------------------------------------
// INPUT PARAMS (for application writevdb)

// size of a voxel cell in units of hair model file
double voxelSize = DEFAULT_VOXEL_SIZE;

// samples to be taken per curve segment. The larger the more accurate the voxel grid becomes.
int samplesPerSegment = DEFAULT_SAMPLES_PER_SEGMENT;
// ----------------------------------------------------------

//
// this is the maximum sample length, stored after hair file is written to voxel grid.
// this length related to the voxel size gives a broad indication of the accuracy
// of the generated voxel grid
double maxSampleLength = 0.0;

// these are the bounding box positions of the hair model.
// after hair model is written to voxel grid, these reflect the bounding box positions
// of the hair model.
Point3 bbMin, bbMax;

/**
 * Computes the step size to be used when sampling the curve
 * @param sampleSize should be >= 2, the amount of samples to be taken for a sampling a curve segment.
 * @return the step size in domain (0, 1] to be used when sampling through the hair fiber.
 */
double getStepSize(int sampleSize) {
    if (samplesPerSegment < 2) {
        cout << "A segment must be sampled at least 2 times, terminating application..." << endl;
        exit(1);
    }

    // step size is the increment to walk through a segment
    return 1.0 / static_cast<double> (sampleSize - 1);
}

void updateBoundingBox(const Point3& p) {
    bbMin.x = std::min(bbMin.x, p.x);
    bbMin.y = std::min(bbMin.y, p.y);
    bbMin.z = std::min(bbMin.z, p.z);

    bbMax.x = std::max(bbMax.x, p.x);
    bbMax.y = std::max(bbMax.y, p.y);
    bbMax.z = std::max(bbMax.z, p.z);
}

/**
 *
 * @param grid
 * @param accessor
 * @param spline
 * @param segment
 */
void sampleSegment(openvdb::FloatGrid::Ptr grid, openvdb::FloatGrid::Accessor& accessor, const BezierSpline& spline, int segment) {
    const double stepSize = getStepSize(samplesPerSegment);

    double prevDistance = 0.0;

    // Sample the starting point of the segment
    Point3 P = spline.sampleSegment(segment, 0.0);
    Point3 nextP;

    for (int i = 0; i < samplesPerSegment; ++i) {
        double distance = 0.0;

        // if not yet at the last control point
        if (i + 1 < samplesPerSegment) {
            nextP = spline.sampleSegment(segment, (i + 1) * stepSize);
            distance = Point3::DistanceBetween(P, nextP);

            updateBoundingBox(nextP);
        }

        double lengthSample = 0.5 * (prevDistance + distance);
        double value = lengthSample / (voxelSize);
        if (value < 0.0) {
            cout << "value: " << value << " (distance: " << lengthSample << ", voxelsize: " << voxelSize << ")" << endl;
        }

        if (lengthSample > maxSampleLength) {
            maxSampleLength = lengthSample;
        }

        // write to voxel grid
        openvdb::Vec3d indexSpace = grid->transform().worldToIndex(openvdb::Vec3d(P.x, P.y, P.z));
        openvdb::Coord ijk(indexSpace.x(), indexSpace.y(), indexSpace.z());

        accessor.setValue(ijk, accessor.getValue(ijk) + value);
        openvdb::FloatGrid::ValueType stored = accessor.getValue(ijk);
        if (stored < 0.0) {
            cout << "Stored value: " << stored << endl;
        }

        // store for next iteration
        prevDistance = distance;
        P = nextP;
    }

}

void write(openvdb::FloatGrid::Ptr grid, const Hair& hair) {
    unsigned int fiberCount = hair.fibers.size();
    unsigned int fiberIndex = 0;


    // set bounding box
    if (hair.fibers.size() > 0) {
        bbMin = bbMax = hair.fibers[0].curve.getControlPoints()[0];
    }

    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    for (const auto& fiber : hair.fibers) {
        double percentageCompleted = 100.0 * static_cast<double> (fiberIndex) / fiberCount;
        cout << "\r" << static_cast<int> (std::round(percentageCompleted)) << "%: " << fiberIndex << " / " << fiberCount << " hair fibers completed      ";

        for (int segment = 0; segment < fiber.curve.getSegmentCount(); ++segment) {

            sampleSegment(grid, accessor, fiber.curve, segment);
        }
        ++fiberIndex;
    }

    cout << "\r" << "100%: " << fiberIndex << " / " << fiberCount << " hair fibers completed        \n";
}

void showProgramInfo() {

    cout << "This program writes PBRT hair models to a voxel grid (using OpenVDB)\n"
            << "The hair models are described in a scene file format, usually ending in *.pbrt\n\n"
            << "This program expects input arguments to be run successfully:\n"
            << "\twritevdb <0> [1] [2] \n"
            << "\t [0] (string) required - input path to hair file\n"
            << "\t [1] (string) optional - output file path of vdb file (default: mygrid.vdb)\n"
            << "\t [2] (int)    optional - samples per curve/hair fiber segment (default: " << DEFAULT_SAMPLES_PER_SEGMENT << ")"
            << "\t [3] (double) optional - voxel size (default: " << DEFAULT_VOXEL_SIZE << ")\n"
            << endl;
}

int main(int argc, const char** argv) {
    cout << "=== Write VDB ===" << endl;

    if (argc <= 1) {
        // ask user for input for hair file
        cout << "ERROR: No input hair file provided, exiting\n\n";
        showProgramInfo();
        return 1;
    }

    string inputHairFile = argv[1];
    string outputVdbFile = argc > 2 ? argv[2] : "mygrid.vdb";

    samplesPerSegment = DEFAULT_SAMPLES_PER_SEGMENT;
    if (argc > 3) {
        try {
            std::string samplesPerSegmentStr = argv[3];
            samplesPerSegment = std::stoi(samplesPerSegmentStr);
            if (samplesPerSegment <= 1) {
                cout << "Samples per segment invalid: must be positive and larger than 2. Setting to default sample size of " << DEFAULT_SAMPLES_PER_SEGMENT << endl;
                samplesPerSegment = DEFAULT_SAMPLES_PER_SEGMENT;
            }
        } catch (exception& e) {
            cout << "Samples per segment must be an integer. Using default samples per segment of " << DEFAULT_SAMPLES_PER_SEGMENT << endl;
        }
    }

    voxelSize = DEFAULT_VOXEL_SIZE;
    if (argc > 4) {
        try {
            std::string voxelSizeStr = argv[4];
            voxelSize = std::stod(voxelSizeStr);
            if (voxelSize <= 0.0) {
                cout << "Voxel size invalid: must be positive value. Setting to default voxel size of " << DEFAULT_VOXEL_SIZE << endl;
                voxelSize = DEFAULT_VOXEL_SIZE;
            }
        } catch (exception& e) {
            cout << "Voxel size must be a double. Using default voxel size of " << DEFAULT_VOXEL_SIZE << endl;
        }
    }

    cout << "Writing voxel grid to " << outputVdbFile << endl;

    // 1. Ask user to enter file name and open file
    ifstream hairFile;
    if (!InputOutputUtil::OpenFile(hairFile, inputHairFile)) {
        cout << "Failed to open hair input file, terminating application..." << endl;
        return -1;
    }

    cout << "Parsing hair file ..." << endl;
    Hair hair;
    hairFile >> hair;

    openvdb::initialize();

    // create voxel grid with background value 0
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/0.0);
    grid->setName("HairDensityGrid");
    grid->insertMeta("Author", openvdb::StringMetadata("Jeffrey Lemein"));

    // assign a transform to define a voxel size
    grid->setIsInWorldSpace(true);
    grid->setTransform(Transform::createLinearTransform(voxelSize));
    grid->insertMeta("VoxelSize", openvdb::FloatMetadata(voxelSize));

    write(grid, hair);

    std::stringstream bbMinStream, bbMaxStream;
    bbMinStream << bbMin.x << " " << bbMin.y << " " << bbMin.z;
    bbMaxStream << bbMax.x << " " << bbMax.y << " " << bbMax.z;

    grid->insertMeta("BoundingBoxMin", openvdb::StringMetadata(bbMinStream.str()));
    grid->insertMeta("BoundingBoxMax", openvdb::StringMetadata(bbMaxStream.str()));

    grid->insertMeta("Maximum Sample Length", openvdb::FloatMetadata(maxSampleLength));
    grid->insertMeta("Updated last", openvdb::Int32Metadata(time(0)));

    cout << "Longest distance between two samples points is: " << maxSampleLength << endl;


    // Create a VDB file object.
    cout << "Writing VDB grid to " << outputVdbFile << endl;
    openvdb::io::File file(outputVdbFile);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();
    cout << "[done]" << endl;

    return 0;
}
