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

#ifdef __APPLE__
using namespace openvdb::v6_0::math;
#else
using namespace openvdb::v3_1::math;
#endif

const int DEFAULT_SAMPLES_PER_SEGMENT = 10;
const double DEFAULT_VOXEL_SIZE = 1.0;

// params of writevdb
double voxelSize = DEFAULT_VOXEL_SIZE;
int samplesPerSegment = DEFAULT_SAMPLES_PER_SEGMENT;

// stats
double max_sample_length = 0.0;
Point3 bbMin, bbMax;

double computeDistance(const Point3& p1, const Point3& p2) {
    double x = (p2.x - p1.x);
    double y = (p2.y - p1.y);
    double z = (p2.z - p1.z);

    return sqrt(x * x + y * y + z * z);
}

double getStepSize(int sampleSize) {
    // step size is the increment to walk through a segment
    return 1.0 / static_cast<double> (sampleSize - 1);
}

void sampleSegment(openvdb::FloatGrid::Accessor& accessor, const BezierSpline& spline, int segment) {
    if (samplesPerSegment < 2) {
        cout << "A segment must be sampled at least 2 times, terminating application..." << endl;
        exit(1);
    }

    double stepSize = getStepSize(samplesPerSegment);

    double prevDistance = 0.0;
    Point3 P = spline.sampleSegment(segment, 0.0);
    Point3 nextP;

    for (int i = 0; i < samplesPerSegment; ++i) {
        double distance = 0.0;
        if (i + 1 < samplesPerSegment) {
            nextP = spline.sampleSegment(segment, (i + 1) * stepSize);
            distance = computeDistance(P, nextP);

            bbMin.x = std::min(bbMin.x, nextP.x);
            bbMin.y = std::min(bbMin.y, nextP.y);
            bbMin.z = std::min(bbMin.z, nextP.z);

            bbMax.x = std::max(bbMax.x, nextP.x);
            bbMax.y = std::max(bbMax.y, nextP.y);
            bbMax.z = std::max(bbMax.z, nextP.z);
        }

        double lengthSample = 0.5 * (prevDistance + distance);
        double value = lengthSample / (voxelSize);

        if (lengthSample > max_sample_length) {
            max_sample_length = lengthSample;
        }

        // write to voxel grid
        openvdb::Coord xyz(P.x / voxelSize, P.y / voxelSize, P.z / voxelSize);
        accessor.setValue(xyz, accessor.getValue(xyz) + value);

        // store for next iteration
        prevDistance = distance;
        P = nextP;
    }

}

void write(openvdb::FloatGrid::Accessor& accessor, const Hair& hair) {
    unsigned int fiberCount = hair.fibers.size();
    unsigned int fiberIndex = 0;

    // set bounding box
    if (hair.fibers.size() > 0) {
        bbMin = bbMax = hair.fibers[0].curve.getControlPoints()[0];
    }

    for (const auto& fiber : hair.fibers) {
        double percentageCompleted = 100.0 * static_cast<double> (fiberIndex) / fiberCount;
        cout << "\r" << static_cast<int> (std::round(percentageCompleted)) << "%: " << fiberIndex << " / " << fiberCount << " hair fibers completed      ";

        for (int segment = 0; segment < fiber.curve.getSegmentCount(); ++segment) {

            sampleSegment(accessor, fiber.curve, segment);
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

    // assign a transform to define a voxel size (voxel size works inverse with transform)
    grid->setTransform(Transform::createLinearTransform(voxelSize));
    grid->insertMeta("VoxelSize", openvdb::FloatMetadata(voxelSize));


    // hair model is read without transformations applied to it (e.g. local space)
    // TODO: check if this can be set to true, and then dont have to self divide by voxel size
    grid->setIsInWorldSpace(false);

    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    write(accessor, hair);

    std::stringstream bbMinStream, bbMaxStream;
    bbMinStream << bbMin.x << " " << bbMin.y << " " << bbMin.z;
    bbMaxStream << bbMax.x << " " << bbMax.y << " " << bbMax.z;

    grid->insertMeta("BoundingBoxMin", openvdb::StringMetadata(bbMinStream.str()));
    grid->insertMeta("BoundingBoxMax", openvdb::StringMetadata(bbMaxStream.str()));

    grid->insertMeta("Maximum Sample Length", openvdb::FloatMetadata(max_sample_length));
    grid->insertMeta("Updated last", openvdb::Int32Metadata(time(0)));

    cout << "Longest distance between two samples points is: " << max_sample_length << endl;


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
