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

#include "InputOutputUtil.h"
#include "../hairstruct.h"
#include "../bezier.h"
using namespace std;

const int SAMPLES_PER_SEGMENT = 10;
const double VOXEL_SIZE = 0.25;

double max_sample_length = 0.0;

const char* ARGUMENT_NAMES[] = {
    "Run directory",
    "Hair input file",
    "Output openvdb file"
};

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
    if (SAMPLES_PER_SEGMENT < 2) {
        cout << "A segment must be sampled at least 2 times, terminating application..." << endl;
        exit(1);
    }

    double stepSize = getStepSize(SAMPLES_PER_SEGMENT);

    double prevDistance = 0.0;
    Point3 P = spline.sampleSegment(segment, 0.0);
    Point3 nextP;

    for (int i = 0; i < SAMPLES_PER_SEGMENT; ++i) {
        double distance = 0.0;
        if (i + 1 < SAMPLES_PER_SEGMENT) {
            nextP = spline.sampleSegment(segment, (i + 1) * stepSize);
            distance = computeDistance(P, nextP);

        }

        double lengthSample = 0.5 * (prevDistance + distance);
        double value = lengthSample / (VOXEL_SIZE);

        if (lengthSample > max_sample_length) {
            max_sample_length = lengthSample;
        }

        // write to voxel grid
        //cout << "Write( " << P << ") = " << value << endl;
        openvdb::Coord xyz(P.x, P.y, P.z);
        accessor.setValue(xyz, accessor.getValue(xyz) + value);

        // store for next iteration
        prevDistance = distance;
        P = nextP;
    }

}

void write(openvdb::FloatGrid::Accessor& accessor, const Hair& hair) {
    unsigned int fiberCount = hair.fibers.size();
    unsigned int fiberIndex = 0;

    for (const auto& fiber : hair.fibers) {
        double percentageCompleted = 100.0 * static_cast<double> (fiberIndex) / fiberCount;
        cout << "\r" << percentageCompleted << " %: " << fiberIndex << " / " << fiberCount << " hair fibers completed      ";

        for (int segment = 0; segment < fiber.curve.getSegmentCount(); ++segment) {
            sampleSegment(accessor, fiber.curve, segment);
        }
        ++fiberIndex;
    }

    cout << "\r" << "100.0%: " << fiberIndex << " / " << fiberCount << " hair fibers completed        \n";
}

int main(int argc, const char** argv) {
    cout << "== Write VDB ==" << endl;

    if (argc <= 1) {
        // ask user for input for hair file
        cout << "No hair file provided, exiting\n";
        return 1;
    }

    string inputHairFile = argv[1];
    string outputVdbFile = argc > 2 ? argv[2] : "mygrid.vdb";

    // 1. Ask user to enter file name and open file
    ifstream hairFile;
    if (!InputOutputUtil::OpenFile(hairFile, inputHairFile)) {
        cout << "Failed to open hair input file, terminating application..." << endl;
        return -1;
    }
    std::cout << "hair : " << (hairFile.fail() ? "failed" : "success") << endl;

    Hair hair;
    hairFile >> hair;
    cout << "Loaded hair model with " << hair.fibers.size() << " fibers" << endl;

    openvdb::initialize();

    // create voxel grid with background value 0
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/0.0);
    grid->setName("HairDensityGrid");
    grid->insertMeta("Author", openvdb::StringMetadata("Jeffrey Lemein"));

    // assign a transform to define a voxel size (for voxelSize 0.25 -> 1 unit in the hair model
    grid->setTransform(openvdb::v3_1::math::Transform::createLinearTransform(VOXEL_SIZE));
    grid->insertMeta("VoxelSize", openvdb::DoubleMetadata(VOXEL_SIZE));

    // hair model is read without transformations applied to it (e.g. local space)
    grid->setIsInWorldSpace(false);

    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    write(accessor, hair);

    grid->insertMeta("Maximum Sample Length", openvdb::FloatMetadata(max_sample_length));
    grid->insertMeta("Updated last", openvdb::Int32Metadata(time(0)));

    // Create a VDB file object.
    cout << "Writing VDB grid to " << outputVdbFile << endl;
    openvdb::io::File file(outputVdbFile);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();

    return 0;
}