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

#include "InputOutputUtil.h"
#include "../hairstruct.h"
#include "../bezier.h"
using namespace std;

const int SAMPLES_PER_SEGMENT = 10;
const double VOXEL_SIZE = 1.0;

const char* ARGUMENT_NAMES[] = {
    "Run directory",
    "Hair input file"
};

void writeSpline(const BezierSpline& spline) {

}

double computeDistance(const Point3& p1, const Point3& p2) {
    double x = (p2.x - p1.x);
    double y = (p2.y - p1.y);
    double z = (p2.z - p1.z);

    double d = sqrt(x * x + y * y + z * z);
    cout << "distance between " << p1 << " and " << p2 << " = " << d << endl;
    return d;
}

double getStepSize(int sampleSize) {
    // step size is the increment to walk through a segment
    double stepSize = 1.0 / static_cast<double> (sampleSize - 1);
    cout << "Step size for sample size " << sampleSize << ": " << stepSize << endl;
    return stepSize;
}

void sampleSegment(const BezierSpline& spline, int segment) {
    if (SAMPLES_PER_SEGMENT < 2) {
        cout << "A segment must be sampled at least 2 times, terminating application..." << endl;
        exit(1);
    }

    cout << "Segment defined by: \n";
    for (auto& pt : spline.getControlPoints()) {
        cout << pt << endl;
    }
    cout << "\n\n";

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

        double value = 0.5 * (prevDistance + distance) / VOXEL_SIZE;

        // write to voxel grid
        cout << "Write( " << P << ") = " << value << endl;

        // store for next iteration
        prevDistance = distance;
        P = nextP;
    }

}

int main(int argc, const char** argv) {
    cout << "== Write VDB ==" << endl;

    for (int i = 0; i < argc; ++i) {
        cout << ARGUMENT_NAMES[i] << ": " << argv[i] << endl;
    }
    cout << "------------------------------------\n\n";

    // 1. Ask user to enter file name and open file
    ifstream hairFile;
    if (!InputOutputUtil::OpenFile(hairFile, argc, argv)) {
        cout << "Failed to open hair input file, terminating application..." << endl;
        return -1;
    }

    //Hair hair;
    //hairFile >> hair;

    //for (auto& fiber : hair.fibers) {

    BezierSpline spline(
            Point3(-11, 49, -14),
            Point3(-11, 49, -14),
            Point3(-11, 49, -14),
            Point3(-11, 49, -14));
    //const BezierSpline& spline = hair.fibers[0].curve;

    for (int segment = 0; segment < spline.getSegmentCount(); ++segment) {
        sampleSegment(spline, segment);
    }



    //}

    //cout << "Loaded hair model with " << hair.fibers.size() << " fibers" << endl;

    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/2.0);
    grid->setName("MySpecialGrid");

    // Create a VDB file object.
    openvdb::io::File file("mygrids.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();

    return 0;
}