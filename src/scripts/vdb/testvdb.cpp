/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <openvdb/openvdb.h>
#include <iostream>
#include <fstream>
#include <string>
#include "../hairstruct.h"
#include "OpenVdbReader.h"
using namespace std;

int main(int argc, const char** argv) {
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

    Point3 leftFrom(bbMin.x - voxelSize, bbMin.y + height / 2.0, bbMin.z + depth / 2.0);
    Point3 leftTo(leftFrom.x + width + 2.0 * voxelSize, leftFrom.y, leftFrom.z);

    Point3 bottomFrom(bbMin.x + width / 2.0, bbMin.y - voxelSize, bbMin.z + depth / 2.0);
    Point3 bottomTo(bottomFrom.x, bottomFrom.y + height + 2.0 * voxelSize, bottomFrom.z);

    Point3 nearFrom(bbMin.x + width / 2.0, bbMin.y + height / 2.0, bbMin.z - voxelSize);
    Point3 nearTo(nearFrom.x, nearFrom.y, nearFrom.z + depth + 2.0 * voxelSize);

    double interpolateX = vdbReader.interpolate(leftFrom, leftTo);
    cout << "Interpolate -x ---> +x: " << interpolateX << " hair strands" << endl;

    double interpolateY = vdbReader.interpolate(bottomFrom, bottomTo);
    cout << "Interpolate -y ---> +y: " << interpolateY << " hair strands" << endl;

    double interpolateZ = vdbReader.interpolate(nearFrom, nearTo);
    cout << "Interpolate -z ---> +z: " << interpolateZ << " hair strands" << endl;


    return 0;
}