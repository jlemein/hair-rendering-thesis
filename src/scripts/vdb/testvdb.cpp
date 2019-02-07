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

    Point3 p1, p2;
    vdbReader.getBoundingBox(&p1, &p2);
    double voxelSize = vdbReader.getVoxelSize();

    double width = p2.x - p1.x;
    double height = p2.y - p1.y;
    double depth = p2.z - p1.z;

    Point3 leftFrom(p1.x - voxelSize, p1.y + height / 2.0, p1.z + depth / 2.0);
    Point3 leftTo(leftFrom.x + width + 2.0 * voxelSize, leftFrom.y, leftFrom.z);

    double interpolateX = vdbReader.interpolate(leftFrom, leftTo);

    cout << "Interpolate -x ---> +x: " << interpolateX << " hair strands" << endl;

    return 0;
}