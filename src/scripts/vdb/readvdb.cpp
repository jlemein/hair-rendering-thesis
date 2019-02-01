/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <openvdb/openvdb.h>
#include <iostream>
#include <fstream>
#include <string>

#include "InputOutputUtil.h"
#include "../hairstruct.h"
#include "OpenVdbReader.h"
using namespace std;

int main(int argc, const char** argv) {
    cout << "== ReadVDB ==" << endl;

    // 1. Ask user to enter file name and open file
    string inputVdbFileName = argc > 1 ? argv[1] : "mygrid.vdb";
    OpenVdbReader vdbReader(inputVdbFileName);

    try {
        vdbReader.initialize();
        vdbReader.printVdbInfo();
    } catch (openvdb::IoError err) {
        cout << "Cannot read voxel grid specified in " << inputVdbFileName << endl;
    }

    return 0;
}