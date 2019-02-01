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
using namespace std;

int main(int argc, const char** argv) {
    cout << "== ReadVDB ==" << endl;

    // 1. Ask user to enter file name and open file
    ifstream inputFile;
    if (!InputOutputUtil::OpenFile(inputFile, argc, argv)) {
        cout << "Failed to open input file, terminating application..." << endl;
        return -1;
    }

    openvdb::initialize();

    string line;
    while (getline(inputFile, line)) {
        cout << line << '\n';
    }

    return 0;
}