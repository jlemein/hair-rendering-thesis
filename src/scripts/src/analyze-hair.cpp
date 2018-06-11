#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>


#include "hairstruct.h"



int main(int argc, char** argv) {
    if (argc <= 1) {
        std::cout << "No input file specified" << std::endl;
        return -1;
    }

    const char* inputFilename = argv[1];
    std::ifstream infile(inputFilename);
    if (infile.fail()) {
        std::cout << "Cannot open file '" << inputFilename << "'" << std::endl;
        return -1;
    }

    Hair hair;
    std::cout << "Parsing hair model '" << inputFilename << "'" << std::endl;
    infile >> hair;
    std::cout << "Finished reading" << std::endl << std::endl;

    std::cout << "Model bounding box information\n-------------------------------" << std::endl;
    std::cout << "Scene extent: " << hair.sceneExtentMin << " x " << hair.sceneExtentMax << std::endl;
    std::cout << "Scene size: " << hair.size << std::endl;
    std::cout << "Center point: " << hair.center << std::endl << std::endl;

    std::cout << "Curve info\n-------------------------------" << std::endl;
    std::cout << "Number of curves: " << hair.curves.size() << std::endl;
    std::cout << "Range curve width0: [ " << hair.minCurveWidth0 << " ; " << hair.maxCurveWidth0 << " ], average: " << hair.avgCurveWidth0 << std::endl;
    std::cout << "Range curve width1: [ " << hair.minCurveWidth1 << " ; " << hair.maxCurveWidth1 << " ], average: " << hair.avgCurveWidth1 << std::endl << std::endl;
    if (abs(hair.avgCurveWidth0 - hair.avgCurveWidth1) < 0.05f*hair.avgCurveWidth0) {
        std::cout << "\tHair strands have circular cross section" << std::endl;
    } else {
        std::cout << "\tHair strands have elliptical cross section" << std::endl;
    }


    return 0;
}
