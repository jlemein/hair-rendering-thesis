#include "hairsimplify.h"

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <sstream>

HairSimplify::HairSimplify(const std::string& inputFileName) : mInputFileName(inputFileName) {
    analyzeHairFile();
    mCurveWidth = -1.0f;
    mCurveType = "flat";
}

void HairSimplify::analyzeHairFile() {
    std::ifstream infile;
    infile.open(mInputFileName.c_str());
    if (infile.fail()) {
        std::cout << "Failed to open input file, aborting ..." << std::endl;
        return;
    }

    std::string line;
    mHairCount = 0;
    while (std::getline(infile, line)) {
        mHairCount++;
    }

    std::cout << "File contains " << mHairCount << " lines" << std::endl;
    infile.close();
}

void HairSimplify::setCurveWidth(float curveWidth) {
    this->mCurveWidth = curveWidth;
}

void HairSimplify::setCurveType(std::string type) {
    this->mCurveType = type;
}

void HairSimplify::reduceToPercentage(const std::string& outputFileName, double percentage, bool sampleRandom) {
    mHairCount = static_cast<int> (mHairCount * (percentage / 100.0));
    // if (sampleRandom) {
    //     sampleHairsRandomly(inputFileName, outputFileName, numberHairs);
    // } else {
    sampleFromBeginning(outputFileName);
    // }
}

std::string replaceCurveWidths(std::string line, float width0, float width1, std::string curveType) {
    std::string newString = "";
    int stringTypeIndex = line.find("\"string type\"");
    int pointTypeIndex = line.find("\"point P\"");


    int index = line.find("width0");
    if (index == std::string::npos) {
        return line;
    } else {
        std::ostringstream os;
        //std::cout << "Found at index: " << index << ": " << line.substr(0, index) << std::endl;

        os << line.substr(0, stringTypeIndex) << "\"string type\" [ " << "\"" << curveType << "\"" << " ] " << line.substr(pointTypeIndex, index - pointTypeIndex) << "width0\" [ " << width0 << " ] \"float width1\" [ " << width1 << " ]";
        return os.str();
    }

}

void HairSimplify::sampleFromBeginning(const std::string& outputFileName) {
    std::cout << "Sample from beginning up to " << mHairCount << " hair strands" << std::endl;

    std::ifstream inFile(mInputFileName.c_str());
    if (inFile.fail()) {
        std::cout << "Cannot read input file " << mInputFileName << std::endl;
        return;
    }

    std::ofstream outFile(outputFileName.c_str());
    if (outFile.fail()) {
        inFile.close();
        std::cout << "Cannot open output file for writing ('" << outputFileName << "')" << std::endl;
        return;
    }

    int lineCount = 0;
    std::string line;
    while (getline(inFile, line) && lineCount < mHairCount) {
        float curve0 = (rand() / static_cast<float> (RAND_MAX)) * 0.1f + 0.05f;
        float curve1 = (rand() / static_cast<float> (RAND_MAX)) * 0.1f + 0.05f;

        if (this->mCurveWidth >= 0.0f) {
            curve0 = mCurveWidth;
            curve1 = mCurveWidth;
        }
        outFile << replaceCurveWidths(line, curve0, curve1, mCurveType) << std::endl;
        lineCount++;
    }
}

void HairSimplify::sampleRandomly() {
    std::cout << "Not supported at the moment" << std::endl;
}
