#include "hairsimplify.h"

#include <string>
#include <fstream>
#include <iostream>

HairSimplify::HairSimplify(const std::string& inputFileName) : mInputFileName(inputFileName) {
    analyzeHairFile();
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

void HairSimplify::reduceToPercentage(const std::string& outputFileName, double percentage, bool sampleRandom) {
    mHairCount = static_cast<int>(mHairCount * (percentage / 100.0));
    // if (sampleRandom) {
    //     sampleHairsRandomly(inputFileName, outputFileName, numberHairs);
    // } else {
        sampleFromBeginning(outputFileName);
    // }
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
    outFile << line << std::endl;
    lineCount++;
  }
}
 
 void HairSimplify::sampleRandomly() {
   std::cout << "Not supported at the moment" << std::endl;
 }
