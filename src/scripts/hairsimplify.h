#ifndef HAIR_SIMPLIFY_H
#define HAIR_SIMPLIFY_H

#include <string>

class HairSimplify {
    std::string mInputFileName;
    int mHairCount;

    void analyzeHairFile();
    void sampleFromBeginning(const std::string& outputFileName);
    void sampleRandomly();

public:
    HairSimplify(const std::string& inputFileName);
    void reduceToPercentage(const std::string& outputFileName, double percentage, bool sampleRandom = false);
};

#endif // HAIR_SIMPLIFY_H
