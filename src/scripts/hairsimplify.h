#ifndef HAIR_SIMPLIFY_H
#define HAIR_SIMPLIFY_H

#include <string>

class HairSimplify {
    std::string mInputFileName;
    int mHairCount;
    float mCurveWidth;
    std::string mCurveType;

    void analyzeHairFile();
    void sampleFromBeginning(const std::string& outputFileName);
    void sampleRandomly();
    void sampleByAxis(const std::string& outputFileName);

public:
    HairSimplify(const std::string& inputFileName);
    void reduceToPercentage(const std::string& outputFileName, double percentage, bool sampleRandom = false);
    void reduceByAxis(const std::string& outputFileName);
    
    void setCurveWidth(float curveWidth);
    void setCurveType(std::string type);
};

#endif // HAIR_SIMPLIFY_H
