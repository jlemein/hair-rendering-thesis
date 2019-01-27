#include <iostream>
#include "sayhello.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "hairsimplify.h"

using namespace std;

enum SamplingMode {
    RANDOM, FROM_BEGINNING
};

string inputFileName, outputFileName;
SamplingMode samplingMode = RANDOM;
int maxHairCount;
bool isHairCountAbsolute;
float curveWidth = -1.0f;
string curveType = "flat";

int hairCountOriginal;

bool isSame(const string word, const string key) {
    cout << "comparing '" << word << "' with '" << key << "'";
    int result = word.compare(key);
    //cout << " result: " << result << endl;
    return result == 0;
}

void parseCommandLineArguments(int argc, char** argv) {
    // parsing parameters
    if (argc <= 1) {
        cout << "Use as follows: \n"
                << "simplify --input <infile> --output <outfile> [--sampling 'random'|'frombegin'] [--haircount '%'|<number_requested_hairs>] --curvetype [flat|cylinder|ribbon] --curvewidth <float or empty (to sample random)>" << endl;
        exit(0);
    }
    for (int i = 1; i + 1 < argc; i += 2) {
        if (isSame(string(argv[i]), string("--input"))) {
            inputFileName = argv[i + 1];
            continue;
        }
        if (isSame(argv[i], "--output")) {
            outputFileName = argv[i + 1];
            continue;
        }
        if (isSame(argv[i], "--curvetype")) {
            curveType = argv[i + 1];
            cout << "Curve type specified: " << curveType;
            continue;
        }
        if (isSame(argv[i], "--curvewidth")) {
            curveWidth = atof(string(argv[i + 1]).c_str());
            cout << "Curve width specified: " << curveWidth;
            continue;
        }
        if (isSame(argv[i], "--sampling")) {
            string strategyName = argv[i + 1];
            if (strategyName == "random") {
                samplingMode = RANDOM;
            } else if (isSame(strategyName, "frombegin")) {
                samplingMode = FROM_BEGINNING;
            } else {
                cout << "Strategy unknown, please provide either 'random' or 'frombegin'";
                return;
            }
            continue;
        }
        if (isSame(argv[i], "--haircount")) {
            string hairCountStr = argv[i + 1];
            if (*hairCountStr.rbegin() == '%') {
                string hairCountInput = hairCountStr.substr(0, hairCountStr.length() - 1);
                cout << "Entered hair count is " << hairCountInput << endl;

                maxHairCount = atoi(hairCountInput.c_str());
                isHairCountAbsolute = false;
            } else {
                maxHairCount = atoi(hairCountStr.c_str());
                isHairCountAbsolute = true;
            }
            continue;
        }
        cout << "Not recognized: " << argv[i] << endl;
    }

    cout << "Hair simplifier\n";
    cout << "Provided parameters: " << endl;
    cout << "input file: " << inputFileName << endl
            << "output file: " << outputFileName << endl
            << "sampling strategy: " << (samplingMode == RANDOM ? "random" : "from beginning") << endl
            << "max hair count: " << maxHairCount << " " << (isHairCountAbsolute ? "strands" : "procent of input file") << endl
            << "curve width: " << curveWidth << endl
            << "curve type: " << curveType << endl;

}

int main(int argc, char** argv) {
    parseCommandLineArguments(argc, argv);

    if (!isHairCountAbsolute) {
        int percentage = maxHairCount;
        HairSimplify hair(inputFileName);

        if (curveWidth >= 0.0f)
            hair.setCurveWidth(curveWidth);

        hair.setCurveType(curveType);
        hair.reduceToPercentage(outputFileName, percentage, samplingMode == RANDOM);
    }

    return 0;
}
