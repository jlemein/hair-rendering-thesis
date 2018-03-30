#include <iostream>
#include "sayhello.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "hairsimplify.h"

using namespace std;

enum SamplingMode { RANDOM, FROM_BEGINNING };

string inputFileName, outputFileName;
SamplingMode samplingMode = RANDOM;
int maxHairCount;
bool isHairCountAbsolute;

int hairCountOriginal;


bool isSame(const string word, const string key) {
  cout << "comparing '" << word << "' with '" << key << "'";
  int result = word.compare(key);
  //cout << " result: " << result << endl;
  return result == 0;
}

void parseCommandLineArguments(int argc, char** argv) {
  // parsing parameters
  for (int i = 1; i+1 < argc; i+=2) {
    if (isSame(string(argv[i]), string("--input"))) {
        inputFileName = argv[i+1];
        continue;
    }
    if (isSame(argv[i], "--output")) {
        outputFileName = argv[i+1];
        continue;
    }
    if (isSame(argv[i], "--sampling")) {
        string strategyName = argv[i+1];
        if (strategyName == "random") {
            samplingMode = RANDOM;
        }
        else if (isSame(strategyName, "frombegin")) {
            samplingMode = FROM_BEGINNING;
        } else {
            cout << "Strategy unknown, please provide either 'random' or 'frombegin'";
            return;
        }
        continue;
    }
    if (isSame(argv[i], "--hairCount")) {
        string hairCountStr = argv[i+1];
        if (*hairCountStr.rbegin() == '%') {
	  string hairCountInput = hairCountStr.substr(0, hairCountStr.length()-1);
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
       << "max hair count: " << maxHairCount << " " << (isHairCountAbsolute ? "strands" : "procent of input file") << endl;

}



int main(int argc, char** argv) {
  parseCommandLineArguments(argc, argv);

  if (!isHairCountAbsolute) {
    int percentage = maxHairCount;
    HairSimplify hair(inputFileName);
    hair.reduceToPercentage(outputFileName, percentage, samplingMode == RANDOM);
  }
    
  return 0;
}
