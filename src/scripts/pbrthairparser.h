#ifndef PBRTHAIRPARSER_H
#define PBRTHAIRPARSER_H

#include <iostream>

struct Hair;

class PbrtHairParser
{
private:
    static bool nextWord(std::istream& istream, std::string& outString);
    static bool nextValues(std::istream& inputStream, std::string& outString);
    PbrtHairParser();

public:

    static void parseFile(std::istream& istream, Hair& outHair);
};

#endif // PBRTHAIRPARSER_H
