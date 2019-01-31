#include "pbrthairparser.h"

#include "hairstruct.h"
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include "bezier.h"
#include "scripts/bezier.h"

PbrtHairParser::PbrtHairParser() {

}

void PbrtHairParser::parseFile(std::istream& istream, Hair& hair) {
    hair.fibers.reserve(10000);

    std::string word;
    float minCurveWidth0 = 9999999, maxCurveWidth0 = -10000, minCurveWidth1 = 999999, maxCurveWidth1 = -10000;
    double sumCurveWidth0 = 0.0, sumCurveWidth1 = 0.0;
    bool firstRunCurves = true;

    while (PbrtHairParser::nextWord(istream, word)) {
        std::string values;

        HairFiber fiber;
        fiber.curve.getControlPoints().reserve(12);
        if (word == "Shape") {
            PbrtHairParser::nextValues(istream, fiber.type);
        } else if (word == "float width") {
            PbrtHairParser::nextValues(istream, values);
            fiber.width0 = fiber.width1 = boost::lexical_cast<float>(values);

            sumCurveWidth0 += fiber.width0;
            sumCurveWidth1 += fiber.width1;
            if (fiber.width0 < minCurveWidth0) {
                minCurveWidth0 = fiber.width0;
            }
            if (fiber.width0 > maxCurveWidth0) {
                maxCurveWidth0 = fiber.width0;
            }
            if (fiber.width1 < minCurveWidth1) {
                minCurveWidth1 = fiber.width1;
            }
            if (fiber.width1 > maxCurveWidth1) {
                maxCurveWidth1 = fiber.width1;
            }

        } else if (word == "float width0") {
            PbrtHairParser::nextValues(istream, values);
            fiber.width0 = boost::lexical_cast<float>(values);

            sumCurveWidth0 += fiber.width0;
            if (fiber.width0 < minCurveWidth0) {
                minCurveWidth0 = fiber.width0;
            }
            if (fiber.width0 > maxCurveWidth0) {
                maxCurveWidth0 = fiber.width0;
            }

        } else if (word == "float width1") {
            PbrtHairParser::nextValues(istream, values);
            fiber.width1 = boost::lexical_cast<float>(values);

            sumCurveWidth1 += fiber.width1;
            if (fiber.width1 < minCurveWidth1) {
                minCurveWidth1 = fiber.width1;
            }
            if (fiber.width1 > maxCurveWidth1) {
                maxCurveWidth1 = fiber.width1;
            }

        } else if (word == "point P") {
            std::string points;
            PbrtHairParser::nextValues(istream, points);
            std::stringstream ss(points);
            Point3 p;

            bool firstRun = true;
            while (ss >> p.x >> p.y >> p.z) {
                if (firstRun) {
                    hair.sceneExtentMax = hair.sceneExtentMin = p;
                    firstRun = false;
                }

                if (p.x < hair.sceneExtentMin.x) hair.sceneExtentMin.x = p.x;
                if (p.y < hair.sceneExtentMin.y) hair.sceneExtentMin.y = p.y;
                if (p.z < hair.sceneExtentMin.z) hair.sceneExtentMin.z = p.z;

                if (p.x > hair.sceneExtentMax.x) hair.sceneExtentMax.x = p.x;
                if (p.y > hair.sceneExtentMax.y) hair.sceneExtentMax.y = p.y;
                if (p.z > hair.sceneExtentMax.z) hair.sceneExtentMax.z = p.z;

                fiber.curve.getControlPoints().push_back(p);
            }
            hair.avgCurveWidth0 = 0;
            hair.avgCurveWidth1 = 0;
            hair.minCurveWidth0 = 0;
            hair.minCurveWidth1 = 0;
            hair.fibers.push_back(fiber);
        } else {
            PbrtHairParser::nextValues(istream, values);
        }

        hair.size.x = hair.sceneExtentMax.x - hair.sceneExtentMin.x;
        hair.size.y = hair.sceneExtentMax.y - hair.sceneExtentMin.y;
        hair.size.z = hair.sceneExtentMax.z - hair.sceneExtentMin.z;

        hair.center.x = hair.sceneExtentMin.x + hair.size.x / 2.0;
        hair.center.y = hair.sceneExtentMin.y + hair.size.y / 2.0;
        hair.center.z = hair.sceneExtentMin.z + hair.size.z / 2.0;

        hair.avgCurveWidth0 = sumCurveWidth0 / hair.fibers.size();
        hair.avgCurveWidth1 = sumCurveWidth1 / hair.fibers.size();
        hair.minCurveWidth0 = minCurveWidth0;
        hair.maxCurveWidth0 = maxCurveWidth0;
        hair.minCurveWidth1 = minCurveWidth1;
        hair.maxCurveWidth1 = maxCurveWidth1;
    }
}

bool PbrtHairParser::nextWord(std::istream& istream, std::string& outString) {
    // find beginning
    while (std::getline(istream, outString, '\"')) {
        boost::algorithm::trim(outString);
        if (outString.length() > 0) {
            return true;
        }
    }

    return false;
}

bool PbrtHairParser::nextValues(std::istream& inputStream, std::string& outString) {
    std::string line;
    inputStream >> line;

    boost::algorithm::trim(line);
    if (line.empty()) {
        return false;
    } else {
        if (line[0] == '\"' && line.length() > 1) {
            if (line[line.length() - 1] == '\"') {
                outString = line.substr(1, line.length() - 1);
            } else {
                outString = line.substr(1);
                boost::algorithm::trim(outString);
                std::string result;
                std::getline(inputStream, result, '\"');
                outString += ' ' + result;
            }
        } else if (line[0] == '[') {
            //std::cout << "READ 2 " << line << std::endl;

            if (line.length() > 1 && line[line.length() - 1] == ']') {
                outString = line.substr(1, line.length() - 1);

                boost::algorithm::trim(outString);
                if (outString.length() > 1 && outString[0] == '\"') {
                    outString = outString.substr(1);

                }
                if (outString.length() > 1 && outString[outString.length() - 1] == '\"') {
                    outString = outString.substr(0, outString.length() - 1);

                }

            } else {
                //std::cout << "READ 22 " << line[line.length()-1] << "     " << line << std::endl;
                outString = line.substr(1);
                boost::algorithm::trim(outString);
                std::string result;
                std::getline(inputStream, result, ']');
                boost::algorithm::trim(result);
                if (result[0] == '\"') {
                    result = result.substr(1);
                }
                if (result[result.length() - 1] == '\"') {
                    result = result.substr(0, result.length() - 1);
                }
                outString += ' ' + result;
            }
        } else {

            if (line.find('\"') != std::string::npos) {
                outString = line.substr(0, line.length() - 1);
            } else {
                outString = line;
            }
        }

        boost::algorithm::trim(outString);
        return true;
    }

    return false;
}
