#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

struct Point {
  float x, y, z;
};

struct Curve {
    std::string type;
    std::vector<Point> points;
    float width0;
    float width1;
};

struct Hair {

    /**
     * @brief nextWord Returns the word between quotes ""
     * @return
     */

    std::vector<Curve> curves;

};

bool nextWord(std::istream& istream, std::string& outString) {
    // find beginning    
    while (std::getline(istream, outString, '\"')) {
        boost::algorithm::trim(outString);
        if (outString.length() > 0) {
            return true;
        }
    }

    return false;
}

bool nextValues(std::istream& inputStream, std::string& outString) {
    std::string line;
    inputStream >> line;

    boost::algorithm::trim(line);
    if (line.empty()) {
        return false;
    } else {
        if (line[0] == '\"' && line.length() > 1) {
            if (line[line.length()-1] == '\"') {
                outString = line.substr(1, line.length()-1);
            } else {
                outString = line.substr(1);
                boost::algorithm::trim(outString);
                std::string result;
                std::getline(inputStream, result, '\"');
                outString += ' ' + result;
            }
        }
        else if (line[0] == '[') {
            //std::cout << "READ 2 " << line << std::endl;

            if (line.length() > 1 && line[line.length()-1] == ']') {
                outString = line.substr(1, line.length()-1);

                boost::algorithm::trim(outString);
                if (outString.length() > 1 && outString[0] == '\"') {
                    outString = outString.substr(1);

                }
                if (outString.length() > 1 && outString[outString.length()-1] == '\"') {
                    outString = outString.substr(0, outString.length()-1);

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
                if (result[result.length()-1] == '\"') {
                    result = result.substr(0, result.length()-1);
                }
                outString += ' ' + result;
            }
        }
        else {

            if (line.find('\"') != std::string::npos) {
                outString = line.substr(0, line.length()-1);
            } else {
                outString = line;
            }
        }

        boost::algorithm::trim(outString);
        return true;
    }

    return false;
}

std::istream& operator>>(std::istream& istream, Hair& data) {

    std::cout << "Parsing hair file" << std::endl;
    data.curves.reserve(10000);

    std::string word;
    while(nextWord(istream, word)) {
        std::string values;

        Curve c;
        c.points.reserve(12);
        if (word == "Shape") {
            nextValues(istream, c.type);
        } else if (word == "float width") {
            nextValues(istream, values);
            c.width0 = boost::lexical_cast<float>(values);
        } else if (word == "point P") {
            std::string points;
            nextValues(istream, points);
            std::stringstream ss(points);
            Point p;
            while (ss >> p.x >> p.y >> p.z) {
                c.points.push_back(p);
            }
            data.curves.push_back(c);
        }
        else {
            nextValues(istream, values);
        }
    }

    std::cout << "Finished reading. Number of curves: " << data.curves.size() << std::endl;

//        if (key == "Shape") {
//            nextWord(istream);
//            getline(istream, value, '"');
//        }

    return istream;
}


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
    infile >> hair;

    return 0;
}
