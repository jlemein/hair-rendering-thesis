#include "hairstruct.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "pbrthairparser.h"

std::ostream& operator<<(std::ostream& out, const Point& p) {
  out << "[" << p.x << " " << p.y << " " << p.z << "]";
}

/**
 * @todo Clean this messy code
 */
std::istream& operator>>(std::istream& istream, Hair& hairData) {

    PbrtHairParser::parseFile(istream, hairData);
    return istream;
}
