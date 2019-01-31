#include "hairstruct.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "pbrthairparser.h"

//std::ostream& operator<<(std::ostream& out, const Point& p) {
//    out << "[" << p.x << " " << p.y << " " << p.z << "]";
//}

Hair::Hair() {
}

void Hair::writeFile(const std::string& outputFileName) {
    std::ofstream out(outputFileName.c_str());
    if (out.fail()) {
        std::cout << "Failed to open " << outputFileName << std::endl;
    }

    for (const auto& fiber : this->fibers) {
        out << "Shape \"" << "curve" << "\" "
                << "\"string type\" [ \"" << "cylinder" << "\" ] "
                << "\"point P\" [ ";
        for (const auto& pt : fiber.curve.getControlPoints()) {
            out << pt.x << " " << pt.y << " " << pt.z << " ";
        }
        out << "] "
                << "\"float width0\" [ " << fiber.width0 << " ] "
                << "\"float width1\" [ " << fiber.width1 << " ] ";
        out << std::endl;
    }
}

/**
 * @todo Clean this messy code
 */
std::istream& operator>>(std::istream& istream, Hair& hairData) {

    PbrtHairParser::parseFile(istream, hairData);
    return istream;
}

std::ostream& operator<<(std::ostream& out, const Hair& hair) {
    for (const auto& fiber : hair.fibers) {
        out << "Shape \"" << "curve" << "\" "
                << "\"string type\" [ \"" << hair.getCurveType() << "\" ] "
                << "\"point P\" [ ";
        for (const auto& pt : fiber.curve.getControlPoints()) {
            out << pt.x << " " << pt.y << " " << pt.z << " ";
        }
        out << "] "
                << "\"float width0\" [ " << fiber.width0 << " ] "
                << "\"float width1\" [ " << fiber.width1 << " ] ";
        out << std::endl;
    }
    return out;
}

std::string Hair::getCurveType() const {
    if (this->curveType == CurveType::CYLINDER) {
        return "cylinder";
    }
    if (this->curveType == CurveType::RIBBON) {
        return "ribbon";
    }
    if (this->curveType == CurveType::FLAT) {
        return "flat";
    }

    //default
    return "flat";
}


/*
 * Shape "curve" "string type" [ "cylinder" ] "point P" [ -11.439552 49.158043 -14.240796 -11.473633 49.230301 -14.286941 -11.509600 49.301491 -14.334203 -11.547457 49.371628 -14.382583  ] "float width0" [ 0.100000 ] "float width1" [ 0.100000 ]
 */
