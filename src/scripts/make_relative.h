#ifndef MAKE_RELATIVE_H
#define MAKE_RELATIVE_H

#include <string>
#include <boost/filesystem.hpp>

// TODO: make proper version
boost::filesystem::path make_relative( boost::filesystem::path from, boost::filesystem::path to) {
    std::cout << "WARNING: using make_relative is only working correctly when from and to are both relative starting from same directory" << std::endl;
    std::cout << "TEST: " << boost::filesystem::path("../..");
    return boost::filesystem::path(std::string("../../") + to.string());
}


#endif // MAKE_RELATIVE_H
