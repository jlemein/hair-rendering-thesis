#include "property.h"
#include <iostream>

Property::Property(const std::string& name)
    : name(name)
{

}

int Property::getKeyFrameCount() const {
    std::cout << name << " has " << 1 << " key frame" << std::endl;
    return 1;
}
