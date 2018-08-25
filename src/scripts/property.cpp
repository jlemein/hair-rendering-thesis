#include "property.h"
#include <iostream>

Property::Property(const std::string& name)
    : name(name)
{

}

int Property::getKeyFrameCount() const {
    return 1;
}
