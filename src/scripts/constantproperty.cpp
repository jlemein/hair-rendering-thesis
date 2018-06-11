#include "constantproperty.h"

ConstantProperty::ConstantProperty(const std::string& name, const std::string& value)
    : Property(name)
{
    this->mValue = value;
}

std::string ConstantProperty::operator()(float t) const {
    return mValue;
}

std::string ConstantProperty::operator[](int keyframe) const {
    return mValue;
}

std::string ConstantProperty::getValue() const {
    return mValue;
}
