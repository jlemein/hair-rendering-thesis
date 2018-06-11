#ifndef CONSTANTPROPERTY_H
#define CONSTANTPROPERTY_H

#include "property.h"
#include <string>

class ConstantProperty : public Property
{
private:
    std::string mValue;

public:
    ConstantProperty(const std::string& name, const std::string& value);

    virtual std::string operator()(float t) const;
    virtual std::string operator[](int keyframe) const;
    virtual std::string getValue() const;
};

#endif // CONSTANTPROPERTY_H
