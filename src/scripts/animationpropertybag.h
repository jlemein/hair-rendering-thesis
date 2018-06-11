#ifndef ANIMATIONPROPERTYBAG_H
#define ANIMATIONPROPERTYBAG_H

#include "propertybag.h"
#include <boost/filesystem.hpp>

class AnimationPropertyBag : public PropertyBag
{
private:
    
public:
    AnimationPropertyBag();
    void addLine();
    
    void replacePropertiesInFileForFrame(const boost::filesystem::path& inputFileName, const boost::filesystem::path& outputFileName, unsigned int frame);
};

#endif // ANIMATIONPROPERTYBAG_H
