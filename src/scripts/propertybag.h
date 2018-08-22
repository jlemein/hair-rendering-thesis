#ifndef PROPERTY_BAG_H
#define PROPERTY_BAG_H

#include <boost/unordered_map.hpp>
#include "property.h"
#include <sstream>

using property_map_t = boost::unordered::unordered_map<std::string, Property*>;

/**
 *  author: Jeffrey Lemein
 */
class PropertyBag {
private:
    property_map_t mProperties;

    void parseAnimatedProperty(std::stringstream& ss);
    void parseConstantProperty(std::stringstream& ss);

public:
    /**
    * Reads a property file and returns the map containing the property-value pairs.
    * A property can be defined in 2 ways: as a static property and an animated property
    *
    * Animated property
    * - Starts with '@', meaning at keyframe #<number>, starting with 0 as first keyframe
    *   Examples: @0: viewRotation=90; @4: viewRotation=360
    */
    void addPropertiesFromFile(const std::string fileName);
    void addProperty(const std::string& key, const std::string& value);
    void addPropertiesFromLine(const std::string& line);

    /**
     * Returns property value from property bag, throws error if not exist
     */
    std::string getProperty(const std::string& key);

    /**
      * Returns the number of keyframes defined for this property map
      */
    int getKeyFrameCount() const;
    
    /**
    * Replaces all properties in inputFileName with the values provided in properties-map.
    * Result is written out to outputFileName
    */
    void replacePropertiesInFile(const std::string& fileName, const std::string& outputFileName);

    /**
    * Writes properties to file with name fileName.
    */
    void writePropertiesToFile(const std::string& fileName);
};

#endif // PROPERTY_BAG_H
