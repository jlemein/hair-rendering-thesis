#ifndef PROPERTY_H
#define PROPERTY_H

#include <string>

/**
 * @brief The Property class
 */
class Property
{
public:
    Property(const std::string& propertyName);

    /**
     * @brief operator () Returns interpolated value at time
     * @param time
     * @return
     */
    virtual std::string operator()(float time) const = 0;

    /**
     * @brief operator [] The keyframe operator returns the property value for the specified keyframe
     * @param keyFrame The keyframe
     * @return The property value at a specific keyframe
     */
    virtual std::string operator[](int keyFrame) const = 0;

    /**
     * @brief getValue Returns default value, when no keyframe or time have been specified
     */
    virtual std::string getValue() const = 0;

    /**
     * @brief getKeyFrameCount returns the number of keyframes defined for this property
     * @return default is 1, but can be any positive number
     */
    virtual int getKeyFrameCount() const;

    std::string name;

};

#endif // PROPERTY_H
