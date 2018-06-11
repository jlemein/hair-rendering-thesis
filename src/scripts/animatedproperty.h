#ifndef ANIMATEDPROPERTY_H
#define ANIMATEDPROPERTY_H

#include "property.h"
#include "animationtype.h"
#include <string>
#include <vector>

template <typename T>
struct KeyFrame {
    int id;
    T value;
    AnimationType animationType;
};

class AnimatedProperty : public Property
{
private:
    std::vector<std::pair<int, float> > mKeyFrames;

public:
    AnimatedProperty(const std::string& key);

    void addKeyFrame(int keyFrame, bool isAbsoluteKeyFrame, AnimationType animationType, float value);

    virtual std::string operator()(float t) const;
    virtual std::string operator[](int keyframe) const;
    virtual std::string getValue() const;
};

#endif // ANIMATEDPROPERTY_H
