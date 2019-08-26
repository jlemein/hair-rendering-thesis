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
    std::vector<std::pair<int, std::string> > mKeyFrames;


public:
    AnimatedProperty(const std::string& key);

    void addKeyFrame(int keyFrame, bool isAbsoluteKeyFrame, AnimationType animationType, std::string value);

    /**
     * @brief findPosition Finds index in key frame list, where the specified key frame index fits between.
     * For example: if key frames are specified [0, 5, 10], then key frame index 7 fits between 5 and 10, corresponding to index 2.
     * @param keyFrame The index of the keyframe
     * @return Index to the key frame followed after the specified index. 0 in case keyframe comes before first specified keyframe, and
     * index = length of keyframe array if key frame comes after the last key frame.
     */
    int findPosition(int keyFrame) const;

    virtual std::string operator()(float t) const;
    virtual std::string operator[](int keyframe) const;
    virtual std::string getValue() const;

    virtual int getKeyFrameCount() const;
};

#endif // ANIMATEDPROPERTY_H
