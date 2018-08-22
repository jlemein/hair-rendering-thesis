#include "animatedproperty.h"
#include <boost/lexical_cast.hpp>
#include <iostream>

AnimatedProperty::AnimatedProperty(const std::string& key)
    : Property(key)
{

}

void AnimatedProperty::addKeyFrame(int keyFrame, bool isAbsoluteKeyFrame, AnimationType animationType, float value) {
    this->mKeyFrames.push_back(std::pair<int, float>(keyFrame, value));
}

float linearInterpolate(float valueA, float valueB, float ratio) {
    return (1.0f - ratio) * valueA + ratio * valueB;
}

std::string AnimatedProperty::operator [](int key) const {
    float returnValue;

    for (int i=0; i<mKeyFrames.size(); ++i) {
        // The case when key frame comes before the sequence of defined key frames (then just return first value)
        if (i == 0 && key <= mKeyFrames[i].first) {
            returnValue = mKeyFrames[0].second;
        } else if (key <= mKeyFrames[i].first) {
            // interpolate
            float valuePrev = mKeyFrames[i-1].second;
            float valueNext = mKeyFrames[i].second;
            float duration = (mKeyFrames[i].first - mKeyFrames[i-1].first);
            float offset = key - mKeyFrames[i-1].first;
            returnValue = linearInterpolate(valuePrev, valueNext, offset / duration);
        } else {
            // otherwise key frame lays after the sequence defined, so return last value of list
            returnValue = mKeyFrames[i].second;
        }
    }

    return boost::lexical_cast<std::string>(returnValue);
}

std::string AnimatedProperty::operator ()(float t = 0.0f) const {
    std::cout << "Not supported yet" << std::endl;
    throw "Not supported yet";
}

std::string AnimatedProperty::getValue() const {
    return this->operator [](0);
}

int AnimatedProperty::getKeyFrameCount() const {
    int maxKeyFrameIndex = 0;
    std::cout << "Requesting keyframe count for @" << this->name << std::endl;
    std::cout << "Key frames: " << mKeyFrames.size() << std::endl;

    for (auto keyFrame : mKeyFrames) {
        std::cout << "Index: " << keyFrame.first << std::endl;
        maxKeyFrameIndex = std::max(maxKeyFrameIndex, keyFrame.first);
    }

    return maxKeyFrameIndex+1;
}
