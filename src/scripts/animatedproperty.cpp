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

float linearExtrapolate(float valueA, float valueB, float ratio) {
    return (1.0f - ratio) * valueA + ratio * valueB;
}

std::string AnimatedProperty::operator [](int key) const {
    float returnValue;
    for (int i=0; i<mKeyFrames.size(); ++i) {
        if (mKeyFrames[i].first >= key) {
            if (i==0) {
                returnValue = mKeyFrames[0].second;
            } else if (i+1 == mKeyFrames.size()) {
                returnValue = mKeyFrames[i].second;
            } else {
                float valuePrev = mKeyFrames[i-1].second;
                float valueNext = mKeyFrames[i].second;
                float duration = (mKeyFrames[i].first - mKeyFrames[i-1].first);
                float offset = key - mKeyFrames[i-1].first;
                returnValue = linearExtrapolate(valuePrev, valueNext, offset / duration);
            }
            return boost::lexical_cast<std::string>(returnValue);
        }
    }
}

std::string AnimatedProperty::operator ()(float t = 0.0f) const {
    std::cout << "Not supported yet" << std::endl;
    throw "Not supported yet";
}

std::string AnimatedProperty::getValue() const {
    return this->operator [](0);
}
