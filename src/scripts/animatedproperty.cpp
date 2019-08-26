#include "animatedproperty.h"
#include <boost/lexical_cast.hpp>
#include <iostream>

AnimatedProperty::AnimatedProperty(const std::string& key)
    : Property(key)
{
    std::cout << "HELLO WORLD" << std::endl;
}

void AnimatedProperty::addKeyFrame(int keyFrame, bool isAbsoluteKeyFrame, AnimationType animationType, std::string value) {
    this->mKeyFrames.push_back(std::pair<int, std::string>(keyFrame, value));
}

float linearInterpolate(float valueA, float valueB, float ratio) {
    return (1.0f - ratio) * valueA + ratio * valueB;
}

int AnimatedProperty::findPosition(int keyFrame) const {
    // assume key frames are sorted
    int foundIndex = 0;
    if (mKeyFrames.empty()) {
        std::cout << "No key frames specified, unable to find position" << std::endl;
        throw "No key frames specified, unable to find position";
    }
    for( int i=0; i<mKeyFrames.size(); ++i ) {
        if (keyFrame < mKeyFrames[i].first) {
            return i;
        }
    }
    return mKeyFrames.size();
}

std::string AnimatedProperty::operator [](int key) const {
    float returnValue;

    try{
        int position = findPosition(key);
        if (position == 0) {
            return mKeyFrames[position].second;
        } else if (position == mKeyFrames.size()) {
            return mKeyFrames[position-1].second;
        } else {
            std::string valuePrev = mKeyFrames[position-1].second;
            std::string valueNext = mKeyFrames[position].second;

            try{
                // interpolate
                float duration = (mKeyFrames[position].first - mKeyFrames[position-1].first);
                float offset = key - mKeyFrames[position-1].first;
                returnValue = linearInterpolate(boost::lexical_cast<float>(valuePrev), boost::lexical_cast<float>(valueNext), offset / duration);
            } catch (std::exception& e) {
                return valueNext;
            }
        }

        return boost::lexical_cast<std::string>(returnValue);
    } catch (const std::exception& e) {
        std::cout << "Animated property could not be read for key = '" << this->name << "', mKeyFrames.size = " << mKeyFrames.size() << std::endl;
    }
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

    for (auto keyFrame : mKeyFrames) {
        maxKeyFrameIndex = std::max(maxKeyFrameIndex, keyFrame.first);
    }

    return maxKeyFrameIndex+1;
}
