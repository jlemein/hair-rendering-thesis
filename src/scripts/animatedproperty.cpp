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

    int position = findPosition(key);
    if (position == 0) {
        returnValue = mKeyFrames[position].second;
    } else if (position == mKeyFrames.size()) {
        returnValue = mKeyFrames[position-1].second;
    } else {
        // interpolate
        float valuePrev = mKeyFrames[position-1].second;
        float valueNext = mKeyFrames[position].second;
        float duration = (mKeyFrames[position].first - mKeyFrames[position-1].first);
        float offset = key - mKeyFrames[position-1].first;
        returnValue = linearInterpolate(valuePrev, valueNext, offset / duration);
    }

//    if (mKeyFrames.size() > 2) {
//        std::cout << "Animated property: " << name << std::endl;
//        for (auto keyframe : mKeyFrames) {
//            std::cout << "#" << keyframe.first << ": " << keyframe.second << std::endl;
//        }
//        std::cout << "Requested frame #" << key << std::endl;

//        if (position <= 0) {
//            std::cout << "Position is before first keyframe" << std::endl;
//        } else if (position >= mKeyFrames.size()) {
//            std::cout << "Position is after last specified keyframe" << std::endl;
//        } else {
//            std::cout << "Position is: " << position << ", meaning between keyframes: " << mKeyFrames[position-1].first << " and " << mKeyFrames[position].first << std::endl;
//        }
//    }

//    std::cout << this->name << " at #" << key << ": " << returnValue << std::endl;

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
