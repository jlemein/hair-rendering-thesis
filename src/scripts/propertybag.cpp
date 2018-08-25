#include "propertybag.h"

#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
//#include <boost/algorithm/string/predicate.hpp>

#include "animatedproperty.h"
#include "constantproperty.h"
#include "property.h"

#include <fstream>
#include <iostream>
#include <sstream>

//#include <boost/smart_ptr/shared_ptr.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/algorithm/string.hpp>

//#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

/**
 * An animated property is a property that describes a discrete function. Two types of animations can be described:
 * - A keyframe animation, meaning the numbers are absolute indices in the animation. #1 means first frame of animation, #2 means second frame.
 *
 *
 *
 * Examples of animation descriptions:
 *
 * For keyframe 0 to keyframe
 * @viewAngle #1: { 10.0 } - #4 { 40.0 }  -- implies that keyframe 1 has orientation 10.0, keyframe 2 and 3 also have orientation 10.0, and from keyframe 4 the orientation is set to 40.0.
 * @viewAngle #1: { 10.0 } -> #4 { 40.0 }  -- implies that keyframe 1 has orientation 10.0, keyframe 2 and 3 are linearly interpolated, meaning keyframe 2 has orientation 20.0, keyframe 3 has orientation 30.0 and keyframe 4 has orientation of 40.0.
 * @viewAngle 0: { 45.0 } -
 * @viewAngle 0: { 45.0 } -> 1: { 90.0 }
 *
 * @viewAngle #1: 10.0 - #4: 20.0
 * @viewAngle #1: 10.0 ~ #4: 20.0
 */

void PropertyBag::parseAnimatedProperty(std::stringstream& ss) {
    std::string line;
    std::string element;
    std::string buffer;
    ss >> element;

    if (element.empty() || !boost::starts_with(element, "@")) {
        std::cout << "Requested to parse an animated property, but it does not start with @" << std::endl;
        return;
    }

    std::string key = element.substr(1);
    AnimatedProperty* property = new AnimatedProperty(key);

    bool stopCondition = false;

    do {
        // read whether keyframe is provided, and if yes, then which frame number
        bool isAbsoluteKeyFrame = true;
        int keyframe;
        std::getline(ss, element, ':');
        boost::trim(element);

        if (element.length() > 1 && element[0] == '#') {
            isAbsoluteKeyFrame = true;
            keyframe = boost::lexical_cast<int>(element.substr(1));
        } else {
            std::cout << "For now only keyframes can be provided using hash (#)" << std::endl;
            keyframe = boost::lexical_cast<int>(element);
        }

        ss >> line;
        boost::trim(line);

        // read value
        float value = boost::lexical_cast<float>(line);

        // read animationType
        ss >> element;
        boost::trim(element);

        AnimationType animationType = AnimationType::Fixed;
        if (!element.empty()) {
            char interpolationSymbol = boost::lexical_cast<char>(element);

            if (interpolationSymbol == '-') {
                animationType = AnimationType::Fixed;
            } else if (interpolationSymbol == '~') {
                animationType = AnimationType::Linear;
            } else if (interpolationSymbol == ';') {
                stopCondition = true;
            } else {
                std::cout << "Unsupported interpolation symbol provided. Using fixed interpolation by default" << std::endl;
            }
        }

        property->addKeyFrame(keyframe, isAbsoluteKeyFrame, animationType, boost::lexical_cast<float>(value));
    } while (!ss.str().empty() && !stopCondition);

    mProperties[key] = property;

}

void PropertyBag::parseConstantProperty(std::stringstream& ss) {
    std::string key, value;
    ss >> key;
    getline(ss, value);
    boost::algorithm::trim(value);

    //std::cout << key << ": '" << value << "'" << std::endl << std::endl;

    mProperties[key] = new ConstantProperty(key, value);
}

void PropertyBag::addPropertiesFromFile(const std::string fileName) {
  std::ifstream infile(fileName.c_str());
  if (infile.fail()) {
    throw "Failed to open input file '" + fileName + "'";
  }

  std::string line;
  while(!infile.eof() && getline(infile, line)) {
      if (!line.empty() && !boost::starts_with(line, "#")) {
          std::string key, value;
          //std::cout << "Read line: '" << line << "'" << std::endl;
          boost::algorithm::trim(line);
          std::stringstream ss(line);

          //ss >> key;

          if (ss.peek() == '@') {
              parseAnimatedProperty(ss);
          } else {
              parseConstantProperty(ss);
          }
      }
  }

  infile.close();
  std::cout << std::endl;
}

void PropertyBag::addProperty(const std::string& key, const std::string& value) {
    mProperties[key] = new ConstantProperty(key, value);
}

void PropertyBag::addPropertiesFromLine(const std::string& line) {
  std::vector<std::string> keyValues;
  boost::split(keyValues, line, boost::is_any_of("\t, "), boost::token_compress_on);

  for (auto keyValue : keyValues) {
    std::vector<std::string> splitted;
    boost::split(splitted, keyValue, boost::is_any_of(":,="));
    if (splitted.size() != 2) {
      throw "Cannot parse key-value pair '" + keyValue + "'";
    } else {
      mProperties[splitted[0]] = new ConstantProperty(splitted[0], splitted[1]);
    }
  }
}

std::string PropertyBag::getProperty(const std::string& key) {
    if (mProperties.find(key) != mProperties.end()) {
        return mProperties[key]->getValue();
    } else {
        std::string msg = "'" + key + "' does not exist in property bag";
        std::cout << "ERROR: " << msg << std::endl;
        throw msg;
    }
//    return "";
}

int PropertyBag::getKeyFrameCount() const {
    int numberKeyFrames = 1;

    for (auto property : mProperties) {
        numberKeyFrames = std::max(numberKeyFrames, property.second->getKeyFrameCount());
    }

    return numberKeyFrames;

}

std::string formatKeyFrame(int keyIndex) {
    std::string keyframeStr;
    if (keyIndex < 10) {
        keyframeStr = "000";
    } else if (keyIndex < 100) {
        keyframeStr = "00";
    } else if (keyIndex < 1000) {
        keyframeStr = "0";
    }
    keyframeStr += boost::lexical_cast<std::string>(keyIndex);
    return keyframeStr;
}

void PropertyBag::replacePropertiesInFile(const std::string& inputFileName, const std::string& outputFileName) {
  namespace fs = boost::filesystem;

  // check how many keyframes there are defined
  int keyFrameCount = getKeyFrameCount();

  for (int keyIndex = 0; keyIndex < keyFrameCount; ++keyIndex) {
      //
      // Opening input file stream
      std::ifstream infile(inputFileName.c_str());
      if (infile.fail()) {
        throw "Failed to open input file '" + inputFileName + "'";
      }

      //
      // Writing output file
      std::string formattedKeyFrameStr = formatKeyFrame(keyIndex);

      fs::path p = fs::path(outputFileName);
      std::string newOutputFileName = (p.parent_path() / p.stem()).string()
              + "_" + formattedKeyFrameStr
              + p.extension().string();

      std::ofstream outfile(newOutputFileName.c_str());
      if (outfile.fail()) {
        infile.close();
        throw "Cannot write to output file '" + outputFileName + "'";
      }

      this->addProperty("KEY_FRAME", formattedKeyFrameStr);

      std::string line;
      while( std::getline(infile, line) ) {
        for (auto prop : mProperties) {
          boost::replace_all(line, "$"+prop.first, (*prop.second)[keyIndex]);
        }

        outfile << line << std::endl;
      }

      outfile.close();
      infile.close();

      std::cout << newOutputFileName << "[Done]" << std::endl;
  }
}

void PropertyBag::writePropertiesToFile(const std::string& fileName) {
//  std::ofstream outfile(fileName.c_str());
//  if (outfile.fail()) {
//    throw "Cannot open file '" + fileName + "' to write properties";
//  }

//  for (auto keyValuePair : mProperties) {
//    outfile << keyValuePair.first << " " << keyValuePair.second << std::endl;
//  }

//  outfile.close();
}
