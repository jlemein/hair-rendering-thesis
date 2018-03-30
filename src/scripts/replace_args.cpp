#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>

#include "propertybag.h"

/********************************************************************************************************************************
 *
 * Replaces all property names in input file with the specified key-values and writes it out to an output file.
 * Property names in the input file are marked by prefixing it with a dollar sign '$', such as $firstName or $sampleCount.
 *
 * author: jeffrey lemein
 * date: 22 march 2018
 *
 ********************************************************************************************************************************/

int main(int argc, char** argv) {
  std::string inputFileName;
  std::string outputFileName;
  PropertyBag propertyBag;

  namespace po = boost::program_options;
  // describe commands the user can enter
  po::options_description desc("Allowed options");
  desc.add_options()
	       ("help", "Replaces variables in input file with properties specified in a properties file, or by properties provided on the command line. Result is written to output file ")
	       ("input", po::value<std::string>(), "unprocessed PBRT scene containing arguments to be replaced")
	       ("output", po::value<std::string>(), "output file to write result to")
	       ("propertyfile", po::value<std::string>(), "add property file containing key values")
	       ("properties", po::value<std::string>(), "specify a list of key value pairs as 'key1:value1 key2:value2 ... keyn:valuen'");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("input")) {
    inputFileName = vm["input"].as<std::string>();
  } else {
    std::cout << "Required input file name --input <filename> is not specified" << std::endl;
    return 1;
  }
  
  if (vm.count("output")) {
    outputFileName = vm["output"].as<std::string>();
  } else {
    std::cout << "Required output file name, using --output <filename> is not specified" << std::endl;
    return 1;
  }

  if (vm.count("propertyfile")) {
    std::string fileName = vm["propertyfile"].as<std::string>();
    std::cout << "Property file specified: " << fileName << "\n";
    try {
      propertyBag.addPropertiesFromFile(fileName);
    } catch(const char* e) {
      std::cout << "Error parsing property file: " << e << std::endl;
    }
  }

  if (vm.count("properties")) {
    std::string line = vm["properties"].as<std::string>();
    propertyBag.addPropertiesFromLine(line);
  }

  try {
    propertyBag.writePropertiesToFile(inputFileName + ".properties");
    propertyBag.replacePropertiesInFile(inputFileName, outputFileName);
  } catch (const char* e) {
    std::cout << "Exception when replacing arguments: " << e << std::endl;
  }
  
  return 0;
}
