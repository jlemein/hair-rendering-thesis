#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "hairstruct.h"

int main (int argc, char** argv) {
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;

  fs::path inputFileName;
  fs::path outputFileName;
  float width0, width1;

  // describe commands the user can enter
  po::options_description desc("Allowed options");
  desc.add_options()
           ("help", "Processes hair ")
           ("input", po::value<std::string>(), "Input hair file (*.pbrt) to use")
           ("output", po::value<std::string>(), "output file to write result to")
           ("width0", po::value<float>(), "set width of hair")
           ("width1", po::value<float>(), "width of second axis for the hair")
           ("random", po::value<std::string>(), "indicates to use random widths for the hair");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("input")) {
    inputFileName = fs::path(vm["input"].as<std::string>());
  } else {
    std::cout << "Required input file name --input <filename> is not specified" << std::endl;
    return 1;
  }

  if (vm.count("output")) {
    outputFileName = fs::path(vm["output"].as<std::string>());
  } else {
      outputFileName = fs::path(inputFileName).replace_extension("pbrt");
      std::cout << "Output file is not specified, default write to: '" << outputFileName << "'" << std::endl;
      if (outputFileName == inputFileName) {
          std::cout << "Input file is equal to output path, cannot continue" << std::endl;
          return 1;
      }
  }

  width0 = 0.001;
  width1 = 0.001;

  if (vm.count("width0")) {
      width0 = vm["width0"].as<float>();
  }

  if (vm.count("width1")) {
      width1 = vm["width1"].as<float>();
  }

  std::cout << "Writing hair file to '" << outputFileName << "' using specified widths: " << width0 << ", " << width1 << std::endl;

  Hair hair;
  std::ifstream input(inputFileName.string().c_str());
  input >> hair;

  for(auto& curve : hair.curves) {
      curve.width0 = width0;
      curve.width1 = width1;
  }

  std::cout << "Writing hair file to output" << std::endl;
  //std::ofstream output(outputFileName.string().c_str());
  //output << hair;
  hair.writeFile(outputFileName.string());
  std::cout << "Done" << std::endl;

}
