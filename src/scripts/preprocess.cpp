#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>

#include "propertybag.h"
#include "hairsimplify.h"
#include "make_relative.h"

using namespace std;
namespace fs = boost::filesystem;


void generatePbrtScene();

/* Output directory contains collection of directories representing render-setups and render-results */
const fs::path OUTPUT_DIRECTORY_NAME("output");

/* Contains generated data, such as simplified models, lookup-data */
const fs::path GENERATED_DIRECTORY_NAME("generated");

/* Contains scene, model and user property files */
const fs::path PROPERTY_DIRECTORY_NAME("properties");

const string PROPERTY_SUFFIX = ".properties";
const string SCENE_PROPERTY_FILE_NAME = "scene.properties";
const string MODEL_PROPERTY_FILE_NAME = "model.properties";
const string USER_PROPERTY_FILE_NAME = "user.properties";


/** User specified properties */
PropertyBag modelProperties, sceneProperties, userProperties;
int simplifyPercentage;
bool userRequestSimplifiedModel = false;

fs::path outputDirectory("./output"), sessionResultsDirectory, generatedContentDirectory, propertyFilesDirectory;
fs::path modelFilePath, sceneFilePath;
fs::path scenePropertiesFilePath, modelPropertiesFilePath, userPropertiesFilePath;


void createDirectoriesIfNotExists(fs::path& directory) {
  if (fs::create_directories(directory)) {
    cout << "Succesfully created directory " << outputDirectory << endl;
  } else {
    cout << "Directory '" << outputDirectory << " already exists\n" << endl;
  }  
}

/** Returns the current time in a readable format */
string getCurrentTime() {
  namespace pt = boost::posix_time;
  return pt::to_simple_string(pt::second_clock::local_time());
}

void parseSceneFile() {
  try {
    cout << "Parsing scene properties in " << sceneFilePath.string() << ".properties" << " ... ";
    sceneProperties.addPropertiesFromFile(sceneFilePath.string() + ".properties");
    cout << "[done]" << endl;
  } catch (const char* e) {
    cout << "[failed]" << endl;
  }
}

/**
 * Parse properties from a property file in the specified properties collection
 */
void readProperties(const fs::path& propertyFilePath, PropertyBag& properties) {
  try {
      cout << "Reading property file " << propertyFilePath.string() << "... ";
      properties.addPropertiesFromFile(propertyFilePath.string());
      cout << "[done]" << endl;
    } catch (string e) {
        cout << e << endl;
    }
}

void simplifyAndWriteModel(const fs::path& simplifiedModelFilePath) {
  cout << "Simplifying model..." << endl;
  if (fs::exists(simplifiedModelFilePath)) {
    cout << "Already exists... [done]" << endl;
  } else {
    cout << "Simplify model to " << simplifyPercentage << " percent of original... ";
    HairSimplify hair(modelFilePath.string());
    hair.reduceToPercentage(simplifiedModelFilePath.string(), static_cast<double>(simplifyPercentage));
    cout << "[done]" << endl;
  }
  
  cout << "Generate voxel grid... [todo]" << endl;
}

/**
 * @brief generatePbrtScene Generates a PBRT scene from the template scene, by substituting all variables by the properties specified in the property files.
 * @param basePath The path to the directory where the PBRT file is to be generated
 */
void generatePbrtScene(fs::path basePath, bool deriveSceneFileFromProperties = true) {
  cout << "Generating PBRT scene..." << endl;

  // Read and merge properties in property files
  PropertyBag properties;
  readProperties((propertyFilesDirectory / SCENE_PROPERTY_FILE_NAME).string(), properties);
  readProperties((propertyFilesDirectory / MODEL_PROPERTY_FILE_NAME).string(), properties);
  readProperties((propertyFilesDirectory / USER_PROPERTY_FILE_NAME).string(), properties);

  // Output resulting PBRT file
  if (deriveSceneFileFromProperties) {
      sceneFilePath = fs::path(fs::path(properties.getProperty("input_hairmodel_filename")).stem().string() + ".pbrt");
  }

  cout << "Reading scene file: " << sceneFilePath << endl;
  string pbrtFileName = sceneFilePath.filename().string();
  fs::path pbrtFilePath(basePath / pbrtFileName);
  cout << "PBRT file path: " << pbrtFilePath << endl;
  
  // Replace properties in scene file
  cout << "replace properties in: " << sceneFilePath << " and writing to: " << pbrtFilePath.string() << endl;
  properties.replacePropertiesInFile(sceneFilePath.string(), pbrtFilePath.string());
  cout << "[done]" << endl;
}

/**
 * @brief startProcessing Sets up a test case, by copying property files from original and
 * possibly simplifying hair models
 */
void startProcessing() {
  cout << "Start processing..." << endl;

  // Create directories
  sessionResultsDirectory = outputDirectory / getCurrentTime();
  propertyFilesDirectory = sessionResultsDirectory / PROPERTY_DIRECTORY_NAME;
  generatedContentDirectory = outputDirectory / GENERATED_DIRECTORY_NAME;

  createDirectoriesIfNotExists(sessionResultsDirectory);
  createDirectoriesIfNotExists(propertyFilesDirectory);
  createDirectoriesIfNotExists(generatedContentDirectory);
  
  // copy original scene file
  fs::copy_file(sceneFilePath, sessionResultsDirectory/sceneFilePath.filename()); // Copy original scene file

   // Read properties
  readProperties(scenePropertiesFilePath, sceneProperties);
  readProperties(modelPropertiesFilePath, modelProperties);
  readProperties(userPropertiesFilePath, userProperties);
  
  // simplify model (or use original) and write to properties
  fs::path originalSceneDirectory = make_relative(sessionResultsDirectory, sceneFilePath.parent_path());
  sceneProperties.addProperty("currentDirectory", originalSceneDirectory.string());
  if (userRequestSimplifiedModel) {
    string simplifiedModelFileName = modelFilePath.filename().string() + boost::lexical_cast<string>(simplifyPercentage);
    fs::path simplifiedModelFilePath = generatedContentDirectory / simplifiedModelFileName;
    simplifyAndWriteModel(simplifiedModelFilePath.string());
    sceneProperties.addProperty("input_hairmodel_filename", simplifiedModelFilePath.string());
    cout << "Using generated model at location: " << simplifiedModelFilePath << endl;
  } else {
      fs::path relativeModelFilePath = make_relative(sessionResultsDirectory, modelFilePath);
      cout << "No simplification needed. Using original model: " << relativeModelFilePath << endl;
      sceneProperties.addProperty("input_hairmodel_filename", relativeModelFilePath.string());
  }

  // Write property files
  sceneProperties.writePropertiesToFile((propertyFilesDirectory / SCENE_PROPERTY_FILE_NAME).string());
  modelProperties.writePropertiesToFile((propertyFilesDirectory / MODEL_PROPERTY_FILE_NAME).string());
  userProperties.writePropertiesToFile((propertyFilesDirectory / USER_PROPERTY_FILE_NAME).string());

  // generate a pbrt scene
  createDirectoriesIfNotExists(sessionResultsDirectory);
  generatePbrtScene(sessionResultsDirectory, false);
}

int main(int argc, char** argv) {
    namespace po = boost::program_options;

    if (argc == 2) {
        fs::path directory(argv[1]);
        cout << "Preprocessing in directory: " << directory << endl;

        propertyFilesDirectory = directory / PROPERTY_DIRECTORY_NAME;
        generatePbrtScene(directory, true);
        return 0;
    }
    
    // describe commands the user can enter
    po::options_description desc("Allowed options");
    desc.add_options()
            ("output", po::value<string>(), "Specify the output directory for generating results")
            ("scene", po::value<string>(), "scene file")
            ("model", po::value<string>(), "model file to use")
            ("simplify", po::value<int>(), "Percentage of model to be rendered, e.g. 10 to indicate 10 percent of original model")
            ("propertyfile", po::value<string>(), "add property file containing key values")
            ("properties", po::value<string>(), "specify a list of key value pairs as 'key1:value1 key2:value2 ... keyn:valuen'")
            ("help", "Replaces variables in input file with properties specified in a properties file, or by properties provided on the command line. Result is written to output file ");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    // If no arguments specified, or help is requested, describe options for calling the program
    if (vm.size() == 0 || vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("output")) {
      outputDirectory = fs::path(vm["output"].as<string>());
    }

    if (vm.count("scene")) {
      sceneFilePath = fs::path(vm["scene"].as<string>());
      scenePropertiesFilePath = fs::path(sceneFilePath.string() + PROPERTY_SUFFIX);
    } else {
      cout << "Required scene file not provided. Please specify 'scene' argument and run again." << endl;
      return 1;
    }

    if (vm.count("propertyfile")) {
      userPropertiesFilePath = fs::path(vm["propertyfile"].as<string>());
    }

    if (vm.count("model")) {
      modelFilePath = fs::path(vm["model"].as<string>());
      modelPropertiesFilePath = modelFilePath.string() + PROPERTY_SUFFIX;
    }

    if (vm.count("simplify")) {
      if (!vm.count("model")) {
        cout << "Argument missing: Please specify a model to be simplified" << endl;
        return -1;
      }
      userRequestSimplifiedModel = true;
      simplifyPercentage = vm["simplify"].as<int>();
    }

    startProcessing();

    return 0;
}
