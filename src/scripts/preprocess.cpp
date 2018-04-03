#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>

#include "propertybag.h"
#include "hairsimplify.h"
#include "make_relative.h"

using namespace std;
namespace fs = boost::filesystem;


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


void createDirectoriesIfNotExists(fs::path& directory) {
  if (fs::create_directories(directory)) {
    cout << "* Succesfully created directory " << directory << endl;
  }
}

/** Returns the current time in a readable format */
string getCurrentTime() {
  namespace pt = boost::posix_time;
  return pt::to_simple_string(pt::second_clock::local_time());
}

/**
 * Parse properties from a property file in the specified properties collection
 */
void readProperties(const fs::path& propertyFilePath, PropertyBag& properties) {
  try {
      cout << "* Reading property file " << propertyFilePath.string() << "... ";
      properties.addPropertiesFromFile(propertyFilePath.string());
      cout << "[done]" << endl;
    } catch (string e) {
        cout << "\n\t" << e << endl;
    }
}

void generateVoxelGrid(const fs::path& inputModelPath, const fs::path& outputVoxelGridPath) {
    cout << "Generating voxel grid for " << inputModelPath << "\n  storing in " << outputVoxelGridPath << endl;
    cout << "[skipped] - Not yet supported" << endl << endl;
}


void simplifyAndWriteModel(const fs::path& inputModelPath, const fs::path& outputModelPath, int percentage) {
  if (fs::exists(outputModelPath)) {
    cout << "Already exists" << endl << endl;
  } else {
    cout << "Simplify model to " << percentage << " percent of original... " << endl;
    HairSimplify hair(inputModelPath.string());
    hair.reduceToPercentage(outputModelPath.string(), static_cast<double>(percentage));
  }
  
  fs::path outputVoxelGridPath = outputModelPath.string() + ".vdb";
  generateVoxelGrid(outputModelPath, outputVoxelGridPath);
}

/**
 * @brief startProcessing Sets up a test case, by copying property files from original and
 * possibly simplifying hair models
 */
void startProcessing(fs::path& baseOutputDirectory, const fs::path& inputScenePath, const fs::path& inputModelPath, int simplifyPercentage, bool shouldSimplify) {
    cout << "\nStart processing..." << endl << endl;

    PropertyBag sceneProperties, modelProperties, userProperties;

    // Create directories for output
    cout << "Creating directories for session output\n"
         << "------------------------------------------------" << endl;
    fs::path generatedContentDirectory = baseOutputDirectory / GENERATED_DIRECTORY_NAME;
    fs::path sessionResultsDirectory = baseOutputDirectory / getCurrentTime();
    fs::path propertyFilesDirectory = sessionResultsDirectory / PROPERTY_DIRECTORY_NAME;
    createDirectoriesIfNotExists(sessionResultsDirectory);
    createDirectoriesIfNotExists(propertyFilesDirectory);
    createDirectoriesIfNotExists(generatedContentDirectory);
    cout << "[done]" << endl << endl;

    // Read properties
    cout << "Reading property files (default attributes)\n"
         << "------------------------------------------------" << endl;
    fs::path inputScenePropertiesPath = inputScenePath.string() + ".properties";
    fs::path inputModelPropertiesPath = inputModelPath.string() + ".properties";
    readProperties(inputScenePropertiesPath, sceneProperties);
    readProperties(inputModelPropertiesPath, modelProperties);
    cout << "[done]" << endl << endl;

    // copy original scene file
    cout << "Copy original scene file... ";
    string sceneName = inputScenePath.filename().stem().string();
    sceneProperties.addProperty("scene_name", sceneName);
    fs::path outputScenePath = sessionResultsDirectory / (sceneName + ".pbrt");
    fs::copy_file(inputScenePath, sessionResultsDirectory / inputScenePath.filename());
    cout << "[done]" << endl << endl;

    // simplify model (or use original) and write to properties
    cout << "Simplifying model...\n"
         << "-----------------------------------------------" << endl;
    fs::path fromResultsToOriginalSceneDirectory = make_relative(sessionResultsDirectory, inputScenePath.parent_path());
    sceneProperties.addProperty("currentDirectory", fromResultsToOriginalSceneDirectory.string());

    string modelName = inputModelPath.filename().stem().string();
    modelProperties.addProperty("model_name", modelName);

    if (shouldSimplify)
    {
        string percentageAsStr = boost::lexical_cast<string>(simplifyPercentage);
        string outputModelFileName = modelName + "." + percentageAsStr + ".pbrt";
        fs::path outputModelPath = generatedContentDirectory / outputModelFileName;
        simplifyAndWriteModel(inputModelPath, outputModelPath, simplifyPercentage);

        fs::path relativeModelFilePath = make_relative(sessionResultsDirectory, outputModelPath);
        sceneProperties.addProperty("hairmodel_path", relativeModelFilePath.string());
        userProperties.addProperty("model_simplification_percentage", percentageAsStr);

        cout << "Using generated model at location: " << relativeModelFilePath << endl;
    }
    else
    {
        fs::path relativeModelFilePath = make_relative(sessionResultsDirectory, inputModelPath);
        sceneProperties.addProperty("hairmodel_path", relativeModelFilePath.string());
        cout << "No simplification needed. Using original model: " << relativeModelFilePath << endl;
    }
    cout << "[done]" << endl << endl;

    fs::path relativePathOriginalModel = make_relative(sessionResultsDirectory, inputModelPath);
    sceneProperties.addProperty("original_hairmodel_path", relativePathOriginalModel.string());

    // Write property files
    cout << "Writing property files...\n"
         << "-----------------------------------------------" << endl;
    fs::path outputScenePropertiesPath = sessionResultsDirectory / "properties" / "scene.properties";
    fs::path outputModelPropertiesPath = sessionResultsDirectory / "properties" / "model.properties";
    fs::path outputUserPropertiesPath = sessionResultsDirectory / "properties" / "user.properties";
    try {
        sceneProperties.writePropertiesToFile(outputScenePropertiesPath.string());
        modelProperties.writePropertiesToFile(outputModelPropertiesPath.string());
        userProperties.writePropertiesToFile(outputUserPropertiesPath.string());
        cout << "[done]" << endl << endl;
    } catch (string e) {
        cout << "\tERROR: " << e << endl << "[fail]" << endl << endl;
    }


    // generate a pbrt scene
    cout << "Generating PBRT scene " << outputScenePath << "...\n"
         << "-----------------------------------------------" << endl;
    PropertyBag properties;
    readProperties(outputScenePropertiesPath, properties);
    readProperties(outputModelPropertiesPath, properties);
    readProperties(outputUserPropertiesPath, properties);
    properties.replacePropertiesInFile(inputScenePath.string(), outputScenePath.string());
    cout << "[done]" << endl << endl;
}

string getProperty(string key, PropertyBag propertyBag) {
    try {
        return propertyBag.getProperty(key);
    } catch (string e) {
        cout << "ERROR: " << e << endl;
    }
}

string simplifyPath(string path) {
    if (boost::starts_with(path, ".")) {
        return path.substr(1);
    } else if (path.length() > 2 && boost::starts_with(path, "./")) {
        return path.substr(2);
    }
    return path;
}

bool isSimplificationPercentageMatching(const fs::path& modelPath, int percentage) {
    vector<string> parts;
    boost::split(parts, modelPath.filename().string(), boost::is_any_of("."));

    if (parts.size() < 3) {
        cout << "Not matching: " << modelPath << " does not contain 3 parts in the filename\n";
        return false;
    } else {
        int percentageDerivedFromFile = boost::lexical_cast<int>(parts[parts.size()-2]);
        cout << modelPath << " derived percentage is " << percentageDerivedFromFile << ", percentage is " << percentage << endl;
        return percentage == percentageDerivedFromFile;
    }
}

/**
 * @brief preprocessExisting Updates specified directory to make sure all changed properties are reflected
 * in scene file. Also it checks if model still exist and if not regenerate it.
 * @param directory The directory where the scene file is located and needs to be updated.
 */
void preprocessExisting(fs::path& directory) {
    directory = fs::path(simplifyPath(directory.string()));

    PropertyBag properties;
    cout << "Reading property files...\n"
         << "----------------------------------------" << endl;
    fs::path userPropertiesPath = directory / "properties" / "user.properties";
    fs::path scenePropertiesPath = directory / "properties" / "scene.properties";
    fs::path modelPropertiesPath = directory / "properties" / "model.properties";
    readProperties(scenePropertiesPath, properties);
    readProperties(modelPropertiesPath, properties);
    readProperties(userPropertiesPath, properties);
    cout << "[done]" << endl << endl;

    string sceneName = properties.getProperty("scene_name");
    fs::path inputScenePath = directory / (sceneName + ".scene");
    fs::path outputScenePath = directory / (sceneName + ".pbrt");

    fs::path originalModelPath = fs::path(properties.getProperty("original_hairmodel_path"));
    string modelName = properties.getProperty("model_name");
    fs::path modelPath = fs::path(properties.getProperty("hairmodel_path"));

    // check if model still exists, if not then recompute the model
    // also check if property file contains the same simplification percentage as derived by filename.
    // If not, then the user changed the percentage in the property file, indicating the model must be recomputed.
    if (!fs::exists(modelPath)) {
        // create generate folder
        fs::path modelDirectory = directory / modelPath.parent_path();
        createDirectoriesIfNotExists(modelDirectory);

        try {
            int simplifyPercentage = boost::lexical_cast<int>(properties.getProperty("model_simplification_percentage"));
            simplifyAndWriteModel(originalModelPath, directory / modelPath, simplifyPercentage);
        } catch (string e) {
            cout << "ERROR: " << e << endl
                 << "Something went wrong, the original model " << modelPath.filename() << " is not there anymore.\n"
                 << "Make sure you didn't delete it or change its file name. Nothing to restore\n";
            return;
        }

    }

    // write scene with replaced variables
    properties.replacePropertiesInFile(inputScenePath.string(), outputScenePath.string());
}


int main(int argc, char** argv) {
    namespace po = boost::program_options;

    // Input properties specified by user
    fs::path baseOutputDirectory("output"), outputDirectory;
    fs::path inputScenePath, outputScenePath;
    fs::path inputModelPath;
    fs::path inputScenePropertiesPath, inputUserPropertiesPath, inputModelPropertiesPath;
    bool isAskingForSimplification = false;
    int inputSimplificationPercentage;


    if (argc == 2) {
        fs::path directory(argv[1]);
        cout << "Preprocessing in directory: " << directory << endl;

        preprocessExisting(directory);
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
      baseOutputDirectory = fs::path(vm["output"].as<string>());
    }

    if (vm.count("scene")) {
      inputScenePath = fs::path(vm["scene"].as<string>());
      inputScenePropertiesPath = fs::path(inputScenePath.string() + PROPERTY_SUFFIX);
    } else {
      cout << "Required scene file not provided. Please specify 'scene' argument and run again." << endl;
      return 1;
    }

    if (vm.count("propertyfile")) {
      inputUserPropertiesPath = fs::path(vm["propertyfile"].as<string>());
    }

    if (vm.count("model")) {
      inputModelPath = fs::path(vm["model"].as<string>());
      inputModelPropertiesPath = inputModelPath.string() + PROPERTY_SUFFIX;
    }

    if (vm.count("simplify")) {
      if (!vm.count("model")) {
        cout << "Argument missing: Please specify a model to be simplified" << endl;
        return 1;
      }
      isAskingForSimplification = true;
      inputSimplificationPercentage = vm["simplify"].as<int>();
    }

    startProcessing(baseOutputDirectory, inputScenePath, inputModelPath, inputSimplificationPercentage, isAskingForSimplification);

    return 0;
}
