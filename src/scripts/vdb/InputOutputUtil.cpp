/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   InputOutputUtil.cpp
 * Author: jeffrey
 *
 * Created on February 1, 2019, 10:08 AM
 */

#include "InputOutputUtil.h"

#include <iostream>

const std::string InputOutputUtil::DEFAULT_FILE_PATH = "./voxelgrid.vdb";

InputOutputUtil::InputOutputUtil() {
}

InputOutputUtil::InputOutputUtil(const InputOutputUtil& orig) {
}

InputOutputUtil::~InputOutputUtil() {
}

bool InputOutputUtil::QueryInputFileUntilValid(std::ifstream& input) {
    std::cout << "Default file input path: " << DEFAULT_FILE_PATH << std::endl;
    int failedAttempts = 0;

    do {
        std::string vdbFileName;
        std::cout << "Enter path to input file (leave empty to use default):\n>";
        getline(std::cin, vdbFileName);
        vdbFileName = vdbFileName.empty() ? DEFAULT_FILE_PATH : vdbFileName;

        input.open(vdbFileName);
        if (input.fail()) {
            std::cout << "Could not open file '" << vdbFileName << "'" << std::endl;
            failedAttempts++;
            continue;
        }
        return true;
    } while (failedAttempts < 3);

    return false;
}

bool InputOutputUtil::QueryOutputFileUntilValid(std::ofstream& output) {
    std::cout << "Default file output path: " << DEFAULT_FILE_PATH << std::endl;
    int failedAttempts = 0;

    do {
        std::string vdbFileName;
        std::cout << "Enter path to output file (leave empty to use default):\n>";
        getline(std::cin, vdbFileName);
        vdbFileName = vdbFileName.empty() ? DEFAULT_FILE_PATH : vdbFileName;

        output.open(vdbFileName);
        if (output.fail()) {
            std::cout << "Could not open file '" << vdbFileName << "' for writing" << std::endl;
            failedAttempts++;
            continue;
        }
        return true;
    } while (failedAttempts < 3);

    return false;
}

bool InputOutputUtil::OpenFile(std::ifstream& input, const std::string& pathToFile) {
    if (!pathToFile.empty()) {
        input.open(pathToFile.c_str());
        if (!input.fail()) {
            return true;
        } else {
            std::cout << "Could not open '" << pathToFile << "'" << std::endl;
        }
    }

    return QueryInputFileUntilValid(input);
}

bool InputOutputUtil::OpenFile(std::ofstream& output, const std::string& pathToFile) {
    if (!pathToFile.empty()) {
        output.open(pathToFile.c_str());
        if (!output.fail()) {
            return true;
        } else {
            std::cout << "Could not open '" << pathToFile << "' for writing" << std::endl;
        }
    }

    return QueryOutputFileUntilValid(output);
}