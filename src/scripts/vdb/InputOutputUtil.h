/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputOutputUtil.h
 * Author: jeffrey
 *
 * Created on February 1, 2019, 10:08 AM
 */

#ifndef INPUTOUTPUTUTIL_H
#define INPUTOUTPUTUTIL_H

#include <string>
#include <fstream>

class InputOutputUtil {
public:
    static bool QueryInputFileUntilValid(std::ifstream& ifstream);
    static bool QueryOutputFileUntilValid(std::ofstream& ofstream);
    static bool OpenFile(std::ifstream& ifstream, const std::string& pathToFile = "");
    static bool OpenFile(std::ofstream& ofstream, const std::string& pathToFile = "");
    
private:
    const static std::string DEFAULT_FILE_PATH;
    
    InputOutputUtil();
    InputOutputUtil(const InputOutputUtil& orig);
    virtual ~InputOutputUtil();
};

#endif /* INPUTOUTPUTUTIL_H */

