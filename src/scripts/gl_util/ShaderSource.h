/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ShaderSource.h
 * Author: jeffrey
 *
 * Created on October 10, 2018, 9:49 PM
 */

#ifndef SHADERSOURCE_H
#define SHADERSOURCE_H

#include <string>
#include <vector>

class ShaderSource {
public:
    static void Source(const std::string& fileName, std::vector<char*>& sourceLines, int& lineCount);
};

#endif /* SHADERSOURCE_H */

