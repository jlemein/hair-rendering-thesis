/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShaderSource.h"
#include <fstream>
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>
using namespace std;

void ShaderSource::Source(const std::string& fileName, std::vector<char*>& sourceLines, int& lineCount) {
    ifstream input(fileName.c_str());
    if (input.fail()) {
        cout << "Cannot read shader source from: " << fileName << endl;
        lineCount = 0;
        return;
    }
    sourceLines.clear();

    string source = "";
    string s;
    while (!input.eof()) {
        getline(input, s);
        s += "\n";
        char* line = new char[s.length() + 1];
        strcpy(line, s.c_str());
        sourceLines.push_back(line);
    }

    input.close();
    lineCount = sourceLines.size();
    return;
}
