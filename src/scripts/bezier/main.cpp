/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "CurveRenderer.h"
#include "../bezier.h"
#include "../hairstruct.h"

using namespace std;

enum UserMode {
    RENDER_SAMPLE_BEZIER, RENDER_USER_INPUT, RENDER_HAIR_MODEL
};

const std::string hairFile = "/home/jeffrey/hair-rendering-thesis/scenes/hair/models/wStraight.10.pbrt";

const double bezier_points[] = {
    -5.0, 0.0, 0.0,
    0.0, 5.0, 0.0,
    5.0, 0.0, 0.0,
    5.0, 5.0, 0.0,
    5.0, 7.0, 5.0,
    7.0, 3.0, 5.0,
    0.0, 0.0, 0.0
};

CurveRenderer* renderer;
BezierSpline spline(bezier_points, 21);

const unsigned int WINDOW_WIDTH = 1024, WINDOW_HEIGHT = 1024;

void queryHairFile(std::string& fileName) {
    std::ifstream input;
    bool isInvalidFile = true;
    int numTries = 0;

    do {
        if (numTries < 3) {
            cout << "Enter hair model file name:\n";
            cin >> fileName;
        } else {
            cout << "Using default one: " << hairFile << endl;
            fileName = hairFile;
            break;
        }

        input.open(fileName.c_str());
        isInvalidFile = input.fail();
        if (isInvalidFile) {
            cout << "Cannot open specified file..." << endl;
            numTries++;
        }
    } while (isInvalidFile);

    if (input.fail()) {
        std::cout << "Giving up" << endl;
        exit(1);
    }

    Hair hair;
    input >> hair;

    cout << "Curve count for hair model: " << hair.fibers.size() << endl;
}

//void menu(UserMode& mode, std::string& fileName) {
//    bool isValidUserInput = false;
//    int number;
//
//    do {
//        cout << "Choose any of the options: " << endl;
//        cout << "\t1.) Use example Bezier curve" << endl
//                << "\t2.) Use hair model" << endl;
//
//
//        cin >> number;
//        isValidUserInput = number > 0 || number < 2;
//        if (!isValidUserInput) {
//            cout << "Invalid input ...";
//        }
//    } while (!isValidUserInput);
//
//    if (number == 1) {
//        mode = UserMode::RENDER_SAMPLE_BEZIER;
//        renderer->addCurve(spline);
//    } else {
//        mode = UserMode::RENDER_HAIR_MODEL;
//        queryHairFile(fileName);
//        renderer->addCurve(spline);
//    }
//}

int main(void) {
    renderer = new CurveRenderer(WINDOW_WIDTH, WINDOW_HEIGHT);
    cout << "Hello Bezier curve" << endl;

    UserMode mode;
    std::string fileName = hairFile;
    std::ifstream input(hairFile.c_str());
    if (input.fail()) {
        cout << "Failed to open hair file" << endl;
        return 0;
    }
    Hair hair;
    input >> hair;

    for (auto& pt : hair.fibers[0].curve.getControlPoints())
        cout << pt << endl;
    cout << "Curve count for hair model: " << hair.fibers.size() << endl;


    //menu(mode, fileName);

    renderer->addCurve(hair.fibers[0].curve);
    renderer->startup();
    return 0;
}