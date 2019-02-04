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

enum class UserMode {
    RENDER_SAMPLE_BEZIER, RENDER_USER_INPUT, RENDER_HAIR_MODEL
};

const std::string hairFile = "/home/jeffrey/hair-rendering-thesis/scenes/hair/models/wStraight.10.pbrt";

CurveRenderer* renderer;

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
//        cout << "1.) Use example Bezier curve" << endl
//             << "2.) Load hair model" << endl
//                << "3.) ";
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

int main(int argc, char** argv) {
    renderer = new CurveRenderer(WINDOW_WIDTH, WINDOW_HEIGHT);
    cout << "Hello Bezier Renderer" << endl;

    UserMode userMode;
    //std::map<std::string, std::string> params;
    //    handleMenu(userMode, params);
    //
    //    switch (userMode) {
    //        case UserMode::RENDER_HAIR_MODEL:
    //            break;
    //        case UserMode::RENDER_SAMPLE_BEZIER:
    //            break;
    //        case UserMode::RENDER_USER_INPUT:
    //            break;
    //    }


    std::string fileName = hairFile;
    std::ifstream input(hairFile.c_str());
    if (input.fail()) {
        cout << "Failed to open hair file" << endl;
        return 0;
    }
    Hair hair;
    input >> hair;
    //
    //    for (auto& pt : hair.fibers[0].curve.getControlPoints())
    //        cout << pt << endl;
    //    cout << "Curve count for hair model: " << hair.fibers.size() << endl;


    //menu(mode, fileName);


    double pt[36] = {-25.045174, 50.892357, 18.287086, -25.260208, 51.178432, 18.687077, -25.510658, 51.287296, 19.049347, -25.796530, 51.218952, 19.373894,
        -25.796530, 51.218952, 19.373894, -26.082401, 51.150608, 19.698441, -26.423557, 50.851059, 19.976362, -26.760401, 50.482285, 20.234373,
        -26.760401, 50.482285, 20.234373, -27.097244, 50.113510, 20.492384, -27.469509, 49.567532, 20.712675, -27.817591, 49.006298, 20.921959};

    BezierSpline bezier_all(pt, 24);
    bezier_all.setUseSharedControlPoints(false);

    double pt1[12] = {-25.045174, 50.892357, 18.287086, -25.260208, 51.178432, 18.687077, -25.510658, 51.287296, 19.049347, -25.796530, 51.218952, 19.373894};
    double pt2[12] = {-25.796530, 51.218952, 19.373894, -26.082401, 51.150608, 19.698441, -26.423557, 50.851059, 19.976362, -26.760401, 50.482285, 20.234373};
    double pt3[12] = {-26.760401, 50.482285, 20.234373, -27.097244, 50.113510, 20.492384, -27.469509, 49.567532, 20.712675, -27.817591, 49.006298, 20.921959};
    BezierSpline bezier_pt1(pt1, 12);
    BezierSpline bezier_pt2(pt2, 12);
    BezierSpline bezier_pt3(pt3, 12);
    bezier_pt1.setUseSharedControlPoints(false);
    bezier_pt2.setUseSharedControlPoints(false);
    bezier_pt3.setUseSharedControlPoints(false);

    cout << "Control points:" << endl;
    for (auto cp : bezier_all.getControlPoints()) {
        cout << cp << endl;
    }
    cout << endl;
    cout << "Sample points" << endl;
    cout << "sampleCurve(" << 0.5f << ") = " << bezier_all.sampleCurve(0.5f) << endl;
    cout << "sampleSegment (0, 1.0) = " << bezier_all.sampleSegment(0, 1.0f) << endl;
    cout << "sampleSegment (1, 0.0) = " << bezier_all.sampleSegment(1, 0.0f) << endl;
    //    for (float t = 0.0; t <= 1.01; t += 0.1) {
    //        cout << "bezier A(" << t << ") = " << bezier_all.sampleCurve(t) << endl;
    //    }

    //    BezierSpline bs1(
    //            Point3(0.0, 5.0, 0.0),
    //            Point3(-1.0, 5.1, 0.0),
    //            Point3(1.0, 5.1, 0.0),
    //            Point3(0.0, 5.0, 0.0));
    //    BezierSpline bs2(
    //            Point3(0.0, 10.0, 0.0),
    //            Point3(-1.0, 10.1, 0.0),
    //            Point3(1.0, 10.0, 0.0),
    //            Point3(0.0, 10.0, 0.0));
    //    BezierSpline bs3(
    //            Point3(0.0, 20.0, 0.0),
    //            Point3(-1.0, 20.0, 0.0),
    //            Point3(1.0, 20.0, 0.0),
    //            Point3(0.0, 20.0, 0.0));
    //    bs1.setUseSharedControlPoints(false);
    //    bs2.setUseSharedControlPoints(false);
    //    bs3.setUseSharedControlPoints(false);

    cout << "Curves before optimize: " << hair.fibers.size() << endl;
    hair.optimizeCurves();
    cout << "Curves after optimize: " << hair.fibers.size() << endl;
    renderer->addHairModel(hair, 100);

    //    renderer->addCurve(bezier_pt1);
    //    renderer->addCurve(bezier_pt2);
    //    renderer->addCurve(bezier_pt3);
    renderer->startup();

    return 0;
}