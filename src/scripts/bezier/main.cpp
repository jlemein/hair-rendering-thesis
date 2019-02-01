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


    double pt[36] = {-25.045174, 50.892357, 18.287086, -25.260208, 51.178432, 18.687077, -25.510658, 51.287296, 19.049347, -25.796530, 51.218952, 19.373894,
        -25.796530, 51.218952, 19.373894, -26.082401, 51.150608, 19.698441, -26.423557, 50.851059, 19.976362, -26.760401, 50.482285, 20.234373,
        -26.760401, 50.482285, 20.234373, -27.097244, 50.113510, 20.492384, -27.469509, 49.567532, 20.712675, -27.817591, 49.006298, 20.921959};

    BezierSpline bezier_all(pt, 36);
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

    //renderer->addHairModel(hair);
    renderer->addCurve(bezier_pt1);
    renderer->addCurve(bezier_pt2);
    renderer->addCurve(bezier_pt3);
    //    renderer->addCurve(bs2);
    //renderer->addCurve(bs3);
    //renderer->addCurve(hair.fibers[0].curve);
    renderer->startup();
    return 0;
}