/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "CurveRenderer.h"
#include "../bezier.h"
#include "../hairstruct.h"

using namespace std;

static CurveRenderer* renderer = 0;

const std::string DEFAULT_HAIR_PATH = "/home/jeffrey/hair-rendering-thesis/scenes/hair/models/wStraight.10.pbrt";

const unsigned int WINDOW_WIDTH = 1024, WINDOW_HEIGHT = 1024;

void queryHairFile(Hair& hairModel, const string& defaultHairFile = "") {
    ifstream inputStream;
    bool isInvalidFile = true;
    string fileName;
    int numTries = 0;

    do {
        cout << "Enter hair model file name:\n(Leave empty to use default hair file: '" << defaultHairFile << "')\n";
        getline(cin, fileName);

        if (fileName.empty()) {
            cout << "Using default hair model" << endl;
            fileName = defaultHairFile;
        }

        cout << "Opening " << fileName << endl;
        inputStream.open(fileName.c_str());
        isInvalidFile = inputStream.fail();
        if (isInvalidFile) {
            cout << "Cannot open specified file..." << endl;
            numTries++;
        }
    } while (isInvalidFile && numTries < 3);

    if (inputStream.fail()) {
        std::cout << "Giving up" << endl;
        exit(1);
    }

    inputStream >> hairModel;

    cout << "Curve count for hair model: " << hairModel.fibers.size() << endl;
}

static void actionUseExampleCurve() {
    double pt[36] = {-25.045174, 50.892357, 18.287086, -25.260208, 51.178432, 18.687077, -25.510658, 51.287296, 19.049347, -25.796530, 51.218952, 19.373894,
        -25.796530, 51.218952, 19.373894, -26.082401, 51.150608, 19.698441, -26.423557, 50.851059, 19.976362, -26.760401, 50.482285, 20.234373,
        -26.760401, 50.482285, 20.234373, -27.097244, 50.113510, 20.492384, -27.469509, 49.567532, 20.712675, -27.817591, 49.006298, 20.921959};

    BezierSpline spline(pt, 24);
    spline.setUseSharedControlPoints(false);
    renderer->addCurve(spline);
}

static void actionSpecifyControlPoints() {
    vector<Point3> controlPoints;
    string input;
    Point3 controlPoint;

    do {
        cout << "Control point #" << controlPoints.size() << " (x y z): ";
        getline(cin, input);
        if (input.empty()) {
            break;
        }

        stringstream istream(input);
        istream >> controlPoint.x >> controlPoint.y >> controlPoint.z;
        controlPoints.push_back(controlPoint);
    } while (true);

    char answer;
    cout << "Does your Bezier curve uses shared control points (y = yes): ";
    cin >> answer;

    BezierSpline spline;
    spline.addControlPoints(controlPoints);
    spline.setUseSharedControlPoints(answer == 'y');
    renderer->addCurve(spline);
}

static void ActionRenderBezierCurve() {
    std::cout << "Enter your choice:" << endl;
    std::cout << "1.) Specify control points" << endl;
    std::cout << "2.) Use example curve" << endl;

    int choice;
    cin >> choice;

    if (choice == 1) {
        actionSpecifyControlPoints();
    } else if (choice == 2) {
        actionUseExampleCurve();
    }
}

static void ActionRenderHairModel() {
    Hair hairModel;
    queryHairFile(hairModel, DEFAULT_HAIR_PATH);

    // add only 100 hair strands to limit memory
    hairModel.optimizeCurves();
    renderer->addHairModel(hairModel, 100);
}

static void ShowMenu() {
    bool isValidUserInput = false;
    string input;

    do {
        cout << "Enter your choice (enter a number): " << endl;
        cout << "1.) Render Bezier curve" << endl
                << "2.) Render hair model" << endl;

        getline(cin, input);
        isValidUserInput = input == "1" || input == "2";
        if (!isValidUserInput) {
            cout << "Invalid input ...";
        }
    } while (!isValidUserInput);

    if (input == "1") {
        ActionRenderBezierCurve();
    } else {
        ActionRenderHairModel();
    }
}

int main(int argc, char** argv) {
    renderer = new CurveRenderer(WINDOW_WIDTH, WINDOW_HEIGHT);
    cout << "Hello Bezier Renderer" << endl;

    //ShowMenu();

    Hair hair;
    ifstream input(DEFAULT_HAIR_PATH);
    input >> hair;
    hair.optimizeCurves();

    renderer->addHairModel(hair, 100);
    renderer->startup();

    return 0;
}