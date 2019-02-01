/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   CurveRenderer.cpp
 * Author: jeffrey
 *
 * Created on January 31, 2019, 8:47 PM
 */

#include "CurveRenderer.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
//#include "plyloader.h"
//#include "ShaderSource.h"
//#include "SimpleGlUtil.h"
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/matrix.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/rotate_normalized_axis.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <fstream>

#include "../gl_util/ShaderSource.h"
#include "../gl_util/SimpleGlUtil.h"
#include "../bezier.h"
#include "../hairstruct.h"
using namespace std;

GLFWwindow* window = nullptr;
const int BEZIER_CURVE = 0;
const int CONTROL_POINTS = 1;
GLuint vbo[2];
GLuint shaderProgram[2];

//BezierSpline spline(bezier_points, 21);
//std::vector<double> curve;

const unsigned int WINDOW_WIDTH = 1024, WINDOW_HEIGHT = 1024;
float aspect = WINDOW_WIDTH / static_cast<float> (WINDOW_HEIGHT);
float fovy = 120;
float zNear = 0.1f, zFar = 1000.0f;

glm::vec3 red(1.0f, 0.0f, 0.0f);
glm::vec3 green(0.0f, 1.0f, 0.0f);
glm::vec3 blue(0.0f, 1.0f, 1.0f);
glm::vec3 colors[3] = {red, green, blue};

glm::vec3 eye(0.0f, 0.0f, 7.0f);
glm::vec3 center(.0f, 2.5f, .0f);
glm::vec3 up(.0f, 1.0f, .0f);

glm::mat4 model = glm::mat4(1.0f);
glm::mat4 view = glm::lookAt(eye, center, up);
glm::mat4 perspective = glm::perspective(fovy, aspect, zNear, zFar);

std::vector<double> dataPoints;
std::vector<int> curveOffsets;

void onScrolled(GLFWwindow* window, double xOffset, double yOffset) {
    eye *= yOffset * -0.02 + 1.0;
    view = glm::lookAt(eye, center, up);
    SimpleGlUtil::setUniform(shaderProgram[BEZIER_CURVE], "view", view);
    SimpleGlUtil::setUniform(shaderProgram[CONTROL_POINTS], "view", view);
}

void onKey(GLFWwindow* window, int key, int scancode, int action, int mods) {
    glm::quat rotation;

    if (action == GLFW_PRESS) {
        switch (key) {
            case GLFW_KEY_UP:
                rotation = glm::angleAxis(0.1f, glm::normalize(glm::cross(eye, up)));
                break;
            case GLFW_KEY_DOWN:
                rotation = glm::angleAxis(-0.1f, glm::normalize(glm::cross(eye, up)));
                break;
            case GLFW_KEY_LEFT:
                rotation = glm::angleAxis(0.1f, up);
                break;
            case GLFW_KEY_RIGHT:
                rotation = glm::angleAxis(-0.1f, up);

                break;
            default:
                return;
        }
        eye = glm::rotate(rotation, eye);
        view = glm::lookAt(eye, center, up);
        SimpleGlUtil::setUniform(shaderProgram[BEZIER_CURVE], "view", view);
        SimpleGlUtil::setUniform(shaderProgram[CONTROL_POINTS], "view", view);
    }
}

CurveRenderer::CurveRenderer(unsigned int windowWidth, unsigned int windowHeight)
: WINDOW_WIDTH(windowWidth), WINDOW_HEIGHT(windowHeight) {
}

CurveRenderer::CurveRenderer(const CurveRenderer& orig) : CurveRenderer() {
}

CurveRenderer::~CurveRenderer() {
}

void CurveRenderer::addCurve(const BezierSpline& curve) {
    mCurves.push_back(&curve);
}

void CurveRenderer::addHairModel(const Hair& hair) {
    //    std::vector<double> mCoords;
    //
    //    for (auto& fiber : hair.fibers) {
    //        for (auto& pt : fiber.curve.getControlPoints()) {
    //            mCoords.push_back(pt.x);
    //            mCoords.push_back(pt.y);
    //            mCoords.push_back(pt.z);
    //        }
    //    }
}

void CurveRenderer::init() {
    glfwSetScrollCallback(window, onScrolled);
    glfwSetKeyCallback(window, onKey);

    glEnable(GL_DOUBLEBUFFER | GL_DEPTH_BUFFER);
    //glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const string assetFolder = "/home/jeffrey/hair-rendering-thesis/assets/";
    cout << "Reading shader files ...\n";
    shaderProgram[BEZIER_CURVE] = SimpleGlUtil::createShaderProgram(
            assetFolder + "bezier-vs.glsl",
            assetFolder + "bezier-fs.glsl");
    shaderProgram[CONTROL_POINTS] = SimpleGlUtil::createShaderProgram(
            assetFolder + "bezier-vs.glsl",
            assetFolder + "bezier-fs.glsl");

    // allocate buffer for all curves
    int nSamples = 10;

    curveOffsets.push_back(0);

    for (auto curve : mCurves) {
        for (int i = 0; i < nSamples; ++i) {
            double t = i / static_cast<double> (nSamples - 1);

            Point3 pt = curve->sampleCurve(t);
            dataPoints.push_back(pt.x);
            dataPoints.push_back(pt.y);
            dataPoints.push_back(pt.z);
        }
        std::cout << "Offset added: " << dataPoints.size() << std::endl;
        curveOffsets.push_back(dataPoints.size());
    }



    glGenBuffers(2, vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo[BEZIER_CURVE]);
    glBufferData(GL_ARRAY_BUFFER, dataPoints.size() * sizeof (double), (void*) &dataPoints[0], GL_STATIC_DRAW);
    //    glBindBuffer(GL_ARRAY_BUFFER, vbo[CONTROL_POINTS]);
    //    glBufferData(GL_ARRAY_BUFFER, mCurves[0]->getControlPointCount() * sizeof (double) * 3, (const void*) &(mCurves[0]->getControlPoints()[0]), GL_STATIC_DRAW);
}

void CurveRenderer::render() {
    // Render bezier
    glUseProgram(shaderProgram[BEZIER_CURVE]);
    SimpleGlUtil::setUniform(shaderProgram[BEZIER_CURVE], "view", view);
    SimpleGlUtil::setUniform(shaderProgram[BEZIER_CURVE], "perspective", perspective);
    SimpleGlUtil::setUniform(shaderProgram[BEZIER_CURVE], "model", model);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vbo[BEZIER_CURVE]);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0 * sizeof (GLdouble), (void*) (0 * sizeof (GLdouble)));

    for (int offsetIndex = 0; offsetIndex < curveOffsets.size() - 1; ++offsetIndex) {
        SimpleGlUtil::setUniform(shaderProgram[BEZIER_CURVE], "color", colors[offsetIndex]);
        glBindBuffer(GL_ARRAY_BUFFER, vbo[BEZIER_CURVE]);
        int start = curveOffsets[offsetIndex] / 3;
        int end = (curveOffsets[offsetIndex + 1] - curveOffsets[offsetIndex]) / 3; //(curveOffsets[offsetIndex + 1] - start) / 3;
        glDrawArrays(GL_LINE_STRIP, start, end);
    }

    //Render control points
    //    glUseProgram(shaderProgram[CONTROL_POINTS]);
    //    SimpleGlUtil::setUniform(shaderProgram[CONTROL_POINTS], "view", view);
    //    SimpleGlUtil::setUniform(shaderProgram[CONTROL_POINTS], "perspective", perspective);
    //    SimpleGlUtil::setUniform(shaderProgram[CONTROL_POINTS], "model", model);
    //
    //    glEnableVertexAttribArray(0);
    //    glBindBuffer(GL_ARRAY_BUFFER, vbo[CONTROL_POINTS]);
    //    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0 * sizeof (GLdouble), (void*) (0 * sizeof (GLdouble)));
    //    SimpleGlUtil::setUniform(shaderProgram[CONTROL_POINTS], "color", green);
    //    glDrawArrays(GL_LINE_STRIP, 0, mCurves[0]->getControlPointCount() / sizeof (double));
}

int CurveRenderer::startup(void) {

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    glfwWindowHint(GLFW_SAMPLES, 16);
    window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Hello Bezier", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    if (GLEW_OK != glewInit()) {
        cout << "GLEW Error" << endl;
        return -1;
    }

    init();

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        render();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}