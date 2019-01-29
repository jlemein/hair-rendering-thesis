/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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
#include <glm/gtc/type_ptr.hpp>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>

#include "../gl_util/ShaderSource.h"
#include "../gl_util/SimpleGlUtil.h"
#include "../bezier.h"

//#include "../bezier.h";

using namespace std;

GLFWwindow* window = nullptr;

const double bezier_points[] = {
    -5.0, 0.0, 0.0,
    0.0, 5.0, 0.0,
    5.0, 0.0, 0.0,
    5.0, 5.0, 0.0
};

GLuint vbo;
GLuint shaderProgram;
BezierSpline spline(bezier_points, 12);
std::vector<double> curve;

const unsigned int WINDOW_WIDTH = 1024, WINDOW_HEIGHT = 1024;
float aspect = WINDOW_WIDTH / static_cast<float> (WINDOW_HEIGHT);
float fovy = 120;
float zNear = 0.1f, zFar = 1000.0f;

glm::vec3 eye(0.0f, 0.0f, 7.0f);
glm::vec3 center(.0f, 2.5f, .0f);
glm::vec3 up(.0f, 1.0f, .0f);

glm::mat4 model = glm::mat4(1.0f);
glm::mat4 view = glm::lookAt(eye, center, up);
glm::mat4 perspective = glm::perspective(fovy, aspect, zNear, zFar);

void onScrolled(GLFWwindow* window, double xOffset, double yOffset) {
    eye *= yOffset * -0.02 + 1.0;
    view = glm::lookAt(eye, center, up);
    SimpleGlUtil::setUniform(shaderProgram, "view", view);
}

void init() {
    glfwSetScrollCallback(window, onScrolled);

    glEnable(GL_DOUBLEBUFFER | GL_DEPTH_BUFFER);
    //glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const string assetFolder = "/home/jeffrey/hair-rendering-thesis/assets/";
    const string vsFile = assetFolder + "vertexShader.glsl";
    const string fsFile = assetFolder + "fragmentShader.glsl";

    cout << "Reading shader files ...\n";



    for (double t = 0.0; t < 1.0; t += 0.05) {
        Point3 pt = spline.sample(t);
        std::cout << "t = " << t << "  [" << pt.x << ", " << pt.y << ", " << pt.z << "]" << std::endl;
        curve.push_back(pt.x);
        curve.push_back(pt.y);
        curve.push_back(pt.z);
    }

    shaderProgram = SimpleGlUtil::createShaderProgram(vsFile, fsFile);


    /* Create plane buffers */
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, curve.size() * sizeof (double), (void*) &curve[0], GL_STATIC_DRAW);

    glUseProgram(shaderProgram);
}

void renderControlPoints() {
    //glUseProgram(controlPointShader);
}

void render() {
    // activate shader and set uniforms
    glUseProgram(shaderProgram);
    SimpleGlUtil::setUniform(shaderProgram, "view", view);
    SimpleGlUtil::setUniform(shaderProgram, "perspective", perspective);
    SimpleGlUtil::setUniform(shaderProgram, "model", model);
    //SimpleGlUtil::setUniform(shaderProgram, "color", color);

    // enable the vertex attribute arrays
    glEnableVertexAttribArray(0);
    //glEnableVertexAttribArray(1);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0 * sizeof (GLdouble), (void*) (0 * sizeof (GLdouble)));
    //glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 6 * sizeof (GLfloat), (void*) (3 * sizeof (GLfloat)));

    /* Render triangles */
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    //glDrawArrays(GL_TRIANGLES, 0, 3);
    glDrawArrays(GL_LINE_STRIP, 0, curve.size() / 3);
}

int main(void) {
    cout << "Hello Bezier curve" << endl;

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