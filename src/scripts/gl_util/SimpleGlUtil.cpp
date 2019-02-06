/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#ifndef __APPLE__
#include <GL/glew.h>
#else
#include <OpenGL/gl3.h> //OS x libs
#endif

#include "SimpleGlUtil.h"
#include <iostream>
#include "ShaderSource.h"
#include <vector>
#include <sstream>
#include <string>

#include <glm/matrix.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
using namespace std;

static GLuint createShaderProgram(const std::string& vsName, const std::string& fsName);

void compileShader(const std::string& fileName, unsigned int shader) {
    std::vector<char*> sourceLines;
    int lineCount;

    ShaderSource::Source(fileName, sourceLines, lineCount);
    glShaderSource(shader, lineCount, &sourceLines[0], NULL);
    glCompileShader(shader);

    GLint compiled = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (compiled) {
        cout << "Shader '" + fileName + "' OK" << endl;
    } else {
        cout << "Shader '" << fileName << "' failed to compile:" << endl;
        GLint logSize = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logSize);
        GLchar* infoLog = new GLchar[logSize];
        glGetShaderInfoLog(shader, logSize, NULL, infoLog);
        cout << infoLog << endl;
        delete[] infoLog;
    }
}

void linkProgram(int program) {
    glLinkProgram(program);
    GLint isLinked = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &isLinked);
    if (isLinked == GL_FALSE) {
        GLint maxLength = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);

        cout << "Failed to link" << endl;
        cout << "Info log size: " << infoLog.size() << endl;
        for (auto line : infoLog) {
            cout << line;
        }

        // The program is useless now. So delete it.
        glDeleteProgram(program);
    } else {
        cout << "Linking OK" << endl;
    }
}

bool SimpleGlUtil::isShadingLanguageSupported(int major, int minor) {
    std::string glslVersion = string((const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
    std::istringstream stream(glslVersion);
    std::string token;
    
    getline(stream, token, '.');
    int supportedMajor = std::stoi(token);
    
    getline(stream, token, '.');
    int supportedMinor = std::stoi(token);
    
    return !(supportedMajor < major || (supportedMajor == major && supportedMinor < 3));
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, const glm::mat4& value) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(value));
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, const glm::mat3& value) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniformMatrix3fv(location, 1, GL_FALSE, glm::value_ptr(value));
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, const glm::vec3& value) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniform3fv(location, 1, glm::value_ptr(value));
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, const glm::vec4& value) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniform4fv(location, 1, glm::value_ptr(value));
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, int value) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniform1i(location, value);
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, float value) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniform1f(location, value);
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, const glm::vec3 values[], unsigned int length) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniform3fv(location, length, glm::value_ptr(values[0]));
}

void SimpleGlUtil::setUniform(GLuint program, const std::string& param, const glm::vec4 values[], unsigned int length) {
    GLint location = glGetUniformLocation(program, param.c_str());
    glUniform4fv(location, length, glm::value_ptr(values[0]));
}

GLuint SimpleGlUtil::createShaderProgram(const std::string& vsName, const std::string& fsName) {
    /* Set up render program */
    GLuint program = glCreateProgram();
    GLint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    compileShader(vsName, vertexShader);
    compileShader(fsName, fragmentShader);

    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    linkProgram(program);

    return program;
}

GLuint SimpleGlUtil::createShaderProgram(const std::string& vsName, const std::string& gsName, const std::string& fsName) {
    /* Set up render program */
    GLuint program = glCreateProgram();
    GLint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLint geometryShader = glCreateShader(GL_GEOMETRY_SHADER);
    GLint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    compileShader(vsName, vertexShader);
    compileShader(gsName, geometryShader);
    compileShader(fsName, fragmentShader);

    glAttachShader(program, vertexShader);
    glAttachShader(program, geometryShader);
    glAttachShader(program, fragmentShader);
    linkProgram(program);

    return program;
}

RenderState SimpleGlUtil::createBufferWithData(const GLfloat* vertices, int verticesSizeBytes, const GLushort* indices, int indicesSizeBytes) {
    RenderState state;

    glGenBuffers(1, &state.vbo);
    glGenBuffers(1, &state.ibo);
    glBindBuffer(GL_ARRAY_BUFFER, state.vbo);
    glBufferData(GL_ARRAY_BUFFER, verticesSizeBytes, vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state.ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indicesSizeBytes, indices, GL_STATIC_DRAW);

    return state;
}

RenderState SimpleGlUtil::createBufferWithData(const GLfloat* vertices, int vertexSize, int attributesPerVertex, const GLushort* indices, int indexCount) {
    RenderState state;

    glGenBuffers(1, &state.vbo);
    glGenBuffers(1, &state.ibo);
    glBindBuffer(GL_ARRAY_BUFFER, state.vbo);
    glBufferData(GL_ARRAY_BUFFER, vertexSize * attributesPerVertex * sizeof (GLfloat), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state.ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof (GLushort) * indexCount, indices, GL_STATIC_DRAW);

    return state;
}

void SimpleGlUtil::setRenderState(const RenderState& rs) {
    glBindBuffer(GL_ARRAY_BUFFER, rs.vbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rs.ibo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof (GLfloat), (void*) (0 * sizeof (GLfloat)));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof (GLfloat), (void*) (3 * sizeof (GLfloat)));
}
