/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GlUtil.h
 * Author: jeffrey
 *
 * Created on October 11, 2018, 11:53 PM
 */

#ifndef GLUTIL_H
#define GLUTIL_H

#include <string>
#include <glm/glm.hpp>

struct RenderState {
    GLuint vbo, ibo;
    unsigned int vertexCount, indexCount;
   // GLfloat* vertices;
   // GLushort* indices;
};

class SimpleGlUtil {
public:
    static RenderState createBufferWithData(const GLfloat* vertices, int verticesSizeBytes, const GLushort* indices, int indicesSizeBytes);
    static RenderState createBufferWithData(const GLfloat* vertices, int vertexCount, int attributeCount, const GLushort* indices, int indexCount);
    static void setRenderState(const RenderState& rs);
    
    
    /**
    * Returns program
    */
    static GLuint createShaderProgram(const std::string& vsName, const std::string& fsName);
    static GLuint createShaderProgram(const std::string& vsName, const std::string& gsName, const std::string& fsName);
    //static int compileShader(const std::string& vsName, const std::string& fsName);
    
    /**
     * Checks if GLSL is supported
     */
    static bool isShadingLanguageSupported(int major, int minor = 0);
    
    /**
     * Sets the uniform parameter of the program
     * @param param Parameter name
     * @param value Value of param
     * @return True when param exists and could be set
     */
    static void setUniform(GLuint program, const std::string& param, const glm::mat4& value);
    static void setUniform(GLuint program, const std::string& param, const glm::mat3& value);
    static void setUniform(GLuint program, const std::string& param, const glm::vec3& value);
    static void setUniform(GLuint program, const std::string& param, const glm::vec4& value);
    static void setUniform(GLuint program, const std::string& param, int value);
    static void setUniform(GLuint program, const std::string& param, float value);
    static void setUniform(GLuint program, const std::string& param, const glm::vec3 value[], unsigned int length);
    static void setUniform(GLuint program, const std::string& param, const glm::vec4 value[], unsigned int length);
    
};

#endif /* GLUTIL_H */

