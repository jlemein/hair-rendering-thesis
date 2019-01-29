/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Renderable.h
 * Author: jeffrey
 *
 * Created on October 11, 2018, 9:32 PM
 */

#ifndef RENDERABLE_H
#define RENDERABLE_H

#include <memory>
#include <GL/glew.h>

struct VertexInfo {
    int stride;
    
};
class Renderable {
private:
    GLint mVbo, mIbo;
    GLint mProgram;
    
public:
    Renderable(const std::string& fileName);
    Renderable(std::shared_ptr<GLfloat*> vertices, shared_ptr<GLushort*> indices, VertexInfo& info);
    void render(GLint program);
};

#endif /* RENDERABLE_H */

