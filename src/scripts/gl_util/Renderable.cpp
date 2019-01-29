/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Renderable.h"
#include <GL/glew.h>
#include <iostream>

Renderable::Renderable(const std::string& fileName) {
    PlyLoader* loader = new PlyLoader("assets/icosphere.ply");
    loader->parse();

    /* Set up vertex buffer */
    glGenBuffers(1, &mVbo);
    glGenBuffers(1, &mIbo);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    for (int i = 0; i < 12; ++i) {
        std::cout << i << ": " << loader->getVertices()[i] << std::endl;
    }
    glBufferData(GL_ARRAY_BUFFER, loader->getVertexCount() * 6 * sizeof (GLfloat), (void*) loader->getVertices(), GL_STATIC_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, 9 * sizeof (GLfloat), (void*) vertices, GL_STATIC_DRAW);


    //glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3, 0); // for triangle
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof (GLfloat), (void*) (0 * sizeof (GLfloat)));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof (GLfloat), (void*) (3 * sizeof (GLfloat)));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, loader->getIndexCount() *3 * sizeof (unsigned short), (void*) loader->getIndices(), GL_STATIC_DRAW);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * sizeof (GLushort), (void*) indices, GL_STATIC_DRAW);

    /* Set up render program */
    program = glCreateProgram();
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    compileShader("assets/blinn_vs.glsl", vertexShader);
    compileShader("assets/blinn_fs.glsl", fragmentShader);

    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    linkProgram(program);
}