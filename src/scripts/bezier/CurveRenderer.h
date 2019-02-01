/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CurveRenderer.h
 * Author: jeffrey
 *
 * Created on January 31, 2019, 8:47 PM
 */

#ifndef CURVERENDERER_H
#define CURVERENDERER_H

#include <vector>
class BezierSpline;
class Hair;

class CurveRenderer {
public:
    CurveRenderer(unsigned int width = 600, unsigned int height = 400);
    CurveRenderer(const CurveRenderer& orig);
    virtual ~CurveRenderer();
    
    void addCurve(const BezierSpline& curve);
    void addHairModel(const Hair& hair, int limitFiberCount = 3);
    //void addCurves(const std::vector<const BezierSpline>& curves);
    
    int startup();
    void init();
    void render();
    
private:
    const unsigned int WINDOW_WIDTH, WINDOW_HEIGHT;
    std::vector<const BezierSpline*> mCurves;


};

#endif /* CURVERENDERER_H */

