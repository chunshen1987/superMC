#ifndef PARTICLE_h
#define PARTICLE_h

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "GaussianDistribution.h"

class Particle
{
 protected:
    
    int nucl;
    double x,y,z;
    int numberOfCollision;

    // added by Kevin Welsh
    GaussianDistribution* quarkDist;
    double quarkWidth;
    double nucleonWidth;
    double xsave,ysave;
    double fluctfactor;
    double quarkFluctfactors[3];

 public:
    std::vector<int> who_hit_me;
    double ValenceQuarks[3][3];
    Particle(double x0,double y0, double z0, double nWidth = 0, double qWidth = 0);
    ~Particle();

    double getXorg() {return xsave;}
    double getYorg() {return ysave;}
    void resetCoordinate() {x = xsave; y = ysave;}

    void setFluctfactor(double fluct) {fluctfactor = fluct;}
    void setQuarkFluctfactor(double f1, double f2, double f3);
    void resetFluctFactors();
    
    double getFluctfactor() {return fluctfactor;}
    double getX() {return x;}
    double getY() {return y;}
    double getZ() {return z;}
    double getWidth() {return nucleonWidth;}
    double getQuarkWidth() {return quarkWidth;}
    void   setX(double a) {x=a;}
    void   setY(double a) {y=a;}
    void   setZ(double a) {z=a;}
    void   setWidth(double nW,double qW)
    {
        nucleonWidth = nW;
        quarkWidth = qW;
    }
    void   setQuarkDist(GaussianDistribution* &dist)
    {
        quarkDist = dist;
        generateQuarkPositions();
    }

    int    getNucl(){return nucl;}
    int    getNumberOfCollision() {return numberOfCollision;}
    void   setNucl(int a){nucl=a;}
    void   setNumberOfCollision() {numberOfCollision += 1;}
    void   setNumberOfCollision(int i) {numberOfCollision=i;}

    // functions for nucleon substructure added by Kevin Welsh
    void generateQuarkPositions();

    double getFluctuatedTn(double xg, double yg);
    double getSmoothTn(double xg, double yg);
    double getFluctuatedDensity(double xg, double yg);
    double getSmoothDensity(double xg, double yg);
};

#endif
