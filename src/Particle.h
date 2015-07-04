#ifndef PARTICLE_h
#define PARTICLE_h

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "GaussianDistribution.h"
#include "MathBasics.h"

class Particle: public Point3D
{
 protected:
    int numberOfCollision;

    // added by Kevin Welsh
    GaussianDistribution* quarkDist;
    double quarkWidth;
    double nucleonWidth;
    double xsave,ysave,zsave;
    double fluctfactor;
    double quarkFluctfactors[3];
    std::vector<Particle*> who_hit_me;

 public:
    double ValenceQuarks[3][3];
    Particle(double x0,double y0, double z0, double nWidth = 0, double qWidth = 0);
    ~Particle();

    double getXorg() {return xsave;}
    double getYorg() {return ysave;}
    void   saveOrg() {xsave = x; ysave = y; zsave = z;}
    void resetCoordinate() {x = xsave; y = ysave; z = zsave;}

    void setFluctfactor(double fluct) {fluctfactor = fluct;}
    void setQuarkFluctfactor(double f1, double f2, double f3);
    void resetFluctFactors();
    
    std::vector<Particle*>& getCollidingParticles(){return who_hit_me;}
    double getFluctfactor() {return fluctfactor;}
    double getWidth() {return nucleonWidth;}
    double getQuarkWidth() {return quarkWidth;}
    void   addCollidingParticle(Particle* colliding){who_hit_me.push_back(colliding);}
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

    int    getNumberOfCollision() {return numberOfCollision;}
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
