#ifndef PARTICLE_h
#define PARTICLE_h

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "GaussianDistribution.h"
#include "MathBasics.h"
#include "Box2D.h"
#include "IGluonSource.h"
#include "Quark.h"

class Quark;

class Particle: private Point3D, public IGluonSource
{
 protected:
    int numberOfCollision;
    
    double xsave,ysave,zsave;
    double fluctfactor;
    double quarkFluctfactors[3];
    std::vector<Quark> ValenceQuarks;
    std::vector<Particle*> who_hit_me;
    Box2D boundingBox;

 public:
    static double width;
    static GaussianDistribution* quarkDist;
    Particle(double x0=0,double y0=0, double z0=0);
    ~Particle();

    double getXorg() {return xsave;}
    double getYorg() {return ysave;}
    void   saveOrg() {xsave = x; ysave = y; zsave = z;}
    void resetCoordinate() {x = xsave; y = ysave; z = zsave;}
    
    vector<Quark>& getQuarks(){return ValenceQuarks;}
    double getX(){return x;}
    double getY(){return y;}
    double getZ(){return z;}
    void   setX(double a);
    void   setY(double a);
    void   setZ(double a);
    
    Box2D  getBoundingBox() const {return boundingBox;}
    void   rotate(double theta, double phi);

    void setFluctfactor(double fluct) {fluctfactor = fluct;}
    void setQuarkFluctfactor(double f1, double f2, double f3);
    void resetFluctFactors();
    
    std::vector<Particle*>& getCollidingParticles(){return who_hit_me;}
    double getFluctfactor() {return fluctfactor;}
    double getWidth() {return width;}
    void   addCollidingParticle(Particle* colliding){who_hit_me.push_back(colliding);}
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
