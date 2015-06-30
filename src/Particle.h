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
  GaussianDistribution* gaussDist;
  double quarkWidth;
  double nucleonWidth;

 public:
  double ValenceQuarks[3][3];
  Particle(double x0,double y0, double z0, double nWidth = 0, double qWidth = 0);
  Particle(Particle* inPart);
  ~Particle();

  double getX() {return x;}
  double getY() {return y;}
  double getZ() {return z;}
  double getWidth() {return nucleonWidth;}
  double getQuarkWidth() {return quarkWidth;}
  int    getNucl() {return nucl;}
  void   setX(double a) {x=a;}
  void   setY(double a) {y=a;}
  void   setZ(double a) {z=a;}
  void   setNucl(int a) {nucl=a;}
  void   setWidth(double nW,double qW)
  {
      nucleonWidth = nW;
      quarkWidth = qW;
      double quark_dist_width = sqrt(nW*nW-qW*qW);
  }
  void   setGaussDist(GaussianDistribution* dist)
  {
      gaussDist = dist;
      generateQuarkPositions();
  }
  
  void   setQuarkWidth(double a) {quarkWidth=a;}

  int    getNumberOfCollision() {return numberOfCollision;}
  void   setNumberOfCollision() {numberOfCollision += 1;}
  void   setNumberOfCollision(int i) {numberOfCollision=i;}

  // functions for nucleon substructure added by Kevin Welsh
  void generateQuarkPositions();
  
  double getFluctuatedTn(double xg, double yg);
  double getSmoothTn(double xg, double yg);

};

#endif
