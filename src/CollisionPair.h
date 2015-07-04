#ifndef CollisionPair_h
#define CollisionPair_h

#include "MathBasics.h"

class CollisionPair : public Point3D
{
 protected:
  double fluctfactor;
 public:
  double additional_weight; // store additional weight used in Uli-Glb model
  CollisionPair(double x0,double y0):Point3D(x0,y0){
    additional_weight = 0.;
  }
  ~CollisionPair() {};

  void setfluctfactor(double fluct) {fluctfactor = fluct;}
  double getfluctfactor() const {return fluctfactor;}
};
#endif
