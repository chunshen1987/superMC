#ifndef CollisionPair_h
#define CollisionPair_h

#include "MathBasics.h"
#include "Box2D.h"
#include "IGluonSource.h"

class CollisionPair : private Point3D, public IGluonSource
{
 protected:
  double fluctfactor;
  Box2D boundingBox;
 public:
  double additional_weight; // store additional weight used in Uli-Glb model
  CollisionPair(double x0,double y0):Point3D(x0,y0){
    additional_weight = 0.;
  }
  ~CollisionPair() {};
    
  Box2D getBoundingBox(){return boundingBox;}
    double getX(){return x;}
    double getY(){return y;}
    double getZ(){return z;}
    void setX(double a)
    {
        x = a;
        boundingBox.setX(a);
    }
    void setY(double a)
    {
        y = a;
        boundingBox.setY(a);
    }
    void setZ(double a)
    {
        z = a;
    }
    
    void rotate(double theta, double phi)
    {
        Point3D::rotate(cos(theta),phi);
    }

  void setfluctfactor(double fluct) {fluctfactor = fluct;}
  double getfluctfactor() const {return fluctfactor;}
};
#endif
