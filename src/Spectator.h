#ifndef Spectator_h
#define Spectator_h

#include <vector>
#include "Particle.h"
#include "MathBasics.h"

class Spectator: public Point3D {
 private:
    double rapidity_Y;

 public:
    Spectator(double xi, double yi, double rapi) {
        x = xi;
        y = yi;
        rapidity_Y = rapi;
    };
    ~Spectator() {};

    double getX() const {return(x);};
    double getY() const {return(y);};
    double getRapidity_Y() const {return(rapidity_Y);};
    void setX(double x_in) {x = x_in;};
    void setY(double y_in) {y = y_in;};
    void rotate(double theta, double phi) {
        Point3D::rotate(cos(theta),phi);
    }
};

#endif
