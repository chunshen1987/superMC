/* 
 * File:   Quark.h
 * Author: kevin
 *
 * Created on July 4, 2015, 8:38 PM
 */

#ifndef QUARK_H
#define	QUARK_H

#include <cmath>
#include "Particle.h"
#include "Box2D.h"
#include "MathBasics.h"
#include "IGluonSource.h"

class Particle;

class Quark: private Point3D, public IGluonSource
{
    Box2D boundingBox;
    double fluctFactor;
    Particle* parent;
    
    public:
        static double width;
        
        Quark(Particle* inParent, double x0 = 0, double y0 = 0, double z0 = 0)
                : Point3D(x0,y0,z0)
        {
            parent = inParent;
            boundingBox.setCenter(x0,y0);
            boundingBox.setSquareDimensions(8*width);
        }

        double getX();
        double getY();
        double getZ();

        double getLocalX(){return x;}
        double getLocalY(){return y;}
        double getLocalZ(){return z;}
        void   setX(double a)
        {x=a; boundingBox.setX(a);}
        void   setY(double a)
        {y=a; boundingBox.setY(a);}
        void   setZ(double a){z=a;}
        void   setFluctFactor(double a){fluctFactor=a;}
        void   rotate(double theta, double phi);
        
        Box2D  getBoundingBox();
        double getSmoothTn(double xg, double yg);
        double getSmoothDensity(double xg, double yg)
        {
            
            return fluctFactor*getSmoothTn(xg,yg);
        }
        
};


#endif	/* QUARK_H */

