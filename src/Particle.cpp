#include <cmath>
#include <fstream>
#include <istream>
#include <iomanip>

#include "Particle.h"
#include "Quark.h"
#include "GaussianDistribution.h"
using namespace std;


double Particle::width = 0;
GaussianDistribution* Particle::quarkDist;

Particle::Particle(double x0, double y0, double z0): Point3D(x0,y0,z0)
{
   numberOfCollision=0;
   
   boundingBox.setCenter(x0,y0);
   boundingBox.setSquareDimensions(5*width);
   
   generateQuarkPositions();
}

Particle::~Particle()
{
}

void Particle::generateQuarkPositions()
{
    // Note* Quarks are generated with a position relative to 
    // their parent nucleon.
   for(int i = 0; i < 3; i++){
       ValenceQuarks.push_back(Quark(this,quarkDist->rand(),quarkDist->rand()));
       boundingBox.overUnion(ValenceQuarks[i].getBoundingBox());
   }
       
}


/* Returns the fluctuated nuclear thickness at a given 
 * point in space. Ignores multiplicity fluctuations. */
double Particle::getFluctuatedTn(double xg, double yg)
{
   double dens = 0;
   
   for(int i(0);i<3;i++)
   {
        //Divide a total of 1 density between 3 quarks
   	dens += ValenceQuarks[i].getSmoothTn(xg,yg)/3; 
   }
   return dens;
}

/* Returns the nuclear thickness at a given 
 * point in space based on a Gaussian representation.
 * Ignores multiplicity fluctuations. */
double Particle::getSmoothTn(double xg, double yg)
{
    double r = sqrt((xg-x)*(xg-x)+(yg-y)*(yg-y));
    return (1/(2*width*width))*
            exp(-r*r/(2*width*width));
}

void Particle::setQuarkFluctfactor(double f1, double f2, double f3)
{
    ValenceQuarks[0].setFluctFactor(f1);
    ValenceQuarks[1].setFluctFactor(f2);
    ValenceQuarks[2].setFluctFactor(f3);
}

void Particle::resetFluctFactors()
{
    fluctfactor = 1;
    setQuarkFluctfactor(1,1,1);
}

/* Returns the gluon density at a 
 * given point in space. Considers sub-nucleonic structure.*/
double Particle::getFluctuatedDensity(double xg, double yg)
{
   if(numberOfCollision == 0)
       return 0;
   double dens = 0;
   for(int i(0);i<3;i++)
   {
        //Divide a total of 1 density between 3 quarks
   	dens += ValenceQuarks[i].getSmoothDensity(xg,yg)/3; 
   }
   return dens;
}

/* Returns the gluon density at a give point,
 * but taken from a Gaussian distribution.
 * Includes multiplicity fluctuations*/
double Particle::getSmoothDensity(double xg, double yg)
{
    if(numberOfCollision == 0)
       return 0;
    double r = sqrt((xg-x)*(xg-x)+(yg-y)*(yg-y));
    return fluctfactor*(1/(2*width*width))*
            exp(-r*r/(2*width*width));
}

void Particle::setX(double a)
{
    x = a;
    boundingBox.setX(a);
}
void Particle::setY(double a)
{
    y = a;
    boundingBox.setY(a);
}
void Particle::setZ(double a)
{
    z = a;
}

void Particle::rotate(double theta, double phi)
{
    Point3D::rotate(cos(theta),phi);
    for(int i = 0; i < ValenceQuarks.size(); i++)
    {
        ValenceQuarks[i].rotate(theta,phi);
    }
}