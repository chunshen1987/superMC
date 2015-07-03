#include <cmath>
#include <fstream>
#include <istream>
#include <iomanip>

#include "Particle.h"
#include "GaussianDistribution.h"
using namespace std;

Particle::Particle(double x0, double y0, double z0, double nWidth, double qWidth)
{
   x = x0; 
   y = y0; 
   z = z0;
   numberOfCollision=0;
   nucleonWidth = nWidth;
   quarkWidth = qWidth;
}

Particle::~Particle()
{
}

void Particle::generateQuarkPositions()
{
   for(int i = 0; i < 3; i++)
   {
     ValenceQuarks[i][1] = quarkDist->rand();
     ValenceQuarks[i][0] = quarkDist->rand();
   }
}


/* Returns the fluctuated nuclear thickness at a given 
 * point in space. Ignores multiplicity fluctuations. */
double Particle::getFluctuatedTn(double xg, double yg)
{
   double dens = 0;
   
   double d;
   for(int i(0);i<3;i++)
   {
        //The squared distance between the Quark and the grid
   	d = pow(ValenceQuarks[i][0]+x-xg,2)+pow(ValenceQuarks[i][1]+y-yg,2); 
        //Divide a total of 1 density between 3 quarks
   	dens += (1./(2*M_PI*quarkWidth*quarkWidth))*exp(-d/(2*quarkWidth*quarkWidth))/3; 
   }
   return dens;
}

/* Returns the nuclear thickness at a given 
 * point in space based on a Gaussian representation.
 * Ignores multiplicity fluctuations. */
double Particle::getSmoothTn(double xg, double yg)
{
    double r = sqrt((xg-x)*(xg-x)+(yg-y)*(yg-y));
    return (1/(2*nucleonWidth*nucleonWidth))*
            exp(-r*r/(2*nucleonWidth*nucleonWidth));
}

void Particle::setQuarkFluctfactor(double f1, double f2, double f3)
{
    quarkFluctfactors[0]=f1;
    quarkFluctfactors[1]=f2;
    quarkFluctfactors[2]=f3;
}

void Particle::resetFluctFactors()
{
    fluctfactor = 1;
    quarkFluctfactors[0] = 1;
    quarkFluctfactors[1] = 1;
    quarkFluctfactors[2] = 1;
}

/* Returns the gluon density at a 
 * given point in space. Considers sub-nucleonic structure.*/
double Particle::getFluctuatedDensity(double xg, double yg)
{
   if(numberOfCollision == 0)
       return 0;
   double dens = 0;
   double d;
   for(int i(0);i<3;i++)
   {
        //The squared distance between the Quark and the grid
   	d = pow(ValenceQuarks[i][0]+x-xg,2)+pow(ValenceQuarks[i][1]+y-yg,2); 
        //Divide a total of 1 density between 3 quarks
   	dens += quarkFluctfactors[i]*(1./(2*M_PI*quarkWidth*quarkWidth))*exp(-d/(2*quarkWidth*quarkWidth))/3; 
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
    return fluctfactor*(1/(2*nucleonWidth*nucleonWidth))*
            exp(-r*r/(2*nucleonWidth*nucleonWidth));
}
