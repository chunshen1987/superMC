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

Particle::Particle(Particle* inPart)
{
    x = inPart->getX();
    y = inPart->getY();
    z = inPart->getZ();
    numberOfCollision=0;
    nucleonWidth = inPart->getWidth();
    quarkWidth = inPart->getQuarkWidth();
    nucl = inPart->getNucl();
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            ValenceQuarks[i][j] = (inPart->ValenceQuarks[i][j]);
}

Particle::~Particle()
{
}

void Particle::generateQuarkPositions()
{
   double x, y;

   int i=0;
   while(i<3)
   {
     x = gaussDist->rand();
     y = gaussDist->rand();
     ValenceQuarks[i][0] = x;
     ValenceQuarks[i][1] = y;
     i++;
   }
   cout << endl;
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