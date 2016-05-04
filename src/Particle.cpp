#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Particle.h"
#include "Quark.h"
#include "GaussianDistribution.h"
using namespace std;


double Particle::width = 0;
double Particle::R = 0;
vector< vector<double> > Particle::quark_pos;

Particle::Particle(double x0, double y0, double z0): Point3D(x0,y0,z0)
{
   gluonField = NULL;
   baseBox.setCenter(x0,y0);
   baseBox.setSquareDimensions(8*width);
   boundingBox = baseBox;
   generateQuarkPositions();
   resetFluctFactors();
}

Particle::~Particle()
{
  if(gluonField)
    delete gluonField;
}

void Particle::generateQuarkPositions()
{
  ValenceQuarks.clear();
  // Note* Quarks are generated with a position relative to 
  // their parent nucleon.
  double r1, r2, z12;
  int index = (int)250000*drand48();
  r1 = quark_pos[index][0];
  r2 = quark_pos[index][1];
  z12 = quark_pos[index][2];
  double r1x, r1y, r1z, r2x, r2y, r2z;
  r1=r1*R;
  r2=r2*R;

  double Theta1,phi1,phi2; // polar, azimuthal, and "conical" angles
  double Theta12 = acos(z12);
  double ux,uy,uz,vx,vy,vz,c,s;
  // randomly generate polar angle theta wrt arbitrary z axis:
  double z1 = 2.*drand48()-1.;
  Theta1=acos(z1);
  // randomly generate azimuthal angle:
  phi1 = 2*M_PI*drand48();
  phi2 = 2*M_PI*drand48();
  ux=sin(Theta1)*cos(phi1);
  uy=sin(Theta1)*sin(phi1);
  uz=z1;
  vx=sin(Theta1+Theta12)*cos(phi1);
  vy=sin(Theta1+Theta12)*sin(phi1);
  vz=cos(Theta1+Theta12);
  c=cos(phi2);
  s=sin(phi2);
  
  r1x = r1*ux;
  r1y = r1*uy;
  r1z = r1*uz;

  r2x = vx*(c+ux*ux*(1-c))+vy*(ux*uy*(1-c)-uz*s)+vz*(ux*uz*(1-c)+uy*s);
  r2y = vx*(ux*uy*(1-c)+uz*s)+vy*(c+uy*uy*(1-c))+vz*(uy*uz*(1-c)-ux*s);
  r2z = vx*(uz*ux*(1-c)-uy*s)+vy*(uz*uy*(1-c)+ux*s)+vz*(c+uz*uz*(1-c));

  r2x = r2x*r2;
  r2y = r2y*r2;
  r2z = r2z*r2;

  ValenceQuarks.push_back(Quark(this,r1x,r1y,r1z));
  ValenceQuarks.push_back(Quark(this,r2x,r2y,r2z));
  ValenceQuarks.push_back(Quark(this,-r1x-r2x,-r1y-r2y,-r1z-r2z));

  /* double cmOfQuarks[] = {0,0};
   for(int i = 0; i < 2; i++){
       ValenceQuarks.push_back(Quark(this,quarkDist->rand(),quarkDist->rand()));
       cmOfQuarks[0] += ValenceQuarks[i].getX();
       cmOfQuarks[1] += ValenceQuarks[i].getY();
       boundingBox.overUnion(ValenceQuarks[i].getBoundingBox());
   }

   ValenceQuarks.push_back(Quark(this,-cmOfQuarks[0],-cmOfQuarks[1]));
   boundingBox.overUnion(ValenceQuarks[2].getBoundingBox());*/
}

void Particle::calculateBounds()
{
    boundingBox = baseBox;
    for(int i = 0; i < 3; i++)
        boundingBox.overUnion(ValenceQuarks[i].getBoundingBox());   
}


/* Returns the fluctuated nuclear thickness at a given 
 * point in space. Ignores multiplicity fluctuations. */
double Particle::getFluctuatedTn(double xg, double yg)
{
   double dens = 0;
   
   for(int i(0);i<3;i++)
   {
        //Divide a total of 1 density between 3 quarks
    double s = ValenceQuarks[i].getSmoothTn(xg,yg)/3.0;
    if(gluonField)
      s*= gluonField->getFactor(xg-x,yg-y);
   	dens += s;
   }
   return dens;
}

/* Returns the nuclear thickness at a given 
 * point in space based on a Gaussian representation.
 * Ignores multiplicity fluctuations. */
double Particle::getSmoothTn(double xg, double yg)
{
    double r = sqrt((xg-x)*(xg-x)+(yg-y)*(yg-y));
    if(r > 5*width)
      return 0;
    double Tn = (1/(2*M_PI*width*width))*exp(-r*r/(2*width*width));
    if(gluonField)
      return Tn*gluonField->getFactor(xg-x,yg-y);
    else
      return Tn;
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
    setQuarkFluctfactor(1.0/3,1.0/3,1.0/3);
}

/* Returns the gluon density at a 
 * given point in space. Considers sub-nucleonic structure.*/
double Particle::getFluctuatedDensity(double xg, double yg)
{
   double dens = 0;
   for(int i(0);i<3;i++)
   {
    //Divide a total of 1 density between 3 quarks
    // The factor of 1/3 is in the multiplicty fluctuation factor
    double s = ValenceQuarks[i].getSmoothDensity(xg,yg);
    if(gluonField)
      s*= gluonField->getFactor(xg-x,yg-y);

   	dens += s; 
   }
   return dens;
}

/* Returns the gluon density at a give point,
 * but taken from a Gaussian distribution.
 * Includes multiplicity fluctuations*/
double Particle::getSmoothDensity(double xg, double yg)
{
  if(gluonField)
    return getSmoothTn(xg,yg)*fluctfactor*gluonField->getFactor(xg-x,yg-y);
  else
    return getSmoothTn(xg,yg)*fluctfactor;

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
    boundingBox.setX(x);
    boundingBox.setY(y);
}

void rotateQuarks(double r1, double r2, double z12, 
        double &r1x, double &r1y, double &r1z,
        double &r2x, double &r2y, double &r2z)
{
  double Theta1,phi1,phi2; // polar, azimuthal, and "conical" angles
  double Theta12 = acos(z12);
  double ux,uy,uz,vx,vy,vz,c,s;
  // randomly generate polar angle theta wrt arbitrary z axis:
  double z1 = 2.*drand48()-1.;
  Theta1=acos(z1);
  // randomly generate azimuthal angle:
  phi1 = 2*M_PI*drand48();
  phi2 = 2*M_PI*drand48();
  ux=sin(Theta1)*cos(phi1);
  uy=sin(Theta1)*sin(phi1);
  uz=cos(Theta1);
  vx=sin(Theta1+Theta12)*cos(phi1);
  vy=sin(Theta1+Theta12)*sin(phi1);
  vz=cos(Theta1+Theta12);
  c=cos(phi2);
  s=sin(phi2);
  
  r1x = r1*ux;
  r1y = r1*uy;
  r1z = r1*uz;
  r2x = vx*(c+ux*ux*(1-c))+vy*(ux*uy*(1-c)-uz*s)+vz*(ux*uz*(1-c)+uy*s);
  r2y = vx*(ux*uy*(1-c)+uz*s)+vy*(c+uy*uy*(1-c))+vz*(uy*uz*(1-c)-ux*s);
  r2z = vx*(uz*ux*(1-c)-uy*s)+vy*(uz*uy*(1-c)+ux*s)+vz*(c+uz*uz*(1-c));
  r2x = r2x*r2;
  r2y = r2y*r2;
  r2z = r2z*r2;
}