#ifndef MCnucl_h
#define MCnucl_h
#include <vector>
#include <fstream>
#include "Particle.h"
#include "Nucleus.h"
#include "KLNModel.h"
#include "GlueDensity.h"
#include "CollisionPair.h"
#include "Largex.h"
#include "Spectator.h"

#include "GaussianNucleonsCal.h"
#include "GaussianDistribution.h"
#include "ParameterReader.h"
#include "NBD.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// subroutines to generate and handle MC configurations of colliding nuclei
class MCnucl
{
protected:
    Nucleus* proj;
    Nucleus* targ;
    std::vector<CollisionPair*> binaryCollision;
    std::vector<Spectator*> spectators;

    KLNModel* kln;
    GaussianNucleonsCal *gaussCal; // for Gaussian shaped nucleons calculations
    Large_x* val;
    double **TA1,**TA2;
    double **rho_binary;
    GlueDensity *rho;
    int tmax, tmaxPt;
    double dT;
    double ***dndyTable;
    double ****dndydptTable;
    double dndy;
    int Maxx,Maxy;
    int isKLN;
    double Xmin, Ymin, Xmax, Ymax;
    double PTinte, PTmax, PTmin, dpt, MaxPT;
    int PT_order;
    double dx,dy;
    double siginNN, siginNN200;
    double rapidity;
    int binRapidity;
    double rapMin, rapMax;
    int iptmax;
    int Npart1, Npart2;
    double Xcm, Ycm, angPart;
    int overSample;
    //double Bnucl;// <R^2>/3 of nucleon for Gaussian shape   // static constant double Bnucl = 0.2959
    double dsq; // sigma_in(s) / pi
    int NpartMax, NpartMin;  // accepted range of Npart (centrality cut)
    double BiLinear(int rapBin, double TA, double TB, int ipt);
    double Alpha;

    // For CC fluctuation
    void fluctuateCurrentDensity(int iy);
    int CCFluctuationModel;
    double CCFluctuationK;
    NBD *nbd;
    const gsl_rng_type * gslRngType;
    gsl_rng * gslRng;
    double ccFluctuationGammaTheta;  // scale parameter for gamma distribution

    ParameterReader* paraRdr;
    int shape_of_nucleons;

    int flag_NN_correlation;

    // add by Kelvin Welsh
    int shape_of_entropy;
    int forceCollisionCriterion;
    //quark_width^2 + quark_dist_width^2 = nucleon_width^2
    GaussianDistribution* gaussDist;

    int which_mc_model;
    int sub_model;
    double entropy_gaussian_width, entropy_gaussian_width_sq;

public:
    MCnucl(ParameterReader*);
    ~MCnucl();

    double lastCx1, lastPh1, lastCx2, lastPh2; // store the nuclei rotation angles from the last call
    void setKLN(KLNModel* k) {kln=k;}  // pointer to small-x gluon class
    void setLgX(Large_x* k) {val=k;}   // pointer to large-x partons

    double getTA1(int x,int y) {return TA1[x][y];}
    double getTA2(int x,int y) {return TA2[x][y];}
    double get_rho_binary(int x, int y) {return rho_binary[x][y];}
    double getRho(int i, int x,int y) {return rho->getDensity(i,x,y);}
    double getRho(int i, int x,int y, int pt) {return rho->getDensity(i,x,y,pt);}
    void setRho(int i, int x,int y, double val) {rho->setDensity(i,x,y,val);}
    Box2D getHotSpots(vector<Box2D> & hotSpots);
    int getNcoll() {return binaryCollision.size();}
    int getNpart1() {return Npart1;}
    int getNpart2() {return Npart2;}
    double getdNdy() 
    {
        if(PTinte) return dndy*dx*dy; 
        else return dndy*dx*dy*dpt;
    }

    void setOverSample(int i) {overSample=i;}
    void setRapidity(double y) {rapidity=y;}
    void generateNuclei(double b);
    void deleteNucleus();
    void setDensity(int iy, int ipt); // ipt<0: no dN/dydpt table used
    void addDensity(Nucleus* nucl,double** density);
    void calculateThickness();
    void setThickness(Nucleus* nucl, double ** TA);
    void calculate_rho_binary();    // calculate binary collision density in the transverse plane
    int  getBinaryCollision();
    void createBinaryCollisions();
    void selectFluctFactors(Particle* part);
    int  CentralityCut();
    void setCentralityCut(int Nmin, int Nmax)
                {NpartMax=Nmax; NpartMin=Nmin;}
    void getEccentricityGrid(double& ecc, double& eccp, double& area, double& areap, int iy);
    double getCMAng();
    void recenterGrid(int iy, int n=2);
    void rotateGrid(int iy, int n=2);
    void makeTable();
    void makeTable(double, double, int);
    void dumpdNdyTable4Col(char filename[], double *** dNdyTable, const int iy);
    void dumpdNdydptTable5Col(char filename[], double **** dNdydptTable, const int iy);    
    double getSigEff();
    int hit(Particle* part1, Particle* part2);
    static double Angle(const double x,const double y);

    void dumpparticipantTable(char filename[]);
    void dumpBinaryTable(char filename[]);

    int getSpectators();
    void dumpSpectatorsTable(int event);

    double sampleFluctuationFactorforParticipant();
    double sampleFluctuationFactorforBinaryCollision();
};
#endif
