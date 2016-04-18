/* 
 * File:   Nucleus.h
 * Author: kevin
 *
 * Created on July 1, 2015, 2:35 PM
 */

#ifndef NUCLEUS_H
#define	NUCLEUS_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include "ParameterReader.h"
#include "HulthenFunc.h"
#include "Particle.h"
#include "GaussianDistribution.h"
#include "Box2D.h"

class Nucleus
{
protected:
    vector<Particle*> nucleons;
    vector<Particle*> woundedNucleons;
    GaussianDistribution* quarkDist;
    int nPart;
    
    double lastCx, lastPh;
    
    int deformed;
    double rad,rmaxCut,rwMax;
    double dr;
    double density0;
    double A;  // mass number of nucleus as double
    int atomic;// same as int
    double ecm;

    double beta2,beta4; //deformation parameters 05032010 by TH
    double ctr, phir; // ctr = cos(theta)

    HulthenFunc sample_deuteron;
    vector< vector<double> > triton_pos;

    int flag_NN_correlation;
    int n_configuration;
    double*** nucleon_pos_array;
    bool GFF;

public:
    
    Nucleus(int a, ParameterReader* paraRdr, int deformed=0);
    virtual ~Nucleus();
    double getLastCx1(){return lastCx;}
    double getLastPh1(){return lastPh;}
    int    getAtomic() {return atomic;}
    int    getNpart()  {return woundedNucleons.size();}

    vector<Particle*>& getNucleons() {return nucleons;}
    vector<Particle*>& getParticipants() {return woundedNucleons;}
    void populate(double xCenter, double yCenter);
    void clearNucleons();
    void getRandomWS(double& x, double& y, double& z);
    void markWounded(Particle* part);
    
    static void Gauss38(double xini,double xfin,double* xn,double* wn);
    
    void dumpParticipants(ofstream& of);
    void dumpNucleons(ofstream& of);
    void dumpQuarks(ofstream& of);
    
    void setGluonFields();

    //Deformation
    void getDeformRandomWS(double& x, double& y, double& z);
    void setRotation(double costheta, double phi) {ctr=costheta; phir=phi;}
    double SphericalHarmonics(int l, double theta);
    void readin_nucleon_positions();
    void get_nucleon_position_with_NN_correlation(double **nucleon_ptr);

    // nucleon positions for light nuclei
    void GetDeuteronPosition(double& x1,double& y1,double& z1,double& x2,double& y2,double& z2);
    void readin_triton_position();
    void GetTritonPosition(double& x1,double& y1,double& z1,double &x2,double& y2,double& z2, double &x3, double &y3, double &z3);
    
    // Efficient Collision Calculations
    static bool sortByXLeft(const Particle* p1,const Particle* p2)
    {return p1->getBoundingBox().getXL() < p2->getBoundingBox().getXL();}

};

#endif	/* NUCLEUS_H */

