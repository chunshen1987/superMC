#ifndef GLUONFIELD_H
#define GLUONFIELD_H

#include "Particle.h"
#include "Box2D.h"
#include <cstdlib>
#include "time.h"

class GluonField
{
	static double ImprintArrayLHC[][600];
	static double ImprintArrayRHIC[][600];
	static double ImprintArrayProton[][600];
	static double ImprintArrayDeuteron[][600];

	int a;
	double energy;
	double xOrigin;
	double yOrigin;
public:
	GluonField(int A,double ecm)
	{
		a = A;
		energy = ecm;
        selectOrigin();
	}
	~GluonField(){};

	void selectOrigin();

	double getFactor(double xg, double yg)
	{
		if(a == 197 || a== 208){
			if(energy == 200)
				return ImprintArrayRHIC[(int)((xg+xOrigin)/0.1)][(int)((yg+yOrigin)/0.1)];
			if(energy == 5020)
				return ImprintArrayLHC[(int)((xg+xOrigin)/0.1)][(int)((yg+yOrigin)/0.1)];
		}
		else if(a == 1)
			return ImprintArrayProton[(int)((xg+xOrigin)/0.1)][(int)((yg+yOrigin)/0.1)];
		else if(a == 2)
			return ImprintArrayDeuteron[(int)((xg+xOrigin)/0.1)][(int)((yg+yOrigin)/0.1)];
	}
};


#endif