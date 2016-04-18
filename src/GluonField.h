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
	const static double binWidth = 0.1;
	const static double distanceScalingFactor = 0.5;
public:
	GluonField()
	{
		a = -1;
		energy = -1;
	}
	GluonField(int A,double ecm)
	{
		a = A;
		energy = ecm;
        selectOrigin();
	}
	~GluonField()
	{
	};

	void selectOrigin();

	double getFactor(double xg, double yg)
	{
		int x = (int)((xg+xOrigin)*distanceScalingFactor/binWidth);
		int y = (int)((yg+yOrigin)*distanceScalingFactor/binWidth);
		if(a == 197 || a== 208){
			if(energy == 200)
				return ImprintArrayRHIC[x][y];
			if(energy == 5020)
				return ImprintArrayLHC[x][y];
		}
		else if(a == 1)
			return ImprintArrayProton[x][y];
		else if(a == 2)
			return ImprintArrayDeuteron[x][y];
		else
			// This is not a correct treatment, but mearly a temporary one.
			return ImprintArrayProton[x][y];
	}
};


#endif