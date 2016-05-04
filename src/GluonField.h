#ifndef GLUONFIELD_H
#define GLUONFIELD_H

#include "Particle.h"
#include "Box2D.h"
#include <cstdlib>
#include "time.h"

class GluonField
{
	static const int gridSize = 600;
	const static double binWidth = 0.1;
	static double ImprintArrayLHC[][gridSize];
	static double ImprintArrayRHIC[][gridSize];
	//static double ImprintArrayProton[][600];
	//static double ImprintArrayDeuteron[][600];

	int a;
	double energy;
	double xOrigin;
	double yOrigin;
public:
	static double distanceScalingFactor;
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
		int x = (int)((xg*distanceScalingFactor+xOrigin)/binWidth);
		int y = (int)((yg*distanceScalingFactor+yOrigin)/binWidth);
		if(energy == 200)
			return ImprintArrayRHIC[x][y];
		if(energy == 5020)
			return ImprintArrayLHC[x][y];
		/*else if(a == 1)
			return ImprintArrayProton[x][y];
		else if(a == 2)
			return ImprintArrayDeuteron[x][y];
		else
			// This is not a correct treatment, but mearly a temporary one.
			return ImprintArrayProton[x][y];*/
	}
};


#endif