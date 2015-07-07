#include "GluonField.h"
#include <iostream>
void GluonField::selectOrigin()
{
	double buffer = 10 * Particle::width;
	std::srand(time(NULL));
	xOrigin = (std::rand()/(double)RAND_MAX)* (60-2*buffer) + buffer;
	yOrigin = (std::rand()/(double)RAND_MAX)* (60-2*buffer) + buffer;
}
