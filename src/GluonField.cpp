#include "GluonField.h"
#include <iostream>

double GluonField::distanceScalingFactor = 1;
void GluonField::selectOrigin() {
    double buffer = 10 * Particle::width * GluonField::distanceScalingFactor;
    xOrigin = (std::rand()/(double)RAND_MAX)* (GluonField::gridSize*GluonField::binWidth-2*buffer) + buffer;
    yOrigin = (std::rand()/(double)RAND_MAX)* (GluonField::gridSize*GluonField::binWidth-2*buffer) + buffer;
}
