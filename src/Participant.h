#ifndef Participant_h
#define Participant_h

#include <vector>
#include "Particle.h"

class Participant: public Particle
{
protected:
    double xsave,ysave;
    double fluctfactor;
    double quarkFluctfactors[3];
public:
    std::vector<int> who_hit_me; // those nucleons collided with this nucleon, their indices in the binaryCollision array is stored here. This is not against the capsulate rule since the interpretation of these vector replied on external programs, and this vector is just a storage space for additional info.
    Participant(Particle* part0): Particle(part0) {
        xsave = x;
        ysave = y;
      resetFluctFactors();
    }
    
    double getXorg() {return xsave;}
    double getYorg() {return ysave;}
    void resetCoordinate() {x = xsave; y = ysave;}

    void setFluctfactor(double fluct) {fluctfactor = fluct;}
    void setQuarkFluctfactor(double f1, double f2, double f3);
    void resetFluctFactors();
    double getFluctfactor() {return fluctfactor;}
    
    double getFluctuatedDensity(double xg, double yg);
    double getSmoothDensity(double xg, double yg);
};

#endif
