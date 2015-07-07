/* 
 * File:   IGluonSource.h
 * Author: kevin
 *
 * Created on July 4, 2015, 5:22 PM
 */

#ifndef IGLUONSOURCE_H
#define	IGLUONSOURCE_H

class IGluonSource
{
public:
    virtual double getSmoothDensity(double x, double y){return 0;}
    virtual double getFluctuatesDensity(double x, double y){return 0;}
    
    virtual double getX() = 0;
    virtual double getY() = 0;
    virtual double getZ() = 0;
    virtual void   setX(double a) = 0;
    virtual void   setY(double a) = 0;
    virtual void   setZ(double a) = 0;
};

#endif	/* IGLUONSOURCE_H */

