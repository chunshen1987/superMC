/* 
 * File:   Box2D.h
 * Author: kevin
 *
 * Created on July 3, 2015, 3:35 PM
 */

#ifndef BOX2D_H
#define	BOX2D_H

#include "arsenal.h"
#include <sstream>
#include <iostream>
#include <string>

class Box2D
{
protected:
    double xLeft, xRight, yLeft, yRight;
    double xCenter, yCenter;
public:
    Box2D(double xL, double xR, double yL, double yR);
    Box2D();
    double getXR() const {return xRight;}
    double getXL() const {return xLeft;}
    double getYR() const {return yRight;}
    double getYL() const {return yLeft;}
    double getX() const  {return xCenter;}
    double getY() const  {return yCenter;}
    
    void setCenter(double x, double y);
    void setX(double x){setCenter(x,yCenter);}
    void setY(double y){setCenter(xCenter,y);}
    void setDimensions(double xLength, double yLength);
    void setSquareDimensions(double size);
    
    Box2D intersection(const Box2D& other);
    
    void printToConsole();
    
    std::string toString()
    {
        std::ostringstream os;
        os << xLeft << " " << xRight << " " << yLeft << " " << yRight;
        //os << "(" << xLeft << "," << xRight << "," << yLeft << "," << yRight << ")";
        return os.str();
    }

    void overUnion(const Box2D& other)
    {
        xLeft = min(other.getXL(),xLeft);
        xRight = max(other.getXR(),xRight);
        yLeft = min(other.getYL(),yLeft);
        yRight = max(other.getYR(),yRight);
        xCenter = (xLeft+xRight)/2.0;
        yCenter = (yLeft+yRight)/2.0;
        
        if(xLeft > other.getXL()
                || xRight < other.getXR()
                || yLeft > other.getYL()
                || yRight < other.getYR())
            cout << toString() << endl;;
    }
    
    bool isPositive()
    {
        return !(xLeft >= xRight || yLeft >= yRight);
    }
    
    double getArea()
    {
        return (xRight-xLeft)*(yRight-yLeft);
    }
};

#endif	/* BOX2D_H */

