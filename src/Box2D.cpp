#include "Box2D.h"
#include "iostream"

Box2D::Box2D(double xL, double xR, double yL, double yR)
{
    xLeft = xL;
    xRight = xR;
    yLeft = yL;
    yRight = yR;
    xCenter = (xR+xL)/2.0;
    yCenter = (yR+yL)/2.0;
}

Box2D::Box2D()
{
    xLeft = 0;
    xRight = 0;
    yLeft = 0;
    yRight = 0;
    xCenter = 0;
    yCenter = 0;
}

void Box2D::setCenter(double x, double y)
{
    xLeft += x-xCenter;
    xRight += x-xCenter;
    yLeft += y-yCenter;
    yRight += y-yCenter;
    
    xCenter = x;
    yCenter = y;
}

void Box2D::setDimensions(double xLength, double yLength)
{
    xLeft = xCenter - xLength/2;
    xRight = xCenter + xLength/2;
    yLeft = yCenter - yLength/2;
    yRight = yCenter + yLength/2;
}

void Box2D::setSquareDimensions(double size)
{
    setDimensions(size,size);
}

/* Gets box that describes the intersection between 2
 * boxes. Does not guard against the non-intersecting case.
 * Use Box2D::isPositive() to check if a box makes sense. */
Box2D Box2D::intersection(const Box2D& other)
{
    double xL, xR, yL, yR;
    xL = max(xLeft,other.xLeft);
    xR = min(xRight,other.xRight);
    yL = max(yLeft,other.yLeft);
    yR = min(yRight,other.yRight);

    return Box2D(xL,xR,yL,yR);
}

void Box2D::printToConsole()
{
    cout << "(" << xLeft << "," << xRight
         << "," << yLeft << "," << yRight << ")" << endl;
}
