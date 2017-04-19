#include "Quark.h"

double Quark::width = 0;

Box2D Quark::getBoundingBox()
{
    Box2D result = boundingBox;
    result.setCenter(parent->getX()+x,parent->getY()+y);
    return result;
}

/* Takes absolute grid coordinates xg, yg.
 * Accounts for its relative position to its parent. */
double Quark::getSmoothTn(double xg, double yg)
{
    double xRel = xg-parent->getX();
    double yRel = yg-parent->getY();
    double d = (xRel-x)*(xRel-x)+(yRel-y)*(yRel-y);
    if(d > 5*width)
        return 0;
    return (1/(2*M_PI*width*width)) * exp(-d/(2*width*width));
}

/* Rotates the particle and all of its children. */
void Quark::rotate(double theta, double phi)
{
    Point3D::rotate(cos(theta),phi);
    boundingBox.setX(x);
    boundingBox.setY(y);
    //Doesn't actually rotate children
}

double Quark::getX()
{
    if(parent)
        return x + parent->getX();
    return x;
}
double Quark::getY()
{
    if(parent)
        return y + parent->getY();
    return y;
}
double Quark::getZ()
{
    if(parent)
        return y + parent->getY();
    return y;
}