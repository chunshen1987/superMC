/* 
 * File:   MathVector.h
 * Author: kevin
 *
 * Created on June 30, 2015, 1:59 PM
 */

#ifndef MATHVECTOR_H
#define	MATHVECTOR_H

struct MathVector
{
public:
    double x;
    double y;
    double z;
    
    MathVector(double xin = 0, double yin = 0, double zin = 0)
    {
        x = xin;
        y = yin;
        z = zin;
    };
    
    double dot(MathVector other)
    {
        return x*other.x+y*other.y+z*other.z;
    }
    
    MathVector cross(MathVector other)
    {
        MathVector result;
        result.x = y*other.z-z*other.y;
        result.y = z*other.x-x*other.z;
        result.z = x*other.y-y*other.x;
        
        return result;
    }
};

inline MathVector operator+(const MathVector lhs, const MathVector rhs)
{
    MathVector result(lhs.x+rhs.x,lhs.y+rhs.y,lhs.z+rhs.z);
    return result;
}

inline MathVector operator-(const MathVector lhs, const MathVector rhs)
{
    MathVector result(lhs.x-rhs.x,lhs.y-rhs.y,lhs.z-rhs.z);
    return result;
}

#endif	/* MATHVECTOR_H */

