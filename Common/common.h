#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <list>
#include <set>
#include "string.h"

using namespace std;

struct POS_2D // POS means postition which can be used to store coordinates, offsets
{
    float x;
    float y;
    POS_2D() { SetZero(); };
    POS_2D(float _x, float _y)
    {
        x = _x;
        y = _y;
    }
    inline void SetZero()
    {
        x = y = 0.0; 
    }
    friend inline std::ostream &operator<<(std::ostream &os, const POS_2D &pos)
    {
        os << "(" << pos.x << "," << pos.y << ")";
        return os;
    }
};

struct VECTOR_3D
{
    float x;
    float y;
    float z;
    VECTOR_3D()
    {
        SetZero();
    }
    inline void SetZero()
    {
        x = y = z = 0.0; //!! 0.0!!!!
    }
    friend inline std::ostream &operator<<(std::ostream &os, const VECTOR_3D &vec)
    {
        os << "[" << vec.x << "," << vec.y << "," << vec.z << "]"; // [] for vectors and () for pos
        return os;
    }
    inline VECTOR_3D operator+(const VECTOR_3D &rhs)
    {
        VECTOR_3D v;
        v.x = this->x + rhs.x;
        v.y = this->y + rhs.y;
        v.z = this->z + rhs.z;
        return v;
    }
    inline VECTOR_3D operator-(const VECTOR_3D &rhs) const
    {
        VECTOR_3D v;
        v.x = this->x - rhs.x;
        v.y = this->y - rhs.y;
        v.z = this->z - rhs.z;
        return v;
    }
    inline VECTOR_3D operator*(float c)
    {
        VECTOR_3D v;
        v.x = this->x * c;
        v.y = this->y * c;
        v.z = this->z * c;
        return v;
    }
    inline float operator*(const VECTOR_3D &rhs) const
    {
        return (this->x * rhs.x + this->y * rhs.y + this->z * rhs.z);
    }
};
#endif