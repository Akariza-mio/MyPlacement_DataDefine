#ifndef MYPLACE_H
#define MYPLACE_H

#include "placedata.h"
#include "common.h"

class Bin_2D;
class MyPlacer;
class Bin_2D
{
public:
    POS_2D center;
    POS_2D ll;
    POS_2D ur;
    float width;
    float height;
    float area;

    float nodeDensity;    
    float fillerDensity;   
    float terminalDensity; 
    float DarkDensity;  
    float phi;   
};

class MyPlacer
{
public:
     MyPlacer(PlaceData *_db)
    {
        db = _db;
    }
    vector<vector<Bin_2D *>> bins;

    float StdCellArea; 
    float MacroArea;
        
    PlaceData *db;

    vector<VECTOR_3D> wirelengthGradient; 
    vector<VECTOR_3D> densityGradient;   
    vector<VECTOR_3D> totalGradient;      
    vector<VECTOR_3D> fillerGradient;
};
#endif