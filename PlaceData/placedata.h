#pragma once
#ifndef PLACEDATA_H
#define PLACEDATA_H
#include "Common/objects.h"
#include <map>
#include <vector>


class PlaceData{
public:
    int moduleCount; 
    int MacroCount;
    int netCount;
    int pinCount;

    vector<Module *> Nodes; 
    vector<Module *> Terminals;
    vector<Pin *> Pins;
    vector<Net *> Nets;
    vector<SiteRow> SiteRows;

    map<string, Module *> moduleMap; 

    void readBookShelf(const string& dir,const string& benchmarkname);
    void readNodes(const string& filename);
    void readPl(const string& filename);
    void readScl(const string& filename);
    void readNets(const string& filename);
};
#endif