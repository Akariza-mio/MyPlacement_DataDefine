#pragma once
#include "objects.h"
#include <map>
#include <vector>
#include <string>

class PlaceData{
public:
    int moduleCount; 
    int MacroCount;
    int netCount;
    int pinCount;

    vector<Module*> Nodes; 
    vector<Module*> Terminals;
    vector<Pin*> Pins;
    vector<Net*> Nets;
    vector<SiteRow> SiteRows;

    map<string, Module*> moduleMap; 
    //网络名称索引表
    map<string, Net*> netMap;

    void readBookShelf(const string& dir,const string& benchmarkname);
    void readNodes(const string& filename);
    void readPl(const string& filename);
    void readScl(const string& filename);
    void readNets(const string& filename);
    void readWts(const string& filename);
    bool readAux(
        string& nodes_filepath,
        string& pl_filepath,
        string& scl_filepath,
        string& nets_filepath,
        string& wts_filepath,
        const string& aux_filepath
    );
};
