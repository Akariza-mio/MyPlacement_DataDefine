#include "placedata.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

void PlaceData::readNodes(const string& filename){
    //读入
    ifstream ifs(filename);
    if(!ifs.is_open()){
        cerr << "failed to open file: " << filename << endl;
        return;
    }
    string line;
    //行读取
    while(getline(ifs,line)){
        //跳过无用行
        if(line.empty() || line[0] == '#' || line[0] == 'U'){
            continue;
        }
        //切割字符流
        stringstream ss(line);
        string token;
        ss >> token;
        if(token == "NumNodes"){
            string colon;
            ss >> colon >> moduleCount;
        }else if(token == "NumTerminals"){
            string colon;
            ss >> colon >> MacroCount;
        }else{
            string node_name = token;
            float w,h;
            ss >> w >> h;
            //构造Module
            Module* tp_mod = new Module();
            tp_mod->name = node_name;
            tp_mod->width = w;
            tp_mod->height = h;
            tp_mod->idx = Nodes.size() + Terminals.size();
            string type = "";
            ss >> type;
            if(type == "terminal"){
                tp_mod->isMacro = true;
                tp_mod->isFixed = true;
                Terminals.push_back(tp_mod);
            }else{
                Nodes.push_back(tp_mod);
            }
            moduleMap[node_name] = tp_mod;
        }
    }
    ifs.close();
}
void PlaceData::readPl(const string& filename){
    ifstream ifs(filename);
    if(!ifs.is_open()){
        cerr << "failed to open file: " << filename << endl;
        return;
    }
    string line;
    while(getline(ifs,line)){
        if(line.empty() || line[0] == '#' || line[0] == 'U'){
            continue;
        }
        stringstream ss(line);
        string node_name;
        float x,y;
        ss >> node_name >> x >> y;
        auto it = moduleMap.find(node_name);
        if(it != moduleMap.end()){
            it->second->center.x = x;
            it->second->center.y = y;
        }else{
            cerr << "unknown node name in .pl: " << node_name << endl;
        }
    }
    ifs.close();
}
void PlaceData::readScl(const string& filename){
    ifstream ifs(filename);
    if(!ifs.is_open()){
        cerr << "failed to open file: " << filename << endl;
        return;
    }
    string line;
    SiteRow sr;
    int numSites = 0;
    while(getline(ifs,line)){
        if(line.empty() || line[0] == '#' || line[0] == 'U'){
            continue;
        }
        stringstream ss(line);
        string token;
        string colon;
        ss >> token;
        if(token == "CoreRow"){
            sr = SiteRow();
        }else if(token == "Coordinate"){
            ss >> colon >> sr.bottom;
            sr.start.y = sr.bottom;
            sr.end.y = sr.bottom;
        }else if(token == "Height"){
            ss >> colon >> sr.height;
        }else if(token == "Sitespacing"){
            ss >> colon >> sr.step;
        }else if(token == "SubrowOrigin"){
            string ns;
            ss >> colon >> sr.start.x >> ns;
            ss >> colon >> numSites;
            sr.end.x = sr.start.x + numSites * sr.step;
        }else if(token == "End"){
            SiteRows.push_back(sr);
        }
    }
    ifs.close();
}
void PlaceData::readBookShelf(const string& dir, const string& benchmarkname){
    //dir = ./test_data/adaptec1    benchmark = adaptec1
    string file_path_prefix = dir + "/" + benchmarkname;
    cout << "开始解析此路径文件：" << dir << "...\n";
    //开始解析BookShelf
    readNodes(file_path_prefix + ".nodes");
    readPl(file_path_prefix + ".pl");
    readScl(file_path_prefix + ".scl");
    //readNets(file_path_prefix + ".nets");
    cout << "[数据解析完成] 总普通单元数: " << Nodes.size() << ", 宏单元/IO数: " << Terminals.size() << endl;
    cout << "[数据解析完成] 总标准单元拜访行数: " << SiteRows.size() << endl;
}
