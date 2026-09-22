#include "placedata.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

void PlaceData::readNodes(const string& filename){
    //读入
    ifstream ifs(filename);
    if(!ifs.is_open()){
        cerr << "fail to open fiel: " << filename << endl;
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
void PlaceData::readBookShelf(const string& dir, const string& benchmarkname){
    //dir = ./test_data/adaptec1    benchmark = adaptec1
    string file_path_prefix = dir + "/" + benchmarkname;
    cout << "开始解析此路径文件：" << dir << "...\n";
    //开始解析BookShelf
    readNodes(file_path_prefix + ".nodes");
    //readPl(file_path_prefix + ".pl");
    //readScl(file_path_prefix + ".scl");
    //readNets(file_path_prefix + ".nets");
    cout << "[数据解析完成] 总普通单元数: " << Nodes.size() << ", 宏单元/IO数: " << Terminals.size() << endl;
}
