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
            tp_mod->area = w * h;
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
        string colon;
        string orientation;
        string fixed_marker;
        float x,y;
        ss >> node_name >> x >> y >> colon >> orientation;
        bool is_fixed = false;
        if(ss >> fixed_marker){
            if(fixed_marker != "/FIXED"){
                cerr << "unknown placement flag in .pl for node " << node_name << ": " << fixed_marker << endl;
                continue;
            }
            is_fixed = true;
        }
        auto it = moduleMap.find(node_name);
        if(it != moduleMap.end()){
            it->second->center.x = x;
            it->second->center.y = y;
            it->second->orientation = 0;
            it->second->isFixed = is_fixed;
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
void PlaceData::readNets(const string& filename){
    ifstream ifs(filename);
    if(!ifs.is_open()){
        cerr << "failed to open file: " << filename << endl;
        return;
    }
    string line;
    Net* nt = nullptr;
    int pins_to_read = 0;
    while(getline(ifs,line)){
        if(line.empty() || line[0] == '#' || line[0] == 'U'){
            continue;
        }
        stringstream ss(line);
        string token;
        string colon;
        ss >> token;
        if(token == "NumNets"){
            ss >> colon >> netCount;
        }else if(token == "NumPins"){
            ss >> colon >> pinCount;
        }else if(token == "NetDegree"){
            string net_name;
            int degree;
            ss >> colon >> pins_to_read >> net_name;
            nt = new Net();
            nt->idx = Nets.size();
            nt->name = net_name;
            Nets.push_back(nt);
            netMap[net_name] = nt;
        }else{
            if(pins_to_read > 0){
                string node_name = token;
                string direction;
                float x,y;
                ss >> direction >> colon >> x >> y;
                auto it = moduleMap.find(node_name);
                if(it != moduleMap.end()){
                    Module* tp_mod = it->second;
                    Pin* p = new Pin();
                    p->idx = Pins.size();
                    p->module = tp_mod;
                    p->net = nt;
                    p->offset.x = x;
                    p->offset.y = y;
                    Pins.push_back(p);
                    nt->netPins.push_back(p);
                    tp_mod->modulePins.push_back(p);
                    tp_mod->nets.push_back(nt);
                }
            }
            --pins_to_read;
        }
    }
    ifs.close();
}
void PlaceData::readWts(const string& filename){
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
        string object_name;
        float weight;
        ss >> object_name >> weight;
        if(ss.fail()){
            cerr << "invalid weight record in .wts: " << line << endl;
            continue;
        }
        auto module_it = moduleMap.find(object_name);
        if(module_it != moduleMap.end()){
            module_it->second->weight = weight;
            continue;
        }
        auto net_it = netMap.find(object_name);
        if(net_it != netMap.end()){
            net_it->second->weight = weight;
        }else{
            cerr << "unknown object name in .wts: " << object_name << endl;
        }
    }
    ifs.close();
}
bool PlaceData::readAux(
    string& nodes_filepath,
    string& pl_filepath, 
    string& scl_filepath, 
    string& nets_filepath, 
    string& wts_filepath,
    const string& aux_filepath
){
    //读取aux文件来获取其他文件路径
    ifstream ifs(aux_filepath);
    if(!ifs.is_open()){
        cerr << "failed to open file: " << aux_filepath << endl;
        return false;
    }
    string line;
    while(getline(ifs,line)){
        if(line.empty()){
            continue;
        }
        stringstream ss(line);
        string token,colon;
        ss >> token;
        if(token == "RowBasedPlacement"){
            ss >> colon;
            if(colon != ":"){
                cerr << "expected ':' after RowBasedPlacement in aux file" << endl;
                return false;
            }
            ss >> token;
            while(!token.empty()){
                string tp;
                for(int j = 0; j < (int)token.size(); ++j){
                    if(token[j] == '.'){
                        tp = token.substr(j);
                    }
                }
                if(tp == ".nodes"){
                    nodes_filepath += token;
                }else if(tp == ".pl"){
                    pl_filepath += token;
                }else if(tp == ".scl"){
                    scl_filepath += token;
                }else if(tp == ".nets"){
                    nets_filepath += token;
                }else if(tp == ".wts"){
                    wts_filepath += token;
                }
                token.clear();
                ss >> token;
            }
        }
    }
    ifs.close();
    return true;
}
void PlaceData::readBookShelf(const string& dir, const string& benchmarkname){
    //dir = ./test_data/adaptec1    benchmarkname = adaptec1
    string file_path_prefix = dir + "/" + benchmarkname;
    cout << "开始解析此路径文件：" << dir << "...\n";
    //开始解析BookShelf
    string nodes_filepath = dir + "/";
    string pl_filepath = dir + "/";
    string scl_filepath = dir + "/";
    string nets_filepath = dir + "/";
    string wts_filepath = dir + "/";
    bool aux_read = readAux(
        nodes_filepath,
        pl_filepath, 
        scl_filepath, 
        nets_filepath, 
        wts_filepath,
        file_path_prefix + ".aux"
    );
    if(!aux_read){
        cerr << "failed to read aux file: " << file_path_prefix + ".aux" << endl;
        return;
    }
    readNodes(nodes_filepath);
    readNets(nets_filepath);
    readWts(wts_filepath);
    readScl(scl_filepath);
    readPl(pl_filepath);
}
