#include "myplace.h"
#include "placedata.h"
//#include "plot.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
    PlaceData *placedata = new PlaceData();
    //MyPlacer *myplacer = new MyPlacer(placedata);
    placedata->readBookShelf("./test_data/adaptec1","adaptec1");
    return 0;
}
