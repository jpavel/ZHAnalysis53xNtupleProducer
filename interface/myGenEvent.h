#ifndef __MYGENEVENT_HH__
#define __MYGENEVENT_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myGenEvent : public TObject {
public:

    myGenEvent() {
        ;
    }

    ~myGenEvent() {
        ;
    }

    //General



    double alphaQCD, alphaQED, qScale, weight;
    bool hasPDF, hasBinningValues;
    unsigned int signalProcessID;
    int id_First, id_Second;
    double scalePDF;
    double x_First, x_Second;
    double xPDF_First, xPDF_Second;
    int binningValueSize;
    double binningValue0;


    ClassDef(myGenEvent, 1)
};
#endif
