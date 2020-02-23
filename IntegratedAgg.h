#ifndef _INTEGRADED_AGG_H_
#define _INTEGRADED_AGG_H_

#include "MICData.h"
#include "CntnsData.h"
#include "DiscrData.h"
#include <vector>
#include "Aggregation.h"
#include "SanderAgg.h"
using namespace std;

class IntegratedAgg : public Aggregation
{
     double** mDisPDir;
     double** mDisMargin;
     double** mCtnInvCorr;
     double aggregate(int pos1, int pos2);

     JointMat* mMatrix; 
    public:
        IntegratedAgg(DiscrData* disData, CntnsData* conData, JointMat* mat, double** inv);
        ~IntegratedAgg();
        void aggregate();

};
#endif
