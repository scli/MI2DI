#ifndef _SANDER_AGG_H_
#define _SANDER_AGG_H_

#include "MICData.h"
#include <vector>
#include "Aggregation.h"
#include "JointMat.h"
using namespace std;

class SanderAgg : public Aggregation
{
        double** mAggregatedMat;
        double** mInvMat;
        double** mMargin;
         

        double* mMu1;
        double* mMu2;
        double* mScale1;
        double* mScale2;
        double* new1;
        double* new2;
        double** mPDir;
        double** mPFac;
        double** mWMatrix; 

        double aggregate(int pos1, int pos2);

        void getMatrixW(int pos1, int pos2, double** w);
        void aggregateByDir();
        void aggregateByMax();
        void aggregateByAverage();
        void aggregateByMinEdge();
        //findMinimum(int v1, int v2);
        void directInformation();
        void allocateSpace4Sander();
        void computeMU(int pos1, int pos2);
        double computeDI(int pos1, int pos2);
	JointMat* mMatrix;

    public:
        SanderAgg(MICData* data, JointMat* mat, double ** inv);
        ~SanderAgg();
        void aggregate();
        double** getPDir() { return mPDir; }
        double** getMargin() { return mMargin; }


};
#endif
