#ifndef _JOINT_MATRIX_H_
#define _JOINT_MATRIX_H_

#include "MICData.h"

class JointMat
{
        MICData* mData;
        
        
        double** mJointMean; //dot product of two columns
        double**   mColMean;  //mean of the margin


        // those are temp vars
        double** mJointSum;
        double** mColSum;
   
        double** mWeights;
        double*  mSingularWeights;

        double*  mReplicateWeight;

        void allocateSpace();
        
        void computeReplicateWeight();
        
        void countNumPairs();
        void sanderCountNumPairs();
        void localCountNumPairs();

        void computeJointMean();
        void computeSubMatrixMean(int pos1, int pos2);
        void pairAndSingularWeights(int row, double* weights, double* singularWeights);
        double getJointMean(int col_id1, int col_id2, int hidden1, int hidden2);

        void computeColMeans();
        double getColMean(int col_id, int hidden_clu);

    public:
        JointMat(MICData* data);
        ~JointMat();
        double** getJointMean()    { return mJointMean; }
        double** getColMean()      { return mColMean;   }


};

#endif
