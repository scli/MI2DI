#ifndef _JOINT_PROB_MATRIX_H_
#define _JOINT_PROB_MATRIX_H_

#include "MICData.h"

class JointProbMatrix
{
        MICData* mData;
        DisDataT** mRawData;
        
        
        double** mJointMat; //dot product of two columns
        double** mMarginMean;  //mean of the margin
        double** mMarginVariance;//variances of the margin


        // those are temp vars
        double** mNumPairMatrix;
        double** mNumSingular;
   
        double** mWeights;
        double*  mSingularWeights;

        void computeSmallProbMatrix(int pos1, int pos2);
        void countNumPairs();
        void allocateSpace();
        void convertCountsToProbMatrix();
        void pairAndSingularWeights(int row, double* weights, double* singularWeights);
        double getProbs(int col_id1, int col_id2, int hidden1, int hidden2);
        void computeProbCols();
        void countSanderNumPairs();
        void countSimpleContinuousNumPairs();

        double getProbs(int col_id, int hidden_clu);
    public:
        JointProbMatrix(MICData* data);
        ~JointProbMatrix();
        double** getJointMatrix() { return mJointMat; }
        double** getMargin()      { return mMargin;       }


};

#endif
