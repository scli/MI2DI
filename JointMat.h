#ifndef _JOINT_MATRIX_H_
#define _JOINT_MATRIX_H_

#include "MICData.h"
//#include "LatentMap.h"
class JointMat
{
        MICData* mData;
        
        int mNumCol;
	int mTotSize;
        double** mJoint; //dot product of two columns
        double**   mMargin;  //mean of the margin


        // those are temp vars
   
        double** mWeights;
        double*  mSingularWeights;

        double*  mReplicateWeight;
	int*     mOffset;
	int*     mDim; 


        void allocateSpace();
        
        void computeReplicateWeight();
        
        void countNumPairs();
        void sanderCountNumPairs();
        void localCountNumPairs();

        void computeJointMean();
        void computeSubMatrixMean(int pos1, int pos2);
        void pairAndSingularWeights(int row, double* weights, double* singularWeights);
        double joint(int col_id1, int col_id2, int hidden1, int hidden2);

        void computeColMeans();
        double margin(int col_id, int hidden_clu);

    public:
        JointMat(MICData* data);
	JointMat(int totSize, int size, int* offset, int* dim);
        JointMat(int totSize, int size, int* offset, int* dim, double** m, double** j);
	~JointMat();
        double** joint()    { return mJoint; }
        double** margin()      { return mMargin;   }
        double get(int pos, int dim);
	void   set(int pos, int dim, double value);
	double get(int pos1, int dim1, int pos2, int dim2);
	void   set(int pos1, int dim1, int pos2, int dim2, double value);
	int    numCols() { return mNumCol;  }
        int    numClus() { return mTotSize; }
        int    offset(int i) { return mOffset[i]; }
	int    dim(int i)    { return mDim[i];    }
	int*   offset()  { return mOffset; }
	int*   dim()     { return mDim;    }


};

#endif
