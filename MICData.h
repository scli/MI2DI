#ifndef _MIC_DATA_H_
#define _MIC_DATA_H_

typedef unsigned char DisDataT;
typedef double        CtnDataT;
enum DataType {_CONTINUOUS_, _DISCRETE_};

double maximum(double a, double b);
/*
class DataType
{
};

class DisDataType : DataType
{
    public:  
      DistDataT mData;
};

class CtnDataType : DataType
{
   public:
     CtnDataTy  mData;
};

*/

class MICData
{

    protected: 
     int mNumCols;
     int mNumRows;
     int mNumHiddenClus; 
     int mNumObserveClus;

     int* mHiddenOffset;
     int* mHiddenDim;
     int* mObserveOffset;
     int* mObserveDim;
     DisDataT** mRawData;
     double** mJointProbMatrix;
     double** mMargin;

     double** mPseudoJointProbMatrix;
     double** mPseudoMargin;
     double** mInvMat;

     double** mCovMat;
     DataType mDataType; 
    public:
     MICData();
     ~MICData() {; }

     int getNumCols() { return mNumCols; }
     int getNumRows() { return mNumRows; }
     DataType getDataType() { return mDataType; } 
     int getObserveOffset(int i) { return mObserveOffset[i]; }
     int getObserveDim(int i)    { return mObserveDim[i];    }
     int getHiddenOffset(int i)  { return mHiddenOffset[i];  }
     int getHiddenDim(int i)     { return mHiddenDim[i];     }
     double mutualInformation(int o1, int o2, double** joint,
                              double* margin1, double* margin2);
     int getNumHiddenClus()      { return mNumHiddenClus;  }
     int getNumObserveClus()     { return mNumObserveClus; }

     int getMaxHiddenDim() { return 150;      }
     DisDataT** getRawData() { return mRawData; }

     void setJointMatrix(double** mat) { mJointProbMatrix=mat;    }
     void setMargin(double** margin)   { mMargin         =margin; }
     double** getJointMatrix() { return mJointProbMatrix; }
     double** getMargin()      { return mMargin;          }

     void setPseudoJointMatrix(double** mat) { mPseudoJointProbMatrix=mat;    }
     void setPseudoMargin(double** margin)   { mPseudoMargin         =margin; }
     double** getPseudoJointMatrix() { return mPseudoJointProbMatrix; }
     double** getPseudoMargin()      { return mPseudoMargin;          }
     
     void     setInverseMat(double** invMat) {  mInvMat=invMat; }
     double** getInverseMat()        { return mInvMat; }

     void     setCovMat(double** covmat) { mCovMat=covmat; }
     double** getCovMat()                { return mCovMat; } 

     virtual double similarity(DisDataT* data1, DisDataT* data2) {; }
     virtual double similarity(int index1, int index2) {; }
     virtual int    getIndex  (int index, DisDataT data) {; }
     virtual int    getIndex  (int index, int row) { ; }
     virtual double getValue  (int index, int row) { ; }
     virtual double localSingularSimilarity(int pos,  int a1, int b1, double identity) { ;}
     virtual double localPairSimilarity(int pos1, int pos2, int a1, int a2, int b1, int b2, double identity) {;}
     virtual double lookupK22(int pos1, int pos2, int aa1, int aa2, int bb1, int bb2) {  if(aa1==bb1 && aa2==bb2) return 1; return 0;}
     virtual double lookupK1(int pos, int aa1, int bb1) { if(aa1==bb1) return 1; return 0; }
};

#endif
