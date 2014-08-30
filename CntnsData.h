#ifndef _CNTNS_DATA_H_
#define _CNTNS_DATA_H_

#include "MICData.h"
class CntnsData : public MICData
{

    void  normalize(double** vect, int size);
    DisDataT** mIndices;
    CtnDataT** mRawData;
    double**   mBounds;
    double**   mInvCorrMat;

  public:
    CntnsData(int numRow, int numCol);
    ~CntnsData();
   
    virtual double similarity(CtnDataT* data1, CtnDataT* data2);
    virtual double similarity(int index1, int index2);
    virtual int    getIndex  (int index, DisDataT data) { return data; }
    virtual double getValue  (int index, int which);
    virtual double localSingularSimilarity(int pos,  int a1, int b1, double identity);
    virtual double localPairSimilarity(int pos1, int pos2, int a1, int a2, int b1, int b2, double identity);
    virtual int    getIndex(int pos, int row);
    void     setRawData(CtnDataT** data);
    void     setIndex(DisDataT** indices);
    void     setData(DisDataT** data);
    void     setObserveDim(int index, int num);
    void     setHiddenDim(int index, int num);
    void     setObserveOffset(int index, int num);
    void     setHiddenOffset(int index, int num);
    void     setNumHiddenClusters(int num);
    void     setNumObserveClusters(int num);
    void     setBounds(double** bounds) { mBounds=bounds; }
    void     setInvCorrMat(double** mat) { mInvCorrMat=mat; }
    double** getInvCorr() { return mInvCorrMat; }
};



#endif
