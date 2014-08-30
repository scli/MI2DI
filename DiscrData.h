#ifndef _DISCR_DATA_H_
#define _DISCR_DATA_H_

#include "MICData.h"
class DiscrData : public MICData
{
  protected:
    int mK22Dim;
    double** mK22Table;
    double** mK1Table;
  public:
    DiscrData(int numRow, int numCols, int k22Dim);
    ~DiscrData();
   
    virtual double similarity(DisDataT* data1, DisDataT* data2);
    virtual double similarity(int index1, int index2);
    virtual int    getIndex  (int index, DisDataT data);
    virtual int    getIndex  (int index, int row);
    virtual double getValue  (int index, int row);
    virtual double localSingularSimilarity(int pos,  int a1, int b1, double identity);
    virtual double localPairSimilarity(int pos1, int pos2, int a1, int a2, int b1, int b2, double identity);
    virtual double lookupK22(int pos1, int pos2, int aa1, int aa2, int bb1, int bb2);
    virtual double lookupK1(int pos, int aa1, int bb1);
    void    setData(DisDataT** data);
    void    setK22Table(double** table);
    void    initK1Table(double** table);
    void    setObserveDim(int index, int num);
    void    setHiddenDim(int index, int num);
    void    setObserveOffset(int index, int num);
    void    setHiddenOffset(int index, int num);
    void    setNumHiddenClusters(int num);
    void    setNumObserveClusters(int num);
   // void initializeData(AlignmentParser* p);

};
#endif
