#ifndef _DATA_PARTITION_H_
#define _DATA_PARTITION_H_
#include "DiscrData.h"
#include "CntnsData.h"
#include <vector>
using namespace std;

class DataPartition 
{
    DiscrData* mDiscrData;
    CntnsData* mCntnsData;
    double  ** mBounds;
    int     *  mNumClusters;
    int  mMaxNumClusters;

    void uniformPartition();
    void setDiscrData();
    void setCntnsData();

    CtnDataT** mRawData;
    DisDataT** mCluData;
    double** mK22Table;
    double** mK1Table;
    int mNumRows; 
    int mNumCols;
   public:

    DataPartition(double** rawData, int n_row, int n_col);
    ~DataPartition();
    void partition(); 
    
    DiscrData* getDiscrData() { return mDiscrData; }
    CntnsData* getCntnsData() { return mCntnsData; }

};

#endif
