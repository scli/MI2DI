#ifndef _ALIGN_DATA_H_
#define _ALIGN_DATA_H_
#include "AlignmentParser.h"
#include "MICData.h"
#include "DiscrData.h"
#include "LatentMap.h"
//typedef unsigned char DisDataT;


class AlignData 
{

    AlignmentParser* mParser;
    int mNumCols;
    int mNumRows;

    //double** mK22Table; 
    //double** mK1Table;//[21][21]; 
    //double** mTransition;
    //double*  mTransMargin;
    DisDataT** mRawData;
    DiscrData* mData;
    LatentMap *mMap;
    int getAAIndex(char C);
    void initData();
    void init21Data();
    void initDataBlock();
   // void initK22Table();
   // void initK1Table();
    //double K22(int aa_i, int aa_j, int bb_i, int bb_j);
    //double K1(int a, int a_prime);
    //double lookupK1(int hidden,int obs);
    //double lookupK22(int h1, int h2, int obs1, int obs2);
    //void  allocateTransition();
    //void  computeTransition();
    //void  normalize(double** vect, int size);
    void setData();
  public:
    AlignData(AlignmentParser* parser);
    ~AlignData();
    DiscrData* getData() { return mData; }
    LatentMap* getMap()  { return mMap; }
   // void initializeData(AlignmentParser* p);
};


#endif
