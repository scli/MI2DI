#include "CntnsData.h"
#include "MICParaMgr.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

CntnsData::CntnsData(int numRows, int numCols) : MICData()
{
   mNumRows = numRows;
   mNumCols = numCols;
   mObserveDim    =new int[numCols];
   mHiddenDim     =new int[numCols];
   mObserveOffset =new int[numCols];
   mHiddenOffset  =new int[numCols];
   mDataType=_CONTINUOUS_;
}

//initialize raw data
void
CntnsData::setRawData(CtnDataT** data)
{
   mRawData=data;
}


void
CntnsData::setIndex(DisDataT** indices)
{
   mIndices=indices;
}

void
CntnsData::setObserveDim(int index, int num)
{
   mObserveDim[index]=num;
}
void
CntnsData::setHiddenDim(int index, int num)
{
   mHiddenDim[index]=num;
}


void
CntnsData::setObserveOffset(int index, int num)
{
   mObserveOffset[index]=num;
}
void
CntnsData::setHiddenOffset(int index, int num)
{
   mHiddenOffset[index]=num;
}

void
CntnsData::setNumHiddenClusters(int num)
{
   mNumHiddenClus  = num;
}
void
CntnsData::setNumObserveClusters(int num)
{
  mNumObserveClus=num;
}

//
//
//=======================================

double
CntnsData::similarity(CtnDataT* data1, CtnDataT* data2)
{
    //this similarity is non sysmemtric
  if(data1==data2)
  {
     return 1;
  }
  return 0;
    
  double rev=0;
  double effectiveIndex=0;
  for(int i=0; i<mNumCols; i++)
  {
     if(data1[i]== data2[i])
      effectiveIndex+=1;
  }
  if(effectiveIndex<=0) return 0;
  return (effectiveIndex/mNumCols); 
}



double
CntnsData::similarity(int index1, int index2)
{
    //this similarity is non sysmemtric
    if(index1==index2)
      return 1;
    return 0;
}



int
CntnsData::getIndex(int pos, int row)
{
  //if 21
    return mIndices[row][pos];
}

double
CntnsData::getValue(int pos, int row)
{
   return mRawData[row][pos];
}

double
CntnsData::localSingularSimilarity(int pos, int a1, int b1, double identity)
{
   double theta=0.7;
   if(identity<theta) return 0;
   return 1; 
}


double 
CntnsData::localPairSimilarity(int pos1, int pos2,  int a1, int a2, int b1, int b2, double identity)
{
   double theta=0.7;
   double tiny=1.0e-100; 
   if(identity<theta) return 0;
   return 1;
}

//=======================================
