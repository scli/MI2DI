#include "DiscrData.h"
#include "MICParaMgr.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

DiscrData::DiscrData(int numRows, int numCols, int k22Dim) : MICData()
{
   mNumRows = numRows;
   mNumCols = numCols;
   mK22Dim  = k22Dim;
   mObserveDim    =new int[numCols];
   mHiddenDim     =new int[numCols];
   mObserveOffset =new int[numCols];
   mHiddenOffset  =new int[numCols];
   mDataType =_DISCRETE_;
}
void
DiscrData::setData(DisDataT** data)
{
   mRawData=data;
}

void
DiscrData::setK22Table(double** table)
{
  mK22Table=table;  
}

void
DiscrData::initK1Table(double** table)
{
  mK1Table=table;
}

void
DiscrData::setObserveDim(int index, int num)
{
   mObserveDim[index]=num;
}
void
DiscrData::setHiddenDim(int index, int num)
{
   mHiddenDim[index]=num;
}


void
DiscrData::setObserveOffset(int index, int num)
{
   mObserveOffset[index]=num;
}
void
DiscrData::setHiddenOffset(int index, int num)
{
   mHiddenOffset[index]=num;
}

void
DiscrData::setNumHiddenClusters(int num)
{
   mNumHiddenClus  = num;
}
void
DiscrData::setNumObserveClusters(int num)
{
  mNumObserveClus=num;
}

//
//
//=======================================

double
DiscrData::similarity(DisDataT* data1, DisDataT* data2)
{
    //this similarity is non sysmemtric
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
DiscrData::similarity(int index1, int index2)
{
    //this similarity is non sysmemtric
  DisDataT* data1=mRawData[index1];
  DisDataT* data2=mRawData[index2];
  return similarity(data1, data2);
}




int
DiscrData::getIndex(int pos, DisDataT data)
{
  //if 21
    return data;
}

int
DiscrData::getIndex(int pos, int row)
{
  //if 21
    return mRawData[row][pos];
}



double
DiscrData::getValue(int pos, int row)
{
   return 1;
}
double
DiscrData::localSingularSimilarity(int pos, int a1, int b1, double identity)
{
   double theta=0.7;
   if(identity<theta) return 0;
   return 1; 
}


double 
DiscrData::localPairSimilarity(int pos1, int pos2,  int a1, int a2, int b1, int b2, double identity)
{

    
   double theta=0.7;
   double tiny=1.0e-100; 
   if(identity<theta) return 0;
   return 1;
   
}



double
DiscrData::lookupK22(int pos1, int pos2, int aa1, int aa2, int bb1, int bb2)
{
   return  mK22Table[mK22Dim*aa1+aa2][mK22Dim*bb1+bb2];
}

double
DiscrData::lookupK1(int pos, int aa1,int bb1)
{
   return  mK1Table[aa1][bb1];
}

