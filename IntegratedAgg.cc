#include "IntegratedAgg.h"
#include "MICData.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

IntegratedAgg::IntegratedAgg(DiscrData* disData, CntnsData* conData) 
    : Aggregation(disData)
{

  SanderAgg* sanderAgg=new SanderAgg(disData);
  mDisPDir     =sanderAgg->getPDir();
  mDisMargin   =sanderAgg->getMargin();
  mCtnInvCorr  =conData->getInvCorr();
  aggregate();
}


IntegratedAgg::~IntegratedAgg()
{


}
    
void
IntegratedAgg::aggregate()
{
  for(int i=0; i<mData->getNumCols(); i++)
  {
    for(int j=0; j<mData->getNumCols(); j++)
    {
       mMIMatrix[i][j] =aggregate(i, j);
    }
  }
}



double
IntegratedAgg::aggregate(int pos1, int pos2)
{
   double rev=0;

   int dim1=mData->getHiddenDim(pos1);
   int off1=mData->getHiddenOffset(pos1);

   int dim2=mData->getHiddenDim(pos2);
   int off2=mData->getHiddenOffset(pos2);
   rev=0;
   for(int i=0; i<dim1; i++)
   {
      for(int j=0; j<dim2; j++)
      {
         double disMI=mDisPDir[off1+i][off2+j]*log(mDisPDir[off1+i][off2+j]/
                                                  (mDisMargin[pos1][i]*mDisMargin[pos2][j]));

         rev+=disMI-0.5*mDisPDir[off1+i][off2+j]*log(1-mCtnInvCorr[off1+i][off2+j]); 
      }
   }
   return rev;
}

