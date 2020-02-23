#include "IntegratedAgg.h"
#include "MICData.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

IntegratedAgg::IntegratedAgg(DiscrData* disData, CntnsData* conData, JointMat* mat, double** inv) 
    : Aggregation(disData)
{

  SanderAgg* sanderAgg=new SanderAgg(disData, mat, inv);
  mMatrix      =mat; 

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
  for(int i=0; i<mMatrix->numCols(); i++)
  {
    for(int j=0; j<mMatrix->numCols(); j++)
    {
       mDIMatrix[i][j] =aggregate(i, j);
    }
  }
}



double
IntegratedAgg::aggregate(int pos1, int pos2)
{
   double rev=0;

   int dim1=mMatrix->dim(pos1);
   int off1=mMatrix->offset(pos1);

   int dim2=mMatrix->dim(pos2);
   int off2=mMatrix->offset(pos2);
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

