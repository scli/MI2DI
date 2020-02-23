#include "CovMat.h"
#include <iostream>
using namespace std;

CovMat::CovMat(JointMat* aData)
{
   mData=aData;
   mMargin=mData->margin();
   mJoint=mData->joint();
   calculateCov();
}



void
CovMat::calculateCov()
{

   int tot_size = mData->numClus();
   int num_cols = mData->numCols();

   mCovMat=new double* [tot_size];

   for(int i=0; i<tot_size; i++)
   {
      mCovMat[i]=new double[tot_size];
   }//for

   for(int i=0; i<num_cols; i++)
   {
       int dim1=mData->dim(i);
       int off1=mData->offset(i);

       for(int j=0; j<num_cols; j++)
       {
         int dim2=mData->dim(j);
         int off2=mData->offset(j);

         for(int k1=0; k1<dim1; k1++)
         {
            for(int k2=0; k2<dim2; k2++)
            {
               mCovMat[off1+k1][off2+k2]=mJoint[off1+k1][off2+k2]-mMargin[i][k1]*mMargin[j][k2];
               //cout<<mJointMat[off1+k1][off2+k2]<<" "<<mMarginMean[i][k1]*mMarginMean[j][k2]<<" lala"<<endl;; 
               if(mCovMat[off1+k1][off2+k2]!=mCovMat[off1+k1][off2+k2])
                  cout<<"error in covarmatrix "<<endl;//=mJointMat[off1+k1][off2+k2]-mMarginMean[i][k1]*mMarginMean[j][k2];
            }
         }
       }
    }
}
