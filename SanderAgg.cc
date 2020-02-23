#include "SanderAgg.h"
#include "MICData.h"
#include "JointMat.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

SanderAgg::SanderAgg(MICData* data, JointMat* mat, double** inv) : Aggregation(data)
{
  mMatrix=mat; 
  mMargin=mat->margin();
  mInvMat=inv;

  directInformation();
}


SanderAgg::~SanderAgg()
{


}
    
void
SanderAgg::aggregate()
{
  for(int i=0; i<mMatrix->numCols(); i++)
  {
    for(int j=0; j<mData->numCols(); j++)
    {
       mDIMatrix[i][j] =aggregate(i, j);
    }
  }
}

double
SanderAgg::aggregate(int pos1, int pos2)
{
   double rev=0;

   int dim1=mMatrix->dim(pos1);
   int off1=mMatrix->offset(pos1);

   int dim2=mMatrix->dim(pos2);
   int off2=mMatrix->offset(pos2);

   for(int i=0; i<dim1; i++)
   {
      for(int j=0; j<dim2; j++)
      {
         rev+=fabs(mInvMat[off1+i][off2+j]); 
      }
   }
   return rev;
}


/*
void
SanderAgg::aggregateByDir()
{

}

void
SanderAgg::aggregateByMax()
{

}

void
SanderAgg::aggregateByAverage()
{

}
*/

//=================================
/*void
SanderAgg::aggregateByMinEdge()
{
  for(int i=0; i<mData->getNumCols(); i++)
  {
    for(int j=0; j<mData->getNumCols(); j++)
    {
       mAggregatedMat[i][j] =findMinimum(i, j)
    }
  }
}*/
/*
DataT
SanderAgg::findMinimum(int v1, int v2)
{
   int len1 = mLens[v1];
   int len2 = mLens[v2];
   int pos1 = mOffset[v1];
   int pos2 = mOffset[v2];
   DataT rev=100000;
   for(int k1=pos1; k1<pos1+len1; k1++)
   {
     for(int k2=pos2; k2<pos2+len2; k2++)
     {
         if(mData[k1][k2]<rev)
            rev = mData[k1][k2];
     }
   }
   return rev;
}*/

void
SanderAgg::allocateSpace4Sander()
{
   int num_clusters=mData->getMaxHiddenDim();
   int tot_size    =mMatrix->numClus();
   mMu1    =new double[num_clusters];
   mMu2    =new double[num_clusters];
   mScale1 =new double[num_clusters];
   mScale2 =new double[num_clusters];
  
   new1    =new double[num_clusters];
   new2    =new double[num_clusters];
  
   mWMatrix=new double* [num_clusters];
   mPFac   =new double* [num_clusters];


   mPDir   =new double* [tot_size];

   for(int i=0; i<tot_size; i++)
   {
      mPDir[i]   =new double[tot_size];
   }

  
   for(int i=0; i < num_clusters; i++)
   {
     mWMatrix[i]=new double[num_clusters];
     mPFac[i]   =new double[num_clusters];
   }

   mPDir   =new double* [tot_size];

   for(int i=0; i<tot_size; i++)
   {
      mPDir[i]   =new double[tot_size];
   }
 
}
void
SanderAgg::directInformation()
{
  allocateSpace4Sander(); 
  for(int i=0; i<mMatrix->numCols(); i++)
  {
      mDIMatrix[i][i]=-1000;
      for(int j=i+1; j<mData->numCols(); j++)
      {
         getMatrixW(i, j, mWMatrix); 
         //fillWMatrix(mWMatrix,ALPHABET_SIZE);
         computeMU(i, j);
         mDIMatrix[i][j]=computeDI(i, j);
         if(mDIMatrix[i][j]!=mDIMatrix[i][j])
         {
            cerr<<"error identified in MI"<<endl;
            exit(0);
         }
         mDIMatrix[j][i]=mDIMatrix[i][j];
         //cerr<<mDIMatrix[j][i]<<" direct"<<endl;
      }
  }
}

/*
void
SanderAgg::sanDiago()
{

   int num_clusters=mData->getMaxHiddenDim();
   int tot_size    =mData->getNumHiddenClus();
   int num_col     =mData->getNumCols();

	for (int i=0; i<tot_size; i++) 
    {
	   for (int j=i; j<tot_size; j++) 
       {
           double sum=0;
			for (k=0; k<ndim; k++) {
				sum += mInvMat[i][k] * mCotMat[k][j];
			}
			beta[i][j] = sum;
            beta[j][i] = sum;
		}		
	}		

    pcor = allocmat(seqlen, seqlen, sizeof(double));
    pcor_apc = allocmat(seqlen, seqlen, sizeof(double));
    marg_di = allocvec(seqlen, sizeof(double));
    
    for (int i=0; i<num_col; i++) 
    {
		for (int j=0; j<num_col; j++)
        {
		     pcor[i][j] = 0;
		}
	}
	
    p=0;
    for (int i=0; i<num_col; i++) 
    {
       int dim1=mData->dim(i);
       int off1=mData->offset(i);
       for (int j=0; j<num_col; j++) 
       {
          int dim2=mData->dim(j);
          int off2=mData->offset(j);
          
          pcor[i][j] =0;
          for(int k1=0; k1<dim1; k1++)
          {
            for(int k2=0; k2<dim2; k2++)
            {
		    	pcor[i][j] = fabs(mInvMat);	
		    }
		  }
          if(i==j) pcor[i][j]=1;
		}	    
	  }
	 
     di_mean=0;
	 for (i=0; i<num_col; i++) 
     {
		sum_di=0;
		for (j=0; j<num_col; j++) 
        {
            if (i==j)
			  continue;              
			sum_di += pcor[i][j];
        }
		marg_di[i] = sum_di/(num_col-1);
			
		for (j=i+1; j<num_col; j++)
        {
		   di_mean += pcor[i][j];
		}
	  }
	  di_mean = (2*di_mean)/(num_col*(num_col-1));
 
	  for (i=0; i<num_col; i++) {
		 for (j=0; j<num_col; j++) {
		 	 if (i==j) {
		 	 	pcor_apc[i][j] = 1;
		 	 }
		 	 else {
			    sum_di = (marg_di[i]*marg_di[j])/di_mean;
			    if (apcflg) {
			        pcor_apc[i][j] = pcor[i][j] - sum_di;
			    }
			    else {
			    	pcor_apc[i][j] = pcor[i][j];
			    }
			 }
		 }
	  }


}*/



void
SanderAgg::getMatrixW(int pos1, int pos2, double** w)
{
   int dim1=mMatrix->dim(pos1);
   int dim2=mMatrix->dim(pos2);

   int off1=mMatrix->offset(pos1);
   int off2=mMatrix->offset(pos2);

    
   for(int i=0; i<dim1; i++)
   {
      for(int j=0; j<dim2; j++)
      {
         w[i][j]=0;
      }
   }
   
   for(int i=0; i<dim1; i++)
   {
      double sum=0;
      for(int j=0; j<dim2; j++)
      {
         w[i][j]=mInvMat[off1+i][off2+j];
      }
   }

   for(int j=0; j<dim1; j++)
   {
      for(int i=0; i<dim2; i++)
      {
         w[i][j]=exp(-w[i][j]);
         if(w[i][j]!=w[i][j])
             cout<<"errors delected from the inv matrix"<<endl;
      }
   }
}
   


void
SanderAgg::computeMU(int pos1, int pos2)
{

    
  double epsilon=1e-4;
  double diff =1.0;
  
  int dim1=mMatrix->dim(pos1);
  int dim2=mMatrix->dim(pos2);

  for(int i=0; i<dim1; i++)
  {
     mMu1[i]=1.0/dim1;
  }
  for(int i=0; i<dim2; i++)
  {
     mMu2[i]=1.0/dim2;
  }

  //pi = P1(i, :);
  //pj = P1(j, :);
  //cerr<<"entering a while loop"<<endl;
  while (diff > epsilon)
  {
     vectorTimesTranspose(mMu2, mWMatrix, dim1, dim2, mScale1);
        vectorTimesMatrix(mMu1, mWMatrix, dim1, dim2, mScale2);
     
     innerDivide(mMargin[pos1], mScale1, new1, dim1);
     normalize(new1, dim1);

     innerDivide(mMargin[pos2], mScale2, new2, dim2);
     normalize(new2, dim2);

     double a=maximum(new1, mMu1, dim1);
     double b=maximum(new2,mMu2, dim2);
     if(a>b) diff=a; else diff=b;
	 //diff = maximum(maximum(new1, mMu1, dim1), maximum(new2,mMu2, dim2));

     for(int i=0; i<dim1; i++)
     {
        mMu1[i]=new1[i];
     }

     for(int i=0; i<dim2; i++)
     {
        mMu2[i]=new2[i];
     }


     //cerr<<diff<<endl;
     //mu1 = new1;
     //mu2 = new2;
 
  }
  //cerr<<"exit the while loop"<<endl;
}
double
SanderAgg::computeDI(int pos1, int pos2)
{
  double tiny = 1.0e-100;
  int dim1=mMatrix->dim(pos1);
  int dim2=mMatrix->dim(pos2);

  int off1=mMatrix->offset(pos1);
  int off2=mMatrix->offset(pos2);



  for(int i=0; i<dim1; i++)
  {
     for(int j=0; j<dim2; j++)
     {
       mPDir[off1+i][off2+j]=mWMatrix[i][j]*(mMu1[i]*mMu2[j]);
     }
  }
  //cout<<"==============================="<<endl;
  double sum=0;
  for(int i=0; i<dim1; i++)
  {
     for(int j=0; j<dim2; j++)
     {
       sum+=mPDir[off1+i][off2+j];
     }
  }

  for(int i=0; i<dim1; i++)
  {
     for(int j=0; j<dim2; j++)
     {
       mPDir[off1+i][off2+j]/=sum;
     }
  }

  return mData->mutualInformation(pos1, pos2, mPDir, mMargin[pos1], mMargin[pos2]);
}


