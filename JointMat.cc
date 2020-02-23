#include "JointMat.h"
#include <iostream>
#include <stdlib.h>
#include "MICParaMgr.h"
using namespace std;


JointMat::JointMat(MICData* data)
{
   mData     = data; 
   mOffset   = data->offset();
   mDim      = data->dim();
    
   mNumCol   = data->numCols();
   mTotSize  = data->numClus();


   allocateSpace();
   computeReplicateWeight();
   
   //count the joint sum
   countNumPairs();
   
   //take the average
   computeJointMean();

   //compute the column means, assume the matrix is sym 
   computeColMeans();
}

JointMat::JointMat(int totSize, int size, int* offset, int* dim)
{
   mTotSize=totSize;
   mNumCol =  size;
   mOffset =  offset;
   mDim    =  dim;
   
   allocateSpace();
   //computeReplicateWeight();
   
   //count the joint sum
   //countNumPairs();
   
   //take the average
   //computeJointMean();

   //compute the column means, assume the matrix is sym 
   //computeColMeans();

     
}


JointMat::JointMat(int totSize, int size, int* offset, int* dim, double** m, double** j)
{
    mTotSize=totSize;
    mNumCol =  size;
    mOffset = offset;
    mDim    = dim;
    mMargin =m;
    mJoint  =j;

}

JointMat::~JointMat()
{
}

double 
JointMat::get(int pos, int dim)
{
   return mMargin[pos][dim];
}

void   
JointMat::set(int pos, int dim, double value)
{
     mMargin[pos][dim]=value;
}

double 
JointMat::get(int pos1, int dim1, int pos2, int dim2)
{
    return mJoint[mOffset[pos1]+dim1][mOffset[pos2]+dim2];
}

void   
JointMat::set(int pos1, int dim1, int pos2, int dim2, double value)
{
    mJoint[mOffset[pos1]+dim1][mOffset[pos2]+dim2] = value;
}


//count number of pairs
//
//
//
//

void
JointMat::countNumPairs()
{
   if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   {
      cout<<"I enter here"<<endl;	   
      sanderCountNumPairs();
   }
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
      localCountNumPairs();
   }
}

//count the number of pairs for the joint matrix. 

void
JointMat::sanderCountNumPairs()
{

   int num_replicates=mData->getNumRows();
   int size=mData->getNumCols();  


   //for each replications, compute the sum 
   for(int i=0; i<num_replicates; i++)
   {
      //for cols pos1 and pos2
      for(int pos1=0; pos1<size; pos1++)
      {
         int index1    = mData->getIndex(pos1, i);
         int off1      = mData->offset(pos1);
         double value1 = mData->getValue(pos1, i);
         
         mMargin[pos1][index1]+=value1*1/mReplicateWeight[i];
         
         for(int pos2=0; pos2<size; pos2++)
         {
            int off2      = mData->offset(pos2);
            int index2    = mData->getIndex(pos2, i);
            double value2 = mData->getValue(pos2, i);
            //cout<<off1<<" "<<index1<<" "<<off2<<" "<<index2<<endl;
            mJoint[off1+index1][off2+index2]+=value1*value2/mReplicateWeight[i]; //the replication has a weight.
            //cout<<value1*value2<<endl;
	 }//for
      }//for
   }

   //setup the defaults weights
   double total=0;
   for(int i=0; i<num_replicates; i++) total+=1/mReplicateWeight[i];

   for(int i=0; i<size; i++)
   {
      mSingularWeights[i]=total; //the total weight
      for(int j=0; j<size; j++)
        mWeights[i][j]=total; //the weeight for each pair
   }
   //exit(0);
}


void
JointMat::localCountNumPairs()
{
   int num_replicates=mData->getNumRows();
   int size=mData->getNumCols();
   
   
   double* pairs   =new double[size*size];
   double* singular=new double[size];


   for(int i=0; i<num_replicates; i++)
   {
      pairAndSingularWeights(i, pairs, singular);
      
      for(int pos1=0; pos1<size; pos1++)
      {
         int index1            = mData->getIndex(pos1, i);
         int off1              = mData->offset(pos1);
         double value1         = mData->getValue(pos1, i);
         mMargin[pos1][index1]+= value1*singular[pos1];
         
         for(int pos2=0; pos2<size; pos2++)
         {
            int off2      = mData->offset(pos2);
            double value2 = mData->getValue(pos2, i);

            int index2    = mData->getIndex(pos2, i);
            mJoint[off1+index1][off2+index2]+=value1*value2*pairs[pos1*size+pos2];
         }
      }//for
   }
   //exit(0);
   delete pairs;
   delete singular;
}

void
JointMat::pairAndSingularWeights(int row, double* weights, double* singularWeights)
{
   int size=mData->getNumCols();  
   for(int j=0; j<size*size; j++) weights[j]=0;
   for(int j=0; j<size;      j++) singularWeights[j]=0;

   int num_replicates=mData->getNumRows(); 

   for(int j=0; j<num_replicates; j++)
   {
      double sim=mData->similarity(row, j);
      
      for(int pos1=0; pos1<size; pos1++)
      {
          int index_a1 =    mData->getIndex(pos1, row );
          int index_b1 =    mData->getIndex(pos1, j   );

          singularWeights[pos1] += mData->localSingularSimilarity(pos1, index_a1, index_b1, sim);

          for(int pos2=0; pos2<size; pos2++)
          {
            int index_a2 = mData->getIndex(pos2, row);
            int index_b2 = mData->getIndex(pos2, j);

            weights[pos1*size+pos2] += mData->localPairSimilarity(pos1, pos2, index_a1, index_a2, 
                                                           index_b1, index_b2, sim);
          }
       }
   }

   for(int i=0; i<size; i++)
   {
      for(int j=0; j<size; j++)
      {
          weights[i*size+j]=1/weights[i*size+j];
          mWeights[i][j]+=weights[i*size+j];
          if(weights[i*size+j]==0)
          {
             cerr<<"Bugs 1 "<<endl;
          }
      }
   }
   for(int i=0; i<size; i++)
   {
      singularWeights[i]=1/singularWeights[i]; 
      mSingularWeights[i]+=singularWeights[i]; 

      if(singularWeights[i]==0)
      {
         cerr<<"Bugs 2 "<<endl;
      }
   }

   return;
}

//
//
//
//end of computing number of pairs
//======================



//begin convert counts to joint matrix
void
JointMat::computeJointMean()
{

   int num_cols = mNumCol;
 
   for(int i=0; i<num_cols; i++)
   {
      for(int j=0; j<num_cols; j++)
      {
         computeSubMatrixMean(i, j);
      }
   }
}


//I think you need a hidden conversion 
void
JointMat::computeSubMatrixMean(int pos1, int pos2)
{
   int off1=mOffset[pos1];
   int off2=mOffset[pos2];

   int dim1=mDim[pos1];
   int dim2=mDim[pos2];
   
   double sum=0;
   for(int h1=0; h1<dim1; h1++)
   {
      for(int h2=0; h2<dim2; h2++)
      {
	 mJoint[off1+h1][off2+h2] =joint(pos1, pos2, h1, h2);
      }//for
   }//for

   //for(int h1=0; h1<dim1; h1++)
   //{
   //   for(int h2=0; h2<dim2; h2++)
   //   {
	//sum+=joint(pos1, pos2, h1, h2);      
     //   mJoint[off1+h1][off2+h2]/=sum; 
     // }//for
   //}//for
}

double
JointMat::joint(int col_id1, int col_id2, int hidden1, int hidden2)
{
   //return mJoint[mOffset[col_id1]+hidden1][mOffset[col_id2]+hidden2];	

   	
   //if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   //{
      int off1 =   mData->offset(col_id1);
      int off2 =   mData->offset(col_id2);
      return mJoint[off1+hidden1][off2+hidden2]/mWeights[col_id1][col_id2];
   /*}
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
      double rev=0;
      int num_obs1=mData->dim(col_id1);
      int num_obs2=mData->dim(col_id2);

      int off1 =   mData->offset(col_id1);
      int off2 =   mData->getObserveOffset(col_id2);

      for(int obs1=0; obs1<num_obs1; obs1++)
      {
        for(int obs2=0; obs2<num_obs2; obs2++)
        {
            rev+=mJointSum[off1+obs1][off2+obs2]/mWeights[col_id1][col_id2]*mData->lookupK22(col_id1, col_id2, hidden1, hidden2, obs1, obs2);
        }
    }
    return rev;
   }
   else
   {
      cerr<<"matrix mode is unspecified!"<<endl;
   }*/
}
//
//
//
//
//end of converting of counts to joint matrix


void
JointMat::computeColMeans()
{
   int size = numCols();
   for(int i=0; i<size; i++)
   {
      int dim=mDim[i];
      double sum=0;
      for(int j=0; j<dim; j++)
      {
        mMargin[i][j]=margin(i, j);
      }

      //for(int j=0; j<dim; j++)
      //{
      //  mMargin[i][j]/=sum;//=margin()(i, j);
     // }
   }
}

double
JointMat::margin(int col_id, int hidden_clu)
{
   //if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   //{
      return mMargin[col_id][hidden_clu]/mSingularWeights[col_id];
   //}
   /*else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
     double rev=0;  
     int num_obs=mData->getObserveDim(col_id);
     
     for(int i=0; i<num_obs; i++)
     {
       rev+=mColSum[col_id][i]/mSingularWeights[col_id]*mData->lookupK1(col_id, hidden_clu, i);
     }
     return rev;
   }*/
}

void
JointMat::allocateSpace()
{

   int size     = mNumCol;
   int tot_size = mTotSize;
  
   //cout<<size<<" "<<tot_size<<endl;
   //exit(0);
   mJoint =new double *[tot_size]; //to score the joint frequencies
   for(int i=0; i< tot_size; i++)
   {
      mJoint[i]=new double[tot_size];
      for(int j=0; j<tot_size; j++)
      {
         mJoint[i][j]=0;
      }
   }

   mMargin=new double *[size]; //to score the single column frequencies
   for(int i=0; i<size; i++)
   {
     mMargin[i]=new double[dim(i)];
     for(int j=0; j<dim(i); j++)
     {
         mMargin[i][j]=0;
     }
   }

 /*  int tot_obs_size=mData->getNumObserveClus();

   mJointSum = new double* [tot_obs_size];
   for(int i=0; i<tot_obs_size; i++)
   {
     mJointSum[i]=new double[tot_obs_size];
     for(int j=0; j< tot_obs_size; j++)
     {
       mJointSum[i][j]=0;
     }
   }

   
   mColSum    =new double * [size];
   for(int i=0; i<size; i++)
   {
      mColSum[i]    =new double [mData->getObserveDim(i)];
      for(int j=0; j<mData->getObserveDim(i); j++)
      {
        mColSum[i][j]=0;
      }
   }
*/
   mWeights    =new double *[size]; //number of pairs for each (i, j)
   for(int i=0; i<size; i++)
   {
      mWeights[i]   =new double[size];
      for(int j=0; j<size; j++)
      {
         mWeights[i][j]   =0;
      }
   }
   
   mSingularWeights=new double [size];
   for(int i=0; i<size; i++)
   {
      mSingularWeights[i]=0;
   }
}

void
JointMat::computeReplicateWeight()
{
   int num_replicates=mData->getNumRows();
   mReplicateWeight =new double[num_replicates];
   
   for(int i=0; i<num_replicates; i++)
   {
      double weight=0;//weight for each replication

      for(int j=0;  j<num_replicates; j++)
      {
         double v=mData->similarity(i, j);
         if(v>0.7)
         weight+=1;//mData->similarity(i, j);
      }//for
      mReplicateWeight[i]=weight;
   }//for

}




/*
void
JointMat::computeColMeans()
{
   int size = mData->getNumCols();
   for(int i=0; i<size; i++)
   {
      int dim=mData->dim(i);
      for(int j=0; j<dim; j++)
      {
        mMargin[i][j]=margin()(i, j);
      }
   }
   
   for(int i=0; i<size; i++)
   { 
       int dim    = mData->dim(i);
       int off    = mData->offset(i);
       int obs_dim= mData->getObserveDim(i);
       for(int a=0; a<dim; a++)
       {
          for(int b=0; b<dim; b++)
          {
              mJoint[off+a][off+b]=0;
          }
       }

       for(int a=0; a<dim; a++)
       {
           for(int b=0; b<obs_dim; b++)
           {
             double v=mColSum[i][b]*mData->lookupK1(i, a, b);
             mJoint[off+a][off+a]+=mColSum[i][b]*mData->lookupK1(i, a, b);
           }
           mJoint[off+a][off+a]/=mSingularWeights[i];
       }
    }
}

*/
