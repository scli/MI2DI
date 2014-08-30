#include "JointMat.h"
#include <iostream>
#include <stdlib.h>
#include "MICParaMgr.h"
using namespace std;


JointMat::JointMat(MICData* data)
{
   mData=data; 

   allocateSpace();
   computeReplicateWeight();
   countNumPairs();
   computeJointMean(); 
   computeColMeans();
}

JointMat::~JointMat()
{
  ;
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
      sanderCountNumPairs();
   }
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
      localCountNumPairs();
   }
}


void
JointMat::sanderCountNumPairs()
{

   int num_replicates=mData->getNumRows();
   int size=mData->getNumCols();  
   //next, count the number of Cia, jb
   for(int i=0; i<num_replicates; i++)
   {
      for(int pos1=0; pos1<size; pos1++)
      {
         int index1=mData->getIndex(pos1, i);
         int off1=mData->getObserveOffset(pos1);
         double value1 = mData->getValue(pos1, i);
         
         mColSum[pos1][index1]+=value1*1/mReplicateWeight[i];
         
         for(int pos2=0; pos2<size; pos2++)
         {
            int off2=mData->getObserveOffset(pos2);
            int index2=mData->getIndex(pos2, i);
            double value2 = mData->getValue(pos2, i);
            //cout<<off1<<" "<<index1<<" "<<off2<<" "<<index2<<endl;
            mJointSum[off1+index1][off2+index2]+=value1*value2/mReplicateWeight[i];
         }//for
      }//for
     // cout<<i<<" "<<(1/pairs[0])<<" "<<(1/singular[0])<<endl;
     //cout<<i<<" "<<mReplicateWeight[i]<<endl;
   }

   //setup the defaults weights
   double total=0;
   for(int i=0; i<num_replicates; i++) total+=1/mReplicateWeight[i];

   for(int i=0; i<size; i++)
   {
      mSingularWeights[i]=total;
      for(int j=0; j<size; j++)
        mWeights[i][j]=total;
   }
   //exit(0);
}


void
JointMat::localCountNumPairs()
{
   int num_replicates=mData->getNumRows();
   int size=mData->getNumCols();
   double* pairs=new double[size*size];
   double* singular=new double[size];


   for(int i=0; i<num_replicates; i++)
   {
      pairAndSingularWeights(i, pairs, singular);
      
      for(int pos1=0; pos1<size; pos1++)
      {
         int index1=mData->getIndex(pos1, i);
         int off1=mData->getObserveOffset(pos1);
         double value1 = mData->getValue(pos1, i);
         mColSum[pos1][index1]+=value1*singular[pos1];
         
         for(int pos2=0; pos2<size; pos2++)
         {
            int off2=mData->getObserveOffset(pos2);
            double value2 = mData->getValue(pos2, i);

            int index2=mData->getIndex(pos2, i);
            mJointSum[off1+index1][off2+index2]+=value1*value2*pairs[pos1*size+pos2];
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

   int num_cols = mData->getNumCols();
 
   for(int i=0; i<num_cols; i++)
   {
      for(int j=0; j<num_cols; j++)
      {
         computeSubMatrixMean(i, j);
      }
   }
}

void
JointMat::computeSubMatrixMean(int pos1, int pos2)
{
   int off1=mData->getHiddenOffset(pos1);
   int off2=mData->getHiddenOffset(pos2);

   int dim1=mData->getHiddenDim(pos1);
   int dim2=mData->getHiddenDim(pos2);
   
   for(int h1=0; h1<dim1; h1++)
   {
      for(int h2=0; h2<dim2; h2++)
      {
        mJointMean[off1+h1][off2+h2]=getJointMean(pos1, pos2, h1, h2);
      }
   }
}


double
JointMat::getJointMean(int col_id1, int col_id2, int hidden1, int hidden2)
{
   if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   {
      int off1 =   mData->getObserveOffset(col_id1);
      int off2 =   mData->getObserveOffset(col_id2);
      return mJointSum[off1+hidden1][off2+hidden2]/mWeights[col_id1][col_id2];
   }
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
      double rev=0;
      int num_obs1=mData->getObserveDim(col_id1);
      int num_obs2=mData->getObserveDim(col_id2);

      int off1 =   mData->getObserveOffset(col_id1);
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
   }

}
//
//
//
//
//end of converting of counts to joint matrix


void
JointMat::computeColMeans()
{
   int size = mData->getNumCols();
   for(int i=0; i<size; i++)
   {
      int dim=mData->getHiddenDim(i);
      for(int j=0; j<dim; j++)
      {
        mColMean[i][j]=getColMean(i, j);
      }
   }
}

double
JointMat::getColMean(int col_id, int hidden_clu)
{
   if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   {
      return mColSum[col_id][hidden_clu]/mSingularWeights[col_id];
   }
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
     double rev=0;  
     int num_obs=mData->getObserveDim(col_id);
     
     for(int i=0; i<num_obs; i++)
     {
       rev+=mColSum[col_id][i]/mSingularWeights[col_id]*mData->lookupK1(col_id, hidden_clu, i);
     }
     return rev;
   }
}

void
JointMat::allocateSpace()
{

   int size     = mData->getNumCols();
   int tot_size = mData->getNumHiddenClus();
  
   //cout<<size<<" "<<tot_size<<endl;
   //exit(0);
   mJointMean =new double *[tot_size]; //to score the joint frequencies
   for(int i=0; i< tot_size; i++)
   {
      mJointMean[i]=new double[tot_size];
      for(int j=0; j<tot_size; j++)
      {
         mJointMean[i][j]=0;
      }
   }

   mColMean=new double *[size]; //to score the single column frequencies
   for(int i=0; i<size; i++)
   {
     mColMean[i]=new double[mData->getHiddenDim(i)];
     for(int j=0; j<mData->getHiddenDim(i); j++)
     {
         mColMean[i][j]=0;
     }
   }

   int tot_obs_size=mData->getNumObserveClus();

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
      int dim=mData->getHiddenDim(i);
      for(int j=0; j<dim; j++)
      {
        mColMean[i][j]=getColMean(i, j);
      }
   }
   
   for(int i=0; i<size; i++)
   { 
       int dim    = mData->getHiddenDim(i);
       int off    = mData->getHiddenOffset(i);
       int obs_dim= mData->getObserveDim(i);
       for(int a=0; a<dim; a++)
       {
          for(int b=0; b<dim; b++)
          {
              mJointMean[off+a][off+b]=0;
          }
       }

       for(int a=0; a<dim; a++)
       {
           for(int b=0; b<obs_dim; b++)
           {
             double v=mColSum[i][b]*mData->lookupK1(i, a, b);
             mJointMean[off+a][off+a]+=mColSum[i][b]*mData->lookupK1(i, a, b);
           }
           mJointMean[off+a][off+a]/=mSingularWeights[i];
       }
    }
}

*/
