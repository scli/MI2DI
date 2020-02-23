#include "LatentMap.h"
#include "MICParaMgr.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

LatentMap::LatentMap()
{
   mSource=NULL; 
   computeTransition();
   
   initData();
   initK22Table();
   initK1Table();


}


//begin convert counts to joint matrix

JointMat*
LatentMap::convert(JointMat* source)
{

   mSource=source;
   //allocateSpace();
   setData();//initilize the matrix
   for(int i=0; i<mTarget->numCols(); i++)
   {
      for(int j=0; j<mTarget->numCols(); j++)
      {
         computeSubMatrixMean(i, j);
      }
   }
   computeColMeans();
   return mTarget;
}



void
LatentMap::allocateSpace()
{

   //int size     = mData->getNumCols();

   //int tot_size = mData->getNumHiddenClus();
  
   //mTarget=new JointMap(tot_size, size, NULL, NULL);
   
   
   //exit(0);
   /*mTargetJointMean =new double *[tot_size]; //to score the joint frequencies
   for(int i=0; i< tot_size; i++)
   {
      mTargetJointMean[i]=new double[tot_size];
      for(int j=0; j<tot_size; j++)
      {
         mTargetJointMean[i][j]=0;
      }
   }

   mTargetColMean=new double *[size]; //to score the single column frequencies
   for(int i=0; i<size; i++)
   {
     mTargetColMean[i]=new double[mData->dim(i)];
     for(int j=0; j<mData->dim(i); j++)
     {
         mTargetColMean[i][j]=0;
     }
   }*/
}



//I think you need a hidden conversion 
void
LatentMap::computeSubMatrixMean(int pos1, int pos2)
{
   int off1=mTarget->offset(pos1);
   int off2=mTarget->offset(pos2);

   int dim1=mTarget->dim(pos1);
   int dim2=mTarget->dim(pos2);
   
   for(int h1=0; h1<dim1; h1++)
   {
      for(int h2=0; h2<dim2; h2++)
      {
        mTarget->set(pos1, h1, pos2, h2, getJointMean(pos1, pos2, h1, h2));
      }
   }
}


double
LatentMap::getJointMean(int col_id1, int col_id2, int hidden1, int hidden2)
{
   if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   {
      //int off1 =   mTarget->offset(col_id1);
      //int off2 =   mTarget->offset(col_id2);
      return mSource->get(col_id1, hidden1, col_id2, hidden2);
   }
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
      double rev=0;
      int num_obs1=mSource->dim(col_id1);
      int num_obs2=mSource->dim(col_id2);

      int off1 =   mSource->offset(col_id1);
      int off2 =   mSource->offset(col_id2);

      for(int obs1=0; obs1<num_obs1; obs1++)
      {
        for(int obs2=0; obs2<num_obs2; obs2++)
        {
            rev+=mSource->get(col_id1, obs1, col_id2, obs2)*lookupK22(col_id1, col_id2, hidden1, hidden2, obs1, obs2);
        }
    }
    return rev;
   }
   else
   {
      cerr<<"matrix mode is unspecified!"<<endl;
   }
   return 0.0;

}
//
//
//
//
//end of converting of counts to joint matrix


void
LatentMap::computeColMeans()
{
   int size = mTarget->numCols();
   for(int i=0; i<size; i++)
   {
      int dim=mTarget->dim(i);
      for(int j=0; j<dim; j++)
      {
        mTarget->set(i, j, getColMean(i, j));
      }
   }
}

double
LatentMap::getColMean(int col_id, int hidden_clu)
{
   if(MICParaMgr::_MATRIX_MODE_==_SANDER_)
   {
      return mSource->get(col_id, hidden_clu); 
   }
   else if(MICParaMgr::_MATRIX_MODE_==_LOCAL_WEIGHT_) 
   {
     double rev=0;  
     int num_obs=mSource->dim(col_id);
     
     for(int i=0; i<num_obs; i++)
     {
       rev+=mSource->get(col_id, i)*lookupK1(col_id, hidden_clu, i);
     }
     return rev;
   }
   return 0.0;
}


//initialize raw data

void
LatentMap::initK22Table()
{

  if(MICParaMgr::_DATA_OPT_==_THE_21_)
  {
    mK22Table=new double*[21*21];
    for(int i=0; i<21*21; i++) mK22Table[i]=new double[21*21];
    int h_dim=21;
    int o_dim=21;

     for(int h1=0; h1 < h_dim; h1++)
     {
       for(int h2=0; h2 < h_dim; h2++)
       {
        for(int obs1=0; obs1<o_dim; obs1++)
        {
           for(int obs2=0; obs2<o_dim; obs2++)
           {
             mK22Table[h_dim*h1+h2][o_dim*obs1+obs2]=K22(h1, h2, obs1, obs2);
           }
        }
       }//for
     }//for
   }
  else
  {
     cerr<<"Please specify a correct data option"<<endl;
  }
}

void
LatentMap::initK1Table()
{
  if(MICParaMgr::_DATA_OPT_==_THE_21_)
  {
   mK1Table=new double*[21];
   for(int i=0; i<21; i++) mK1Table[i]=new double[21];
 
    int h_dim=21;
    int o_dim=21;
    
    for(int h=0; h<h_dim; h++)
    {
        for(int obs=0; obs<o_dim; obs++)
        {
          mK1Table[h][obs]=K1(h,obs);
        }
    }
  }
  else
  {
     cerr<<"Please specify a correct data option"<<endl;
  }

}




void
LatentMap::initData()
{
   if(MICParaMgr::_DATA_OPT_==_THE_21_)
   {
      init21Data();
   }
}

void
LatentMap::init21Data()
{
}


    
double 
LatentMap::K22(int aa_i, int aa_j, int bb_i, int bb_j)
{
   double f1;
   f1=mTransition[aa_i][bb_i]/mTransMargin[bb_i];
   double f2;
   f2=mTransition[aa_j][bb_j]/mTransMargin[bb_j];
   return f1*f2;
;
}

double 
LatentMap::K1(int a, int a_prime)
{
   double f1;
   f1=mTransition[a][a_prime]/mTransMargin[a_prime];
   return f1;  
    ;
}




double
LatentMap::lookupK22(int pos1, int pos2, int aa1, int aa2, int bb1, int bb2)
{
   return  mK22Table[21*aa1+aa2][21*bb1+bb2];
}

double
LatentMap::lookupK1(int pos, int aa1,int bb1)
{
   return  mK1Table[aa1][bb1];
}



void
LatentMap::allocateTransition()
{
   if(MICParaMgr::_DATA_OPT_==_THE_21_)
   {
      mTransition=new double*[21];
      for(int i=0; i<21; i++)
      {
         mTransition[i]=new double[21];
         for(int j=0; j<21; j++)
         {
             mTransition[i][j]=20;
         } 
      }//
   }

}


void
LatentMap::computeTransition()
{
   allocateTransition();
  
   if(MICParaMgr::_DATA_OPT_==_THE_21_)
   {
      for(int j1=0; j1<=20; j1++)
      {
        for(int j2=0; j2<=20; j2++)
       {
         if(j1==20 && j2==20)
         {
              mTransition[j1][j2]=1;//0.000000001;//0.000000;beta;//mLocalBlosumJoint[j1][j2];
         }
         else if(j1==20)
         {
             mTransition[j1][j2]=0.00001;//1.0/1240000;//aaprob90[j2]/42;
             //mLocalBlosumJoint[j1][j2];//0.000;//aaprob90[j2]/42;//mLocalBlosumJoint[j1][j2];
         }
         else if(j2==20)
         {
             mTransition[j1][j2]=0.00001;//0.05;//1.0/1240000;//aaprob90[j2]/42;
             //mLocalBlosumJoint[j1][j2];//0.000;//aaprob90[j1]/42;
         }
         else
         {
           mTransition[j1][j2]=blosum100joint[j1][j2];//Q0[j1][j2];//0*mLocalBlosumJoint[j1][j2]+1*blosum90joint[j1][j2];

         }

         if(j1==j2) mTransition[j1][j2]=1; else mTransition[j1][j2]=0;
       }
      }
      normalize(mTransition, 21);
      mTransMargin=new double[21];

      for(int col=0; col<21; col++)
      {
        double sum=0;
        for(int r=0; r<21; r++)
        {
           sum+=mTransition[r][col];
           //if(j1!=j2) mLocalBlosumJoint[j2][j1]+=counters[j2]*counters[j1];
        }
        mTransMargin[col]=sum;
      }
   }
}




//=======================================
//
//
//
//=======================================


void
LatentMap::normalize(double** vect, int size)
{
   double sum=0;
   for(int i=0; i<size; i++)
   {
       for(int j=0; j<size; j++)
       {
          sum+=vect[i][j];
       }
   }

   for(int i=0; i<size; i++)
   {
       for(int j=0; j<size; j++)
       {
          vect[i][j]/=sum;
       }
   }
}


void
LatentMap::setData()
{

    
   if(MICParaMgr::_DATA_OPT_==_THE_21_)
   {

     
     //mData->setData(mRawData);
     //mData->setK22Table(mK22Table);
     //mData->initK1Table(mK1Table);
     int numCols=mSource->numCols();
     int*    dim=new int[numCols];
     int* offset=new int[numCols];     
     for(int i=0; i<numCols; i++)
     {
        dim[i]    = 21;
        offset[i] = 21*i;	
        //mData->setObserveDim(i, 21);
        //mData->setHiddenDim(i, 21);
        //mData->setObserveOffset(i, 21*i);
        //mData->setHiddenOffset(i,  21*i);
     }
     //mData->setNumHiddenClusters( 21*mNumCols);
     //mData->setNumObserveClusters(21*mNumCols);
     //cout<<mNumCols<<" "<<21*mNumCols<<endl;
     mTarget=new JointMat(21*numCols, numCols, offset, dim);
   }
   
}
