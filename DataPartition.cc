#include "DataPartition.h"
#include "MICParaMgr.h"
#include <iostream>
#include <stdlib.h>
using namespace std;
int compare (const void * a, const void * b)
{
  if( *(double*)a - *(double*)b <0) return -1;
  return 1;
}

DataPartition::DataPartition(double** rawData, int n_row, int n_col)
{
   mRawData=rawData;
   mNumRows=n_row;
   mNumCols=n_col;
   //cout<<n_row<<" "<<n_col<<endl;

}

DataPartition::~DataPartition()
{

}

void
DataPartition::partition()
{

    if(MICParaMgr::_DATA_PARTITION_OPT_==_UNIFORM_PARTITION_)
    {
        uniformPartition();
    }
    else
    {
    
       cerr<<"Data partition option is unspecified "<<endl;
       exit(0);
    }
    setDiscrData();
    setCntnsData();
}


void
DataPartition::uniformPartition()
{
  
  int numOfClusters=MICParaMgr::_NUM_PARTITION_CLUSTERS_;

  mNumClusters=new int[mNumCols];
  for(int i=0; i<mNumCols; i++)
  {
    mNumClusters[i]=numOfClusters;
  }



  mBounds=new double *[mNumCols];
  //compute dissect bounds; 
  //in uniform, each 

  double * col=new double[mNumRows];
  
  for(int i=0; i<mNumCols; i++)
  {
     //cout<<i<<" "<<numOfClusters<<endl;
     mBounds[i]=new double[numOfClusters];
     
     
     for(int j=0; j<mNumRows; j++)
         col[j] = mRawData[j][i];

     qsort(col, mNumRows, sizeof(double), compare); 
     for(int j=0; j<numOfClusters; j++)     
     {
       mBounds[i][j]=col[j*mNumRows/numOfClusters];
       //cout<<mBounds[i][j]<<" ";
     }
     //cout<<endl;
  }
  //exit(0);
  mCluData=new DisDataT*[mNumRows];
  for(int i=0; i<mNumRows; i++)
  {
     mCluData[i] =new DisDataT[mNumCols];
     
     for(int j=0; j<mNumCols; j++)
     {
       int which=0;
       double v=mRawData[i][j];
       int c_id=mNumClusters[j]-1;
       while(v<mBounds[j][c_id])
       {
          c_id--; 
       }
       if(c_id<0) 
       {
         cout<<v<<" "<<mBounds[j][0]<<endl;  
         cerr<<"please debug "<<c_id<<endl;
         exit(0);
       }
       mCluData[i][j]=c_id;
     }
  }

  mK22Table=new double*[numOfClusters*numOfClusters];
  for(int i=0; i<numOfClusters*numOfClusters; i++)
  {
     mK22Table[i] = new double [numOfClusters*numOfClusters];
     for(int j=0; j<numOfClusters*numOfClusters; j++)
     {
        mK22Table[i][j]=0;
        if(i%numOfClusters==j%numOfClusters && i/numOfClusters==j/numOfClusters)
            mK22Table[i][j]=1;
     
     }//for
  }
  mK1Table =new double* [numOfClusters];
  for(int i=0; i<numOfClusters; i++)
  {
     mK1Table[i] = new double [numOfClusters];
     for(int j=0; j<numOfClusters; j++)
     {
        mK1Table[i][j]=0;
        if(i==j)
            mK1Table[i][j]=1;
     
     }//for
  }
  
  mMaxNumClusters=numOfClusters;

  delete [] col;
}

void
DataPartition::setDiscrData()
{

   
   mDiscrData=new DiscrData(mNumRows, mNumCols, mMaxNumClusters);
   mDiscrData->setData(mCluData);
   mDiscrData->setK22Table(mK22Table);
   mDiscrData->initK1Table(mK1Table);
   
   int offset1=0;
   int offset2=0;
   for(int i=0; i<mNumCols; i++)
   {
      mDiscrData->setObserveDim(i, mNumClusters[i]);
      mDiscrData->setHiddenDim(i,  mNumClusters[i]);
      mDiscrData->setObserveOffset(i, offset1);
      mDiscrData->setHiddenOffset(i,  offset2);
      offset1+=mNumClusters[i];
      offset2+=mNumClusters[i];
   }
   mDiscrData->setNumHiddenClusters(offset1);
   mDiscrData->setNumObserveClusters(offset2);
   //cout<<offset1<<" "<<offset2<<endl;
   //cout<<mNumCols<<" "<<21*mNumCols<<endl;
}
void
DataPartition::setCntnsData()
{
   mCntnsData=new CntnsData(mNumRows, mNumCols);
   mCntnsData->setRawData(mRawData);
   mCntnsData->setIndex(mCluData);
   //mCntnsData->setK22Table(mK22Table);
   //mCntnsData->initK1Table(mK1Table);


   int offset1=0;
   int offset2=0;
   for(int i=0; i<mNumCols; i++)
   {
      mCntnsData->setObserveDim(i, mNumClusters[i]);
      mCntnsData->setHiddenDim(i,  mNumClusters[i]);
      mCntnsData->setObserveOffset(i, offset1);
      mCntnsData->setHiddenOffset(i,  offset2);
      offset1+=mNumClusters[i];
      offset2+=mNumClusters[i];
   }
   mCntnsData->setNumHiddenClusters(offset1);
   mCntnsData->setNumObserveClusters(offset2);
}

