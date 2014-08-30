#include "AlignData.h"
#include "MICParaMgr.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

AlignData::AlignData(AlignmentParser* parser) 
{
   mParser=parser;

   initDataBlock();
   computeTransition();
   
   initData();
   initK22Table();
   initK1Table();

   setData();

}



int
AlignData::getAAIndex(char C)
{
   if(MICParaMgr::_INPUT_TYPE_==_PROTEIN_)
   {
     if(C<'A' || C>'Z')
     {
       return 20; 
     }
     else
     {
        int rev=AA2Index[C-'A'];
          if(rev<0) return 20;
        return rev;
     }
   }
   else 
   {
       cerr<<"Please specify your input type, "<<endl;
   }
}

//initialize raw data
void
AlignData::initDataBlock()
{
   vector<string>* seqs=mParser->getAlignStr(); 
   
   mRawData=new DisDataT *[seqs->size()];
   for(int i=0; i<seqs->size(); i++)
   {
      string& s=(*seqs)[i];
      mRawData[i]=new DisDataT[s.size()];
      
      for(int j=0; j<s.size(); j++)
      {
         mRawData[i][j]=getAAIndex( s[j] ); 
      }//for
      //cout<<s.c_str()<<endl;
   }
}



void
AlignData::initK22Table()
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
AlignData::initK1Table()
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
AlignData::initData()
{
   vector<string>* seqs=mParser->getAlignStr(); 
   if(seqs->size()==0)
   {
       cerr<<"No sequences are available "<<endl;
       exit(0);
   }
   if(MICParaMgr::_DATA_OPT_==_THE_21_)
   {
      init21Data();
   }
}

void
AlignData::init21Data()
{
  vector<string>* seqs=mParser->getAlignStr(); 
  mNumRows=seqs->size();
  mNumCols=(*seqs)[0].size();
 /* 
  mObserveDim     = new int[mNumCols];
  mObserveOffset  = new int[mNumCols];

  mHiddenDim      = new int[mNumCols];
  mHiddenOffset   = new int[mNumCols];*/
  //mNumHiddenClus  = 21*mNumCols;
 // mNumObserveClus = 21*mNumCols;

/*
  for(int i=0; i<mNumCols; i++)
  {
    mObserveDim[i]   =21;
    mHiddenDim[i]    =21;
    mObserveOffset[i]=21*i;
    mHiddenOffset[i] =21*i;
  }*/
}


    
double 
AlignData::K22(int aa_i, int aa_j, int bb_i, int bb_j)
{
   double f1;
   f1=mTransition[aa_i][bb_i]/mTransMargin[bb_i];
   double f2;
   f2=mTransition[aa_j][bb_j]/mTransMargin[bb_j];
   return f1*f2;
;
}

double 
AlignData::K1(int a, int a_prime)
{
   double f1;
   f1=mTransition[a][a_prime]/mTransMargin[a_prime];
   return f1;  
    ;
}

///
//
//
//
// compute the transition matrix, from obs-->hidden
//
//

void
AlignData::allocateTransition()
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
AlignData::computeTransition()
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
AlignData::normalize(double** vect, int size)
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
AlignData::setData()
{

   mData=NULL;
    
   if(MICParaMgr::_DATA_OPT_==_THE_21_)
   {

     mData=new DiscrData(mNumRows, mNumCols, 21);
     
     mData->setData(mRawData);
     mData->setK22Table(mK22Table);
     mData->initK1Table(mK1Table);
     for(int i=0; i<mNumCols; i++)
     {
        mData->setObserveDim(i, 21);
        mData->setHiddenDim(i, 21);
        mData->setObserveOffset(i, 21*i);
        mData->setHiddenOffset(i,  21*i);
     }
     mData->setNumHiddenClusters( 21*mNumCols);
     mData->setNumObserveClusters(21*mNumCols);
     //cout<<mNumCols<<" "<<21*mNumCols<<endl;
   }
}
