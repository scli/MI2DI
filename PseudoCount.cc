#include "PseudoCount.h"
#include <iostream>
#include "MICParaMgr.h"
using namespace std;
PseudoCount::PseudoCount(MICData* data)
{
  mData=data;
  mJointMatrix=data->getJointMatrix();
  mMargin         =data->getMargin();
  allocateSpace();
  psuedoCount();
}

PseudoCount::~PseudoCount()
{

}
void
PseudoCount::psuedoCount()
{
   
   //PseudoOpt opt=ParaMgr::PSEUDO_COUNT_OPT;
   //if(opt==UNIFORM)
   if(mData->getDataType()==_DISCRETE_)
   {
      uniformPseudoCount();
   }
   else
   {
     uniformContinuousPseudoCount();
   }
}


void
PseudoCount::allocateSpace()
{
   int tot_dim=mData->getNumHiddenClus();

   mPseudoJointMatrix=new double *[tot_dim];
   for(int i=0; i<tot_dim; i++)
   {
        mPseudoJointMatrix[i]=new double[tot_dim];
   }

   int  num_vars=mData->getNumCols();

   mPseudoMargin    =new double* [num_vars];
   for(int i=0; i<num_vars; i++)
   {
        mPseudoMargin[i]=new double[mData->getHiddenDim(i)];
   }
}

void
PseudoCount::uniformPseudoCount()
{
    int size=mData->getNumHiddenClus();
    double pseudoCount=MICParaMgr::_PSEUDO_COUNT_WEIGHT_;

    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)
        {
          mPseudoJointMatrix[i][j]=(1.0-pseudoCount)*mJointMatrix[i][j];//+pseudoCount/ALPHABET_SIZE.0/ALPHABET_SIZE.0;
        }
    }
    int cols=mData->getNumCols();
    for(int i=0; i<cols; i++)
    {
       int h_dim1=mData->getHiddenDim(i);
       int o_dim1=mData->getObserveDim(i);
       for(int j=0; j<cols; j++)
       {
           if(i==j)  continue;
           int h_dim2=mData->getHiddenDim(j);
           int o_dim2=mData->getObserveDim(j);

           for(int aa1=0; aa1<o_dim1; aa1++)
           {
               for(int aa2=0; aa2<o_dim2; aa2++)
               {
                   for(int bb1=0; bb1<h_dim1; bb1++)
                   {
                       for(int bb2=0; bb2<h_dim2; bb2++)
                       {
                          mPseudoJointMatrix[mData->getHiddenOffset(i)+bb1][mData->getHiddenOffset(j)+bb2]
                              +=pseudoCount/o_dim1/o_dim2*mData->lookupK22(i, j, bb1, bb2, aa1, aa2);
                       }
                   }
               }
           }
       }
    }

/*
     for(int i=0; i<cols; i++)
     {
       for(int j=0; j<cols; j++)
       {
           if(i==j) continue;
           double sum=0;
           for(int aa1=0; aa1<ALPHABET_SIZE; aa1++)
           {
               for(int aa2=0; aa2<ALPHABET_SIZE; aa2++)
               {
                    sum+=mPseudoJointMatrix[i*ALPHABET_SIZE+aa1][j*ALPHABET_SIZE+aa2];
               }
           }
           //    cout<<sum<<endl;
           
           for(int aa1=0; aa1<ALPHABET_SIZE; aa1++)
           {
               for(int aa2=0; aa2<ALPHABET_SIZE; aa2++)
               {
                  mPseudoJointMatrix[i*ALPHABET_SIZE+aa1][j*ALPHABET_SIZE+aa2]/=sum;
               }
           }
       }
    }
*/

    for(int i=0; i<cols; i++)
    {
       int h_dim=mData->getHiddenDim(i);
       int o_dim=mData->getObserveDim(i);
       
       int h_off=mData->getHiddenOffset(i);
       for(int a1=0; a1<h_dim; a1++)
       {
           for(int b=0; b<o_dim; b++)
           {
               mPseudoJointMatrix[h_off+a1][h_off+a1]+=pseudoCount/h_dim*mData->lookupK1(i, a1, b);
           }
       }//for

        /*double sum=0;
        for(int aa1=0; aa1<ALPHABET_SIZE; aa1++)
        {
             for(int aa2=0; aa2<ALPHABET_SIZE; aa2++)
             {
                sum+=mPseudoJointMatrix[i*ALPHABET_SIZE+aa1][i*ALPHABET_SIZE+aa2];
             }
        }
        cout<<sum<<endl;*/
        /*
        for(int aa1=0; aa1<ALPHABET_SIZE; aa1++)
        {
             for(int aa2=0; aa2<ALPHABET_SIZE; aa2++)
             {
                  mPseudoJointMatrix[i*ALPHABET_SIZE+aa1][i*ALPHABET_SIZE+aa2]/=sum;
             }
        }
        */
    }//for
    
    for(int i=0; i<cols; i++)
    {
        int h_dim=mData->getHiddenDim(i);
        int h_off=mData->getObserveOffset(i);
        for(int j=0; j<h_dim; j++)
        {
          mPseudoMargin[i][j]=mPseudoJointMatrix[h_off+j][h_off+j];//+pseudoCount/ALPHABET_SIZE.0;
        }//f
    //}
     /*   for(int j=0; j<ALPHABET_SIZE; j++)
        {
          mPseudoMargin[i][j]=(1-pseudoCount)*mMargin[i][j];//+pseudoCount/ALPHABET_SIZE.0;
        }//for
    
        for(int aa=0; aa<ALPHABET_SIZE; aa++)
        {
           for(int bb=0; bb<ALPHABET_SIZE; bb++)
           {
              mPseudoMargin[i][aa]+=pseudoCount/ALPHABET_SIZE.0*lookupK1(aa, bb);
           }
        }*/
        double sum=0;
        for(int j=0; j<h_dim; j++)
        {
           sum+=mPseudoMargin[i][j];//=(1-pseudoCount)*mMargin[i][j];//+pseudoCount/ALPHABET_SIZE.0;
        }
        //cout<<sum<<" la la "<<endl;
        for(int j=0; j<h_dim; j++)
        {
           mPseudoMargin[i][j]/=sum;//=(1-pseudoCount)*mMargin[i][j];//+pseudoCount/ALPHABET_SIZE.0;
         //  cout<<mPseudoMargin[i][j]<<" ";
        }
       // cout<<endl;
    }//for*/
}
//


void
PseudoCount::uniformContinuousPseudoCount()
{
    int size=mData->getNumHiddenClus();
    double pseudoCount=0;MICParaMgr::_PSEUDO_COUNT_WEIGHT_;

    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)
        {
          mPseudoJointMatrix[i][j]=(1.0-pseudoCount)*mJointMatrix[i][j];//+pseudoCount/ALPHABET_SIZE.0/ALPHABET_SIZE.0;
        }
    }
    int cols=mData->getNumCols();
    for(int i=0; i<cols; i++)
    {
       int h_dim1=mData->getHiddenDim(i);
       for(int j=0; j<cols; j++)
       {
           if(i==j)  continue;
           int h_dim2=mData->getHiddenDim(j);
           for(int bb1=0; bb1<h_dim1; bb1++)
           {
              for(int bb2=0; bb2<h_dim2; bb2++)
             {
                    mPseudoJointMatrix[mData->getHiddenOffset(i)+bb1][mData->getHiddenOffset(j)+bb2]
                         +=pseudoCount/h_dim1/h_dim2*mMargin[i][bb1]*mMargin[j][bb2];//mData->lookupK22(i, j, bb1, bb2, aa1, aa2);
             }
           }
       }
    }

    for(int i=0; i<cols; i++)
    {
       int h_dim=mData->getHiddenDim(i);
       int h_off=mData->getHiddenOffset(i);
       for(int a=0; a<h_dim; a++)
       {
          mPseudoJointMatrix[h_off+a][h_off+a]+=pseudoCount/h_dim*mMargin[i][a]*mMargin[i][a];//
       }//for
    }//for
    
    for(int i=0; i<cols; i++)
    {
        int h_dim=mData->getHiddenDim(i);

        for(int j=0; j<h_dim; j++)
        {
          mPseudoMargin[i][j]=mMargin[i][j];//+pseudoCount/ALPHABET_SIZE.0;
        }//f
    }//for*/
}
