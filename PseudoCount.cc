#include "PseudoCount.h"
#include <iostream>
#include "MICParaMgr.h"
using namespace std;
PseudoCount::PseudoCount(JointMat* data, LatentMap* map, DataType t)
{
  mData       = data;
  mMap        = map;
  mJointMatrix= data->joint();
  mMargin     = data->margin();
  mDataType   = t; 
  allocateSpace();
  psuedoCount();
  mPSMat      =new JointMat(data->numClus(),
		            data->numCols(),
			    data->offset(),
			    data->dim(),
			    mPseudoMargin, 
			    mPseudoJointMatrix);

}





PseudoCount::~PseudoCount()
{

}
void
PseudoCount::psuedoCount()
{
   
   DataOpt opt=MICParaMgr::_DATA_OPT_;
   
   //if(opt==_THE_21_) 
   if(mDataType==_DISCRETE_)
   {
      uniformPseudoCount();
   }
   else if(mDataType==_CONTINUOUS_)
   {
     uniformContinuousPseudoCount();
   }
   else
   {
       cerr<<"Wrong pseudo count data type"<<endl;
   }
}


void
PseudoCount::allocateSpace()
{
   int tot_dim=mData->numClus();

   mPseudoJointMatrix=new double *[tot_dim];
   for(int i=0; i<tot_dim; i++)
   {
        mPseudoJointMatrix[i]=new double[tot_dim];
   }

   int  num_vars=mData->numCols();

   mPseudoMargin    =new double* [num_vars];
   for(int i=0; i<num_vars; i++)
   {
        mPseudoMargin[i]=new double[mData->dim(i)];
   }
}

void
PseudoCount::uniformPseudoCount()
{
    int size=mData->numClus();
    double pseudoCount=MICParaMgr::_PSEUDO_COUNT_WEIGHT_;

    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)
        {
          mPseudoJointMatrix[i][j]=(1.0-pseudoCount)*mJointMatrix[i][j];//+pseudoCount/ALPHABET_SIZE.0/ALPHABET_SIZE.0;
        }
    }
    JointMat* source=mMap->source();
    int cols=mData->numCols();
    for(int i=0; i<cols; i++)
    {
       int h_dim1= mData->dim(i);
       int o_dim1=source->dim(i);
       for(int j=0; j<cols; j++)
       {
           if(i==j)  continue;
           int h_dim2=mData->dim(j);
           int o_dim2=source->dim(j);

           for(int aa1=0; aa1<o_dim1; aa1++)
           {
               for(int aa2=0; aa2<o_dim2; aa2++)
               {
                   for(int bb1=0; bb1<h_dim1; bb1++)
                   {
                       for(int bb2=0; bb2<h_dim2; bb2++)
                       {
                          mPseudoJointMatrix[mData->offset(i)+bb1][mData->offset(j)+bb2]
                              +=pseudoCount/o_dim1/o_dim2*mMap->lookupK22(i, j, bb1, bb2, aa1, aa2);
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
       int h_dim=mData->dim(i);
       int o_dim=source->dim(i);
       
       int h_off=mData->offset(i);
       for(int a1=0; a1<h_dim; a1++)
       {
           for(int b=0; b<o_dim; b++)
           {
               mPseudoJointMatrix[h_off+a1][h_off+a1]+=pseudoCount/h_dim*mMap->lookupK1(i, a1, b);
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
        int h_dim=mData->dim(i);
        int h_off=mData->offset(i);
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
    }//for
}
//


void
PseudoCount::uniformContinuousPseudoCount()
{
    int size=mData->numClus();
    double pseudoCount=0;MICParaMgr::_PSEUDO_COUNT_WEIGHT_;

    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)
        {
          mPseudoJointMatrix[i][j]=(1.0-pseudoCount)*mJointMatrix[i][j];//+pseudoCount/ALPHABET_SIZE.0/ALPHABET_SIZE.0;
        }
    }
    int cols=mData->numCols();
    for(int i=0; i<cols; i++)
    {
       int h_dim1=mData->dim(i);
       for(int j=0; j<cols; j++)
       {
           if(i==j)  continue;
           int h_dim2=mData->dim(j);
           for(int bb1=0; bb1<h_dim1; bb1++)
           {
              for(int bb2=0; bb2<h_dim2; bb2++)
             {
                    mPseudoJointMatrix[mData->offset(i)+bb1][mData->offset(j)+bb2]
                         +=pseudoCount/h_dim1/h_dim2*mMargin[i][bb1]*mMargin[j][bb2];//mData->lookupK22(i, j, bb1, bb2, aa1, aa2);
             }
           }
       }
    }

    for(int i=0; i<cols; i++)
    {
       int h_dim=mData->dim(i);
       int h_off=mData->offset(i);
       for(int a=0; a<h_dim; a++)
       {
          mPseudoJointMatrix[h_off+a][h_off+a]+=pseudoCount/h_dim*mMargin[i][a]*mMargin[i][a];//
       }//for
    }//for
    
    for(int i=0; i<cols; i++)
    {
        int h_dim=mData->dim(i);

        for(int j=0; j<h_dim; j++)
        {
          mPseudoMargin[i][j]=mMargin[i][j];//+pseudoCount/ALPHABET_SIZE.0;
        }//f
    }//for*/
}
