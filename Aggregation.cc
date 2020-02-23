#include "Aggregation.h"
#include "MICData.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

Aggregation::Aggregation(MICData* data)
{
  mData=data;
  mDIMatrix=new double* [data->getNumCols()];
  for(int i=0; i<data->getNumCols(); i++)
  {
     mDIMatrix[i]=new double[data->getNumCols()];
  }
}


Aggregation::~Aggregation()
{


}
    


//=================================
/*void
Aggregation::aggregateByMinEdge()
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
Aggregation::findMinimum(int v1, int v2)
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

/*
void
Aggregation::sanDiago()
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
       int dim1=mData->getHiddenDim(i);
       int off1=mData->getHiddenOffset(i);
       for (int j=0; j<num_col; j++) 
       {
          int dim2=mData->getHiddenDim(j);
          int off2=mData->getHiddenOffset(j);
          
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


//vector times a matrix of dimeneion dim1 X dim2, result stores the results
void
Aggregation::vectorTimesMatrix(double* vect, double** mat, int dim1, int dim2, double * result)
{
   for(int i=0; i<dim2; i++)
   {
      double value=0;
      for(int j=0; j<dim1; j++)
      {
         value+=vect[j]*mat[j][i];
      }
      result[i]=value;
   }//for
}


//vector vect times the transpose of the a matrix mat: result=vect X Mat^t, 
void
Aggregation::vectorTimesTranspose(double* vect, double** mat, int dim1, int dim2, double * result)
{
   for(int i=0; i<dim1; i++)
   {
      double value=0;
      for(int j=0; j<dim2; j++)
      {
         value+=vect[j]*mat[i][j];
      }
      result[i]=value;
   }//for
}

//such as the sum of the vect values will be 1, 
void
Aggregation::normalize(double* vect, int size)
{
   double rev=0;
   for(int i=0; i<size; i++)
   {
      rev+=vect[i];
   }

   for(int i=0; i<size; i++)
   {
      vect[i]/=rev;
   }
}

//such as the sume of all the elements in the mat will be 1
void
Aggregation::normalize(double** mat, int dim1, int dim2)
{
   double sum=0;
   for(int i=0; i<dim1; i++)
   {
       for(int j=0; j<dim2; j++)
       {
          sum+=mat[i][j];
       }
   }

   for(int i=0; i<dim1; i++)
   {
       for(int j=0; j<dim2; j++)
       {
          mat[i][j]/=sum;
       }
   }
}




// result res will store v1 dot 1/v2, element wise 
void
Aggregation::innerDivide(double* v1, double* v2, double* res, int size)
{
   for(int i=0; i<size; i++)
   {
      res[i]=v1[i]/v2[i];
   }
}


//find the maximum absolute difference between vectors v1 and v2, element wise 
double
Aggregation::maximum(double* v1, double* v2, int size)
{
   double value=0;
   for(int i=0; i<size; i++)
   {
      double value2=fabs(v1[i]-v2[i]);
      if(value2>value) value=value2;
   }
   return value;
}

//compute two ScoreItems, by the filed score

int itemCmp(const void *_a, const void *_b)
{   
    struct ScoreItem * c1 = (struct ScoreItem *) _a;
    struct ScoreItem * c2 = (struct ScoreItem *) _b;

    if(c1->score>c2->score) return -1;
    return 1;
}


//sort the pairs into decresing order 
void
Aggregation::topMIPairs(vector<int>* offset)
{

   struct ScoreItem * items=new struct ScoreItem[mData->getNumCols()*(mData->getNumCols()-1)/2];

   int which=0;
   for(int i=0; i<mData->getNumCols(); i++)
   {
      for(int j=i+1; j<mData->getNumCols(); j++)
      {
          items[which].score=mDIMatrix[i][j];
          items[which].i1   =i;
          items[which].i2   =j;
          which++;
      }
   }
   
   int size=mData->getNumCols()*(mData->getNumCols()-1)/2;
   int top=size;
   //if(size>10000) top=10000;
   //else top=size-3;
   
   
   //topK(items, 0, size-1, top, items[0].score);
   qsort(items, top, sizeof(struct ScoreItem), itemCmp); 

   for(int i=0; i<mData->getNumCols()*(mData->getNumCols()-1)/2; i++)
   {
      if(offset!=NULL)
      {
         cout<<(*offset)[items[i].i1]<<" "<<(*offset)[items[i].i2]<<" "<<items[i].score<<endl;
      }
      else
      {
         cout<<items[i].i1<<" "<<items[i].i2<<" "<<items[i].score<<endl;
      }

   }
}




//find the top K ranked pairs. 
void
Aggregation::topK(ScoreItem * items, int start, int end, int k, double pivot)
{
   if(k<=0) return;
   if(end<=start) return; 
   int partitionIndex=start;

   double maximum=-10000000;
   double minimum=100000000;
   //cout<<pivot<<" "<<start<<" "<<end<<" "<<k<<endl;
   for(int i=start; i<=end; i++)
   {
      double s=items[i].score;
      if(s>maximum) maximum=s;
      if(s<minimum) minimum=s; 
      if(s>pivot)
      {
         double score=items[partitionIndex].score;
         int       i1=items[partitionIndex].i1;
         int       i2=items[partitionIndex].i2;
         items[partitionIndex].score=items[i].score;
         items[partitionIndex].i1=items[i].i1;
         items[partitionIndex].i2=items[i].i2;

         items[i].score=score;
         items[i].i1=i1;
         items[i].i2=i2;
         partitionIndex++;
      }
   }
   //exit(0);
   int found=partitionIndex-start;
   
   if(found<k)
   {
      topK(items, partitionIndex, end, k-found, minimum+(pivot-minimum)*(k-found)/(end-partitionIndex+1));
   }
   else if(found>k)
   {
      topK(items, start, partitionIndex-1, k, pivot+(maximum-pivot)*(found-k)/(partitionIndex-start)); 
   }

}


