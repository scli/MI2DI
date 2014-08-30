#include "MICData.h"
#include <math.h>
#include <iostream>
using namespace std;
MICData::MICData()
{

}

double
MICData::mutualInformation(int pos1, int pos2, double** joint, double* margin1, double* margin2)
{
   
  double tiny = 1.0e-100;
  
  int dim1=getHiddenDim(pos1);
  int dim2=getHiddenDim(pos2);

  int off1=getHiddenOffset(pos1);
  int off2=getHiddenOffset(pos2);

  double rev=0;

  for(int i=0; i<dim1; i++)
  {
      for(int j=0; j<dim2; j++)
      {
        if(joint[off1+i][off2+j]!=0)
        rev+=(joint[off1+i][off2+j]*log((tiny+joint[off1+i][off2+j])/(tiny+margin1[i]* margin2[j])));

  //      cout<<joint[off1+i][off2+j]<<" "<<" "<<margin1[i]<<" "<<margin2[j]<<endl;
      }
  }
  return rev;
}
double
maximum(double a, double b)
{
   if(a>b)
       return a;
   return b;
}


