#include "InvMat.h"
#include <iostream>
#include <stdlib.h>
using namespace std;
int main(int argc, char** argv)
{

   if(argc!=5 && argc!=4)
   {
      cerr<<argv[0]<<" matrix_name output_name dimension eta"<<endl;
      exit(0);
   }

   InvMat* mat=new InvMat(atoi(argv[3]));
   mat->read(argv[1]);
   mat->deconvolution(NULL);
   //mat->inverse(atof(argv[4]));
   mat->write(argv[2]);
}
