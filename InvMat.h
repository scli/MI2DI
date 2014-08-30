#ifndef _INVERSION_MAT_H_
#define _INVERSION_MAT_H_
#include "MICData.h"
class InvMat
{

     double** mCovMat;
     double** mInvMat;
     void invertByLasso();
     void invertByLSRegularize();
 
     void eigenDecompose(double** inputMat, double** eigenVec, double* eigenVal);
     void cholesky(double **a, double **l, int n);
     double** inverseByRegulizing(double** mat);
     double mEta;
     int mDim;
     int mNumTrans;
   public:
     InvMat(int n);
     ~InvMat();
     void read(char* filename);
     void inverse(double eta);
     void write(char* filename);
     void deconvolution(double** inputMat);

};
#endif
