#ifndef _INVERSION_H_
#define _INVERSION_H_
#include "MICData.h"
class Inversion
{

     double** mCovMat;
     double** mInvMat;
     MICData* mData;
     double** mMargin;
     double** mJointProbMatrix;
     double** mInvCorrMat;
     void invertByRemovingLastColumns();
     void invertByLasso();
     double** inverseMat(double** mat, int size);
     bool test_cholesky(double *a, const int n);

     void cholesky(double **a, double **l, int n);
     double** inverseByRegulizing(double** mat, int size, double eta);
     void invertByLSRegularize();

   public:
     Inversion(MICData* data);
     ~Inversion();
     void inverse();
     double** getInvMat() { return mInvMat; }
     double** getInvCorr();



};
#endif
