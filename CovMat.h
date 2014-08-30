#ifndef _COV_MAT_H_
#define _COV_MAT_H_
#include "MICData.h"
class CovMat
{
      MICData* mData;
      double** mMarginMean;
      double** mJointMat; 
      double** mCovMat;
      void calculateCov(); 
   public:
     CovMat(MICData *aData);
     double** getCovMat() { return mCovMat; }
};

#endif
