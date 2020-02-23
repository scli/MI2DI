#ifndef _COV_MAT_H_
#define _COV_MAT_H_
#include "JointMat.h"


class CovMat
{
      JointMat* mData;
      double** mMargin;
      double** mJoint; 
      double** mCovMat;
      void calculateCov(); 
   public:
     CovMat(JointMat *aData);
     double** getCovMat() { return mCovMat; }
};

#endif
