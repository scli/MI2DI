#ifndef _PSEUDO_COUNT_H_
#define _PSEUDO_COUNT_H_
#include "MICData.h"
class PseudoCount
{

     void allocateSpace();
     void uniformPseudoCount();
     void uniformContinuousPseudoCount();
     double** mPseudoJointMatrix;
     double** mPseudoMargin;
     double** mJointMatrix;
     double** mMargin;

     MICData* mData;
  public:
     PseudoCount(MICData* data);
     ~PseudoCount();
     void psuedoCount();
     double** getPSJoint() { return  mPseudoJointMatrix; }
     double** getPSMargin() { return mPseudoMargin; }



};



#endif
