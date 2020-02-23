#ifndef _PSEUDO_COUNT_H_
#define _PSEUDO_COUNT_H_
#include "JointMat.h"
#include "LatentMap.h"
class PseudoCount
{

     void allocateSpace();
     void uniformPseudoCount();
     void uniformContinuousPseudoCount();
     double** mPseudoJointMatrix;
     double** mPseudoMargin;
     double** mJointMatrix;
     double** mMargin;
     LatentMap* mMap;
     JointMat*  mData;
     JointMat*  mPSMat;
     DataType   mDataType;
  public:
     PseudoCount(JointMat* data, LatentMap* map, DataType t);
     ~PseudoCount();
     void psuedoCount();
     double** getPSJoint()  { return  mPseudoJointMatrix; }
     double** getPSMargin() { return mPseudoMargin; }
     JointMat* pseudoMat()  { return mPSMat;        }
     //DataType  dataType()   { return mDataType; }
     //DataType  dataType(DataType t)   { return mDataType=t; }

};



#endif
