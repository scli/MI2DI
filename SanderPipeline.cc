#include <stdlib.h>
#include <iostream>
#include "SanderPipeline.h"
using namespace std;

SanderPipeline::SanderPipeline(const char* seqs, const char* theProtein)
{


  SeqPool* sp=new SeqPool();
  sp->readSeqs(seqs); 
  mParser=new AlignmentParser(sp);
  mParser->detectUsefulColumn(theProtein);

}

void
SanderPipeline::solve()
{

  AlignData*       align  =new AlignData(mParser);
  DiscrData*       data =align->getData();
  JointMat* mat   =new JointMat(data);
  data->setMargin(mat->getColMean());
  data->setJointMatrix(mat->getJointMean());
  
  PseudoCount*     pseudo=new PseudoCount(data);
  data->setPseudoMargin(pseudo->getPSMargin());
  data->setPseudoJointMatrix(pseudo->getPSJoint());
 
  CovMat*          cov=new CovMat(data);
  data->setCovMat(cov->getCovMat());

  Inversion*      inverse=new Inversion(data);
  inverse->inverse();
  data->setInverseMat(inverse->getInvMat());

  SanderAgg*     agg   =new SanderAgg(data);
  agg->topMIPairs(mParser->getOffsets());

}
