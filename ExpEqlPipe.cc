#include <stdlib.h>
#include <iostream>
#include "ExpEqlPipe.h"
#include "DataPartition.h"
#include "CovMat.h"
#include "IntegratedAgg.h"
using namespace std;

ExpEqlPipe::ExpEqlPipe(const char* file)
{

  mParser=new ExpressionParser(file);
}

void
ExpEqlPipe::solve()
{

  DataPartition*   partition  =new DataPartition(mParser->getExpressions(), 
                                                 mParser->getNumReplications(),
                                                 mParser->getNumGenes());
  partition->partition();
  DiscrData*       disData =partition->getDiscrData();
  CntnsData*       conData =partition->getCntnsData();
  
  //cout<<"partitioned the data!!!"<<endl;
  JointMat*       disMat   =new JointMat(disData);
  JointMat*       conMat   =new JointMat(conData);
  ///cout<<"Computed the distance matrices"<<endl;
  //cout.flush();
  disData->setMargin(disMat->margin());
  conData->setMargin(conMat->margin());
  
  disData->setJointMatrix(disMat->joint());
  conData->setJointMatrix(conMat->joint());
  
  PseudoCount*     disPseudo=new PseudoCount(disMat, NULL, _DISCRETE_);
  //disData->setPseudoMargin(disPseudo->getPSMargin());
  //disData->setPseudoJointMatrix(disPseudo->getPSJoint());

  PseudoCount*     conPseudo=new PseudoCount(conMat, NULL, _CONTINUOUS_);
  //conData->setPseudoMargin(conPseudo->getPSMargin());
  //conData->setPseudoJointMatrix(conPseudo->getPSJoint());
  
  //CovMat*            disCov=new CovMat(disPseudo->pseudoMat());
  //CovMat*            conCov=new CovMat(conPseudo->pseudoMat());
  
  //disData->setCovMat(disCov->getCovMat());         
  //conData->setCovMat(conCov->getCovMat());         
  
  Inversion*      disInverse=new Inversion(disPseudo->pseudoMat());
  disInverse->inverse();
  disData->setInverseMat(disInverse->getInvMat());

  Inversion*      conInverse=new Inversion(conPseudo->pseudoMat());
//  conInverse->inverse();
//  conData->setInverseMat(conInverse->getInvMat());
  cout<<"I reach here"<<endl;
  cout.flush();
  conData->setInvCorrMat(conInverse->getInvCorr());
   
  IntegratedAgg*     agg   =new IntegratedAgg(disData, conData, disPseudo->pseudoMat(), disInverse->getInvMat());
  agg->topMIPairs(NULL);

}
