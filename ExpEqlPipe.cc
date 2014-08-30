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
  disData->setMargin(disMat->getColMean());
  conData->setMargin(conMat->getColMean());
  
  disData->setJointMatrix(disMat->getJointMean());
  conData->setJointMatrix(conMat->getJointMean());
  
  PseudoCount*     disPseudo=new PseudoCount(disData);
  disData->setPseudoMargin(disPseudo->getPSMargin());
  disData->setPseudoJointMatrix(disPseudo->getPSJoint());

  PseudoCount*     conPseudo=new PseudoCount(conData);
  conData->setPseudoMargin(conPseudo->getPSMargin());
  conData->setPseudoJointMatrix(conPseudo->getPSJoint());
  
  CovMat*            disCov=new CovMat(disData);
  CovMat*            conCov=new CovMat(conData);
  
  disData->setCovMat(disCov->getCovMat());         
  conData->setCovMat(conCov->getCovMat());         
  
  Inversion*      disInverse=new Inversion(disData);
  disInverse->inverse();
  disData->setInverseMat(disInverse->getInvMat());

  Inversion*      conInverse=new Inversion(conData);
//  conInverse->inverse();
//  conData->setInverseMat(conInverse->getInvMat());
  cout<<"I reach here"<<endl;
  cout.flush();
  conData->setInvCorrMat(conInverse->getInvCorr());
   
  IntegratedAgg*     agg   =new IntegratedAgg(disData, conData);
  agg->topMIPairs(NULL);

}
