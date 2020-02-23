#include <stdlib.h>
#include "MICParaMgr.h"
#include <iostream>
#include "ExpEqlPipe.h"
#include "DataPartition.h"
#include "CovMat.h"
#include "IntegratedAgg.h"
#include <math.h>
using namespace std;

int main(int argc, char** argv)
{

  if(argc!=2)
  {
    cerr<<"Usage: "<<argv[1]<<" file_name"<<endl;
    exit(0);
  }
  MICParaMgr::_NUM_PARTITION_CLUSTERS_=1;  
  ExpressionParser* mParser    =new ExpressionParser(argv[1]);
  DataPartition*    partition  =new DataPartition(mParser->getExpressions(), 
                                                  mParser->getNumReplications(),
                                                  mParser->getNumGenes());
  partition->partition();
  CntnsData*       conData =partition->getCntnsData();
  JointMat*        conMat  =new JointMat(conData);
  conData->setMargin(conMat->margin());
  double** c = conMat->joint();

  conData->setJointMatrix(conMat->joint());
  
  MICParaMgr::_PSEUDO_COUNT_WEIGHT_=0;

  PseudoCount*     conPseudo=new PseudoCount(conMat, NULL, _CONTINUOUS_);
  //conData->setPseudoMargin(conPseudo->getPSMargin());
  //conData->setPseudoJointMatrix(conPseudo->getPSJoint());
  
  
  //JointMat conData=new JointMat(conPseudo->pseudoMat());
  
  CovMat*            conCov=new CovMat(conPseudo->pseudoMat());
  cout<<"i reach here "<<endl;
  
  double** mat = conCov->getCovMat();
  int ncol=conData->getNumCols();
  double* values=new double[ncol];
  for(int i=0; i<ncol; i++)
  {
     for(int j=0; j<ncol; j++)
     {
         cout<<fabs(mat[i][j]/sqrt(mat[i][i]*mat[j][j]))<<" \t";
     }
     cout<<endl;
  }

}
