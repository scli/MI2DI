#include <stdlib.h>
#include <iostream>
#include "SanderPipeline.h"
#include <math.h>
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
  //read in alignment data, and it will produce a k by n matrix. 
  AlignData*       align  =new AlignData(mParser);
  
  
  LatentMap*       map    =align->getMap();
  DiscrData*       data   =align->getData();
  //  
  cout<<"step 1"<<endl; 

  //create the joint data matrix, it will create a n by n matrix,  
  JointMat*   mat2   =new JointMat(data);
  //check if mat2 is sym
  int numclus=mat2->numClus();
  for(int i=0; i<numclus; i++)
  {
     for(int j=0; j<numclus; j++)
     {
	 if(fabs(mat2->joint()[i][j]-mat2->joint()[j][i])>=0.0001)
	 {
            cerr<<"errors ...."<<mat2->joint()[i][j]<<" "<<mat2->joint()[j][i]<<endl;
	    exit(0);
	 }
     }
  }
  
  
  cout<<"step 2"<<endl; 
  cout.flush();
  JointMat*    mat   =map->convert(mat2);
  cout<<"step 3"<<endl; 
  cout.flush();


  data->setMargin(mat->margin());
  cout<<"step 4"<<endl; 
  data->setJointMatrix(mat->joint());
  cout<<"step 5"<<endl; 
  
  //add the pseudo count the matrix
  PseudoCount*     pseudo=new PseudoCount(mat, map, _DISCRETE_);
  //data->setPseudoMargin(pseudo->getPSMargin());
  //data->setPseudoJointMatrix(pseudo->getPSJoint());
 
  //compute the covariance matrix
  //CovMat*          cov=new CovMat(pseudo->pseudoMat());
  //data->setCovMat(cov->getCovMat());

  //compute the inversion
  Inversion*      inverse=new Inversion(pseudo->pseudoMat());
  cout<<"step 6"<<endl; 
  inverse->inverse();
  data->setInverseMat(inverse->getInvMat());
  cout<<"step 7"<<endl; 

  //compute the direct information
  SanderAgg*     agg   =new SanderAgg(data, pseudo->pseudoMat(), inverse->getInvMat());
  cout<<"step 8"<<endl; 
  
  
  //output the result
  agg->topMIPairs(mParser->getOffsets());

}
