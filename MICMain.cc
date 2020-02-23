#include "SanderPipeline.h"
#include "ExpEqlPipe.h"

#include <stdlib.h>
#include <iostream>
using namespace std;

int main(int argc, char** argv)
{

  cout<<argc<<endl; 
  if(argc!=4 && argc!=3)
  {
    cerr<<"Usage: "<<argv[0]<<" .fasta_file pdbid output"<<endl;
    exit(0);
  } 
  if(argc==4)
  {
    SanderPipeline* sp=new SanderPipeline(argv[1], argv[2]);
    sp->solve(); 
    exit(0);
  }
  else
  {
    ExpEqlPipe*  exp=new ExpEqlPipe(argv[1]);
    exp->solve();
    exit(0);
  }
  
  //vector<string>* alignedStr=parser->detectUsefulColumn(argv[2]);
  //DirectKernel* dk=new DirectKernel(alignedStr, parser->getOffsets());

}
