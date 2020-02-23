#include "ExpressionParser.h"
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
using namespace std;
ExpressionParser::ExpressionParser(const char* name)
{
  mGeneNames=new vector<string>(); 
  readFile(name);
}

ExpressionParser::~ExpressionParser()
{
}

void
ExpressionParser::readFile(const char* filename)
{
  ifstream input(filename);

  if(!input)
  {
    cerr<<"Can not open file "<<filename<<endl;
    exit(0);
  }

  char buf[100000];
  char* token;
  while(!input.eof())
  {
    input.getline(buf, 100000);
    token=strtok(buf, ", \t\n");
    if(token==NULL) 
        continue;
    while(token!=NULL)
    {
       mGeneNames->push_back(string(token));
       token=strtok(NULL, ", \t\n");
    }
    break;
  }
  vector<double*>* expressions=new vector<double *>(); 
   
  while(!input.eof())
  {
    input.getline(buf, 100000);
    token=strtok(buf, ", \t\n");
    if(token==NULL) 
        continue;
    double* expression=new double[mGeneNames->size()];
    int index=0;
    
    for(int index=0; index<mGeneNames->size(); index++)
    {
      expression[index]=atof(token);
      token=strtok(NULL, ", \t\n");
    }
    expressions->push_back(expression);
  }


  //store the data into a 2D array
  mExpressions=new double *[expressions->size()];
  for(int i=0; i<expressions->size(); i++)
  {
     mExpressions[i]=(*expressions)[i];
  }
  mNumReplicates=expressions->size();
  delete expressions;
}

int
ExpressionParser::getNumGenes()
{
   return mGeneNames->size();
}
int 
ExpressionParser::getNumReplications()
{
   return mNumReplicates;
}
double**
ExpressionParser::getExpressions()
{
  return mExpressions;
}
