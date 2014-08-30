#include "Seq.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
using namespace std;
Seq::Seq(string s, string name)
{
   mSeq=s;
   mName=name;
}

Seq::Seq()
{
  mSeq="";
  mName="";
}

Seq::~Seq()
{
}


SeqPool::SeqPool()
{
   mSeqs=new vector<Seq *>();
}

SeqPool::~SeqPool()
{
  for(int i=0; i<mSeqs->size(); i++)
  {
     delete (*mSeqs)[i];
  }
  delete mSeqs;
}

void
SeqPool::insert(Seq* s)
{
  mSeqs->push_back(s);
}

//read fasts format
void
SeqPool::readSeqs(const char* filename)
{
  ifstream input(filename);

  if(!input)
  {
    cerr<<"Can not open file "<<filename<<endl;
    exit(0);
  }

  char buf[1000];
  char* token;
  while(!input.eof())
  {
    input.getline(buf, 1000);
    token=strtok(buf, " \t\n");
    if(token==NULL) 
        continue;
    
    if(token[0]=='>')
    {   
       mSeqs->push_back(new Seq());
      (*mSeqs)[mSeqs->size()-1]->mName = string(token+1);      
      continue; 
    }
    else 
    {
       (*mSeqs)[mSeqs->size()-1]->mSeq=(*mSeqs)[mSeqs->size()-1]->mSeq+string(token);
    }
   }
}

