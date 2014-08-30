#include <stdlib.h>
#include <iostream>
#include "AlignmentParser.h"
#include "MICParaMgr.h"


int* letter2number_map=new int[300];
char num2letter(int n)
{
  return number2letter_map[n];
}

int letter2num(char c)
{
  return letter2number_map[c];
}


AlignmentParser::AlignmentParser(SeqPool* sp)
{
   letter2number_map[256]= 0; 
   letter2number_map['-']= 1;
   letter2number_map['A']= 2;
   letter2number_map['C']= 3;
   letter2number_map['D']= 4;
   letter2number_map['E']= 5;
   letter2number_map['F']= 6;
   letter2number_map['G']= 7;
   letter2number_map['H']= 8;
   letter2number_map['I']= 9;
   letter2number_map['K']= 10;
   letter2number_map['L']= 11;
   letter2number_map['M']= 12;
   letter2number_map['N']= 13;
   letter2number_map['P']= 14;
   letter2number_map['Q']= 15;
   letter2number_map['R']= 16;
   letter2number_map['S']= 17;
   letter2number_map['T']= 18;
   letter2number_map['V']= 19;
   letter2number_map['W']= 20;
   letter2number_map['Y']= 21;
   letter2number_map['B']= -1; 
   letter2number_map['Z']= -1;
   letter2number_map['J']= -1; 
   letter2number_map['X']= -1;
   letter2number_map['U']= -1; 
   letter2number_map['O']= -1;
   letter2number_map['a']= -2; 
   letter2number_map['c']= -2; 
   letter2number_map['d']= -2; 
   letter2number_map['e']= -2; 
   letter2number_map['f']= -2; 
   letter2number_map['g']= -2; 
   letter2number_map['h']= -2; 
   letter2number_map['i']= -2; 
   letter2number_map['k']= -2; 
   letter2number_map['l']= -2; 
   letter2number_map['m']= -2; 
   letter2number_map['n']= -2; 
   letter2number_map['p']= -2; 
   letter2number_map['q']= -2; 
   letter2number_map['r']= -2; 
   letter2number_map['s']= -2; 
   letter2number_map['t']= -2; 
   letter2number_map['v']= -2; 
   letter2number_map['w']= -2; 
   letter2number_map['y']= -2; 
   letter2number_map['b']= -2; 
   letter2number_map['z']= -2; 
   letter2number_map['j']= -2; 
   letter2number_map['x']= -2; 
   letter2number_map['u']= -2; 
   letter2number_map['o']= -2; 
   letter2number_map['.']= -3;
  
 mSeqs = sp->getSeqs();
}


int
AlignmentParser::findIdentifier(string& pid)
{
   for(int i=0; i<mSeqs->size(); i++)
   {
     Seq* seq=(*mSeqs)[i];
     string &name=seq->mName;
     int pos=name.find('/');
     string pname=name.substr(0, pos);
     
     
     if(pid==pname)
     { 
       int pos2    =name.find('-', pos+1);
       
       mRangeStart =atoi((name.substr(pos+1, pos2-(pos+1))).c_str());
       mRangeEnd   =atoi((name.substr(pos2+1, name.size()-(pos2+1))).c_str());
       
       return i;
     }
   }

   cerr<<"Protein "<<pid<<" is not found in the input alignment block!"<<endl;
   return -1;
}

vector<string>*
AlignmentParser::detectUsefulColumn(const char* proteinIdentifier)
{
  string id=string(proteinIdentifier);
  int int_id=findIdentifier(id);
  string& seq=(*mSeqs)[int_id]->mSeq;

  int len=seq.size();
  mColumns=new vector<int>();
  mOffsets=new vector<int>();
  int offset=mRangeStart;
  int firstSkipped=0;
  
  for(int i=0; i<len; i++)
  {
       int residuecode=letter2num(seq[i]);
       if (residuecode == 0)
       {
       }
       if (residuecode == -1)
       {
       }
       if (residuecode > 1)
       {
          if(MICParaMgr::_ALIGN_SKIP_END_==false)
          {
             mColumns->push_back(i);
             mOffsets->push_back(offset);
          }
          else if(firstSkipped>MICParaMgr::_ALIGN_SKIP_END_COUNTER_)
          {
          
            mColumns->push_back(i);
            mOffsets->push_back(offset);
          }
       }
       if (residuecode == -2 || residuecode > 1)
       {
         firstSkipped++;  
         offset++;
       }
   }
    
   if(MICParaMgr::_ALIGN_SKIP_END_==true)
   {
      for(int i=0; i<MICParaMgr::_ALIGN_SKIP_END_COUNTER_; i++)
      {
         mColumns->pop_back();
         mOffsets->pop_back();
      }
   }
   mAlignedStr = projectAlginedRegion();  
   return NULL;
}

vector<Seq* >* 
AlignmentParser::filterSeqs()
{

  vector<Seq* >* filteredSeqs=new vector<Seq* >();

  for(int i=0; i<mSeqs->size(); i++)
  {
     bool toSkip=false;
     for(int c=0; c<mColumns->size(); c++)
     {
        char cc=((*mSeqs)[i]->mSeq)[(*mColumns)[c]];
        if(letter2number_map[cc]==-1)
        {
            toSkip=true; break;
        }
     }
     toSkip=false;
     if(!toSkip)
     {
       filteredSeqs->push_back((*mSeqs)[i]);
     }
  }
  return filteredSeqs;

}


void
AlignmentParser::filterColumns(vector<Seq*>* filteredSeqs, bool* isUseful)
{

   if(!(MICParaMgr::_ALIGN_SKIP_END_))
       return;


  bool headGaps=true;
  double thres=MICParaMgr::_ALIGN_FILTER_END_THRES_;
  //0.10;
  
  for(int i=0; i<mColumns->size(); i++)
  {
     double variance=0;
     for(int j=0; j<filteredSeqs->size(); j++)
     {
       char c=((*filteredSeqs)[j]->mSeq)[(*mColumns)[i]];
       if(c<'A' || c>'Z') variance+=1;
       else
       {
          int index=AA2Index[c-'A'];
          if(index<0) variance+=1;
       }
     }
     //if(variance/filteredSeqs->size()>=0.2)
     //cerr<<(variance/filteredSeqs->size())<<" hahaha "<<endl;
     if(headGaps && variance>=thres*filteredSeqs->size())
     {
         isUseful[i]=false;
     }
     else 
         headGaps=false;
     if(variance>=MICParaMgr::_ALIGN_FILTER_COL_THRES_*filteredSeqs->size())  
         isUseful[i]=false;
  }//for
  
  
  
  
  
  bool tailGaps=true;
  for(int i=mColumns->size()-1; i>=0; i--)
  {
     double variance=0;
     for(int j=0; j<filteredSeqs->size(); j++)
     {
       char c=((*filteredSeqs)[j]->mSeq)[(*mColumns)[i]];
       if(c<'A' || c>'Z') variance+=1;
       else
       {
          int index=AA2Index[c-'A'];
          if(index<0) variance+=1;
       }
     }
     if(tailGaps && variance>=thres*filteredSeqs->size())
     {
         isUseful[i]=false;
     }
     else 
        tailGaps=false;
  }//for




}

vector<string>*
AlignmentParser::projectAlginedRegion()
{
  vector<Seq* >* filteredSeqs=filterSeqs();

  //int counter=0;
  bool* isUseful=new bool[mColumns->size()];
  for(int i=0; i<mColumns->size(); i++)
  {
      isUseful[i]=true;
  }

  filterColumns(filteredSeqs, isUseful);
  
  vector<int>* cols=new vector<int>();
  vector<int>* offs=new vector<int>();
  
  for(int i=0; i<mColumns->size(); i++)
  {
     if(isUseful[i])
     {
        cols->push_back((*mColumns)[i]);
        offs->push_back((*mOffsets)[i]);
     }
  }
  mColumns=cols;
  mOffsets=offs;
  mAlignedStr=new vector<string>();//mSeqs->size());
  
  for(int i=0; i<filteredSeqs->size(); i++)
  {
     string seq=string(mColumns->size(), 'A');
     for(int c=0; c<mColumns->size(); c++)
     {
        seq[c]=((*filteredSeqs)[i]->mSeq)[(*mColumns)[c]];
     }
     mAlignedStr->push_back(seq);
  //   cout<<((*mAlignedStr)[mAlignedStr->size()-1]).c_str()<<endl;
  }
  return mAlignedStr;  
}


