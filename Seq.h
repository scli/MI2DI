#ifndef _SEQ_H_
#define _SEQ_H_
#include <string>
#include <vector>
using namespace std;
class Seq
{

   public:
    string mSeq;
    string mName;
    Seq(string s, string name);
    Seq();
    ~Seq();
};

class SeqPool
{
    vector<Seq*> *mSeqs;
   public:
     SeqPool(); 
     ~SeqPool();

     void insert(Seq* s);
     vector<Seq *>* getSeqs() { return mSeqs; } 
     void readSeqs(const char* filename);

};
#endif
