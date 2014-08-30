#ifndef _PARTIAL_CORR_
#define _PARTIAL_CORR_
#include "Seq.h"

#include <string>
#include <vector>

using namespace std;

const static int AA2Index[26] =
{ 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16,
 20, 19, 17, 20, 18, 20};


const static char number2letter_map[22]={
   '-',
   'A',
   'C',
   'D',
   'E',
   'F',
   'G',
   'H',
   'I',
   'K',
   'L',
   'M',
   'N',
   'P',
   'Q',
   'R',
   'S',
   'T',
   'V',
   'W',
  'Y'};

char num2letter(int n);
int letter2num(char c);

class AlignmentParser
{

     int findIdentifier(string& pid);
     vector<int>* mColumns;
     vector<int>* mOffsets;
     vector<Seq* >* mSeqs;
     int mRangeStart;
     int mRangeEnd;
     vector<string>* mAlignedStr;
     vector<string>* projectAlginedRegion();
     vector<Seq* >* filterSeqs();
     void filterColumns(vector<Seq*>* filteredSeqs, bool* isUseful);
   public:
     AlignmentParser(SeqPool* sp); 
     ~AlignmentParser();

    
     vector<string>* detectUsefulColumn(const char* proteinIdentifier); 
     vector<int>* getOffsets() { return mOffsets; }    
     vector<string>* getAlignStr() { return mAlignedStr; }

};
#endif
