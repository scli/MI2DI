#ifndef _EXPRESSION_PARSER_H_
#define _EXPRESSION_PARSER_H_
#include <vector>
#include <string>
using namespace std;


class ExpressionParser
{

      vector<string >* mGeneNames;
      int mNumReplicates;
      double** mExpressions;
      void readFile(const char* file);
    public:
      ExpressionParser(const char* filename);
      ~ExpressionParser();
      int getNumGenes();
      int getNumReplications();
      double** getExpressions();
};

#endif
