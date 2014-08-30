#ifndef _EXPRESSION_EQUAL_PIPE_
#define _EXPRESSION_EQUAL_PIPE_

#include "ExpEqlPipe.h"
#include "MICData.h"
#include "AlignData.h"
#include "JointMat.h"
#include "PseudoCount.h"
#include "Inversion.h"
#include "SanderAgg.h"
#include "ExpressionParser.h"
class ExpEqlPipe
{
    ExpressionParser* mParser;

  public:
    ExpEqlPipe(const char* filename);
    void solve();


};

#endif
