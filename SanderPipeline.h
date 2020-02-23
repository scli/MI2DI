#ifndef _SANDER_PIPELINE_H_
#define _SANDER_PIPELINE_H_

#include "SanderPipeline.h"
#include "MICData.h"
#include "AlignData.h"
#include "JointMat.h"
#include "PseudoCount.h"
#include "Inversion.h"
#include "SanderAgg.h"
#include "CovMat.h"
#include "LatentMap.h"
class SanderPipeline
{
    AlignmentParser* mParser;

  public:
    SanderPipeline(const char* seqs, const char* theProtein);
    void solve();


};

#endif
