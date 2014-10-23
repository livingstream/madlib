#ifndef _H_GIBBS_SCHEDULER_H
#define _H_GIBBS_SCHEDULER_H
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "../../common/timer.h"
#include "gibbs_cgla.h"
#include "../../state/groundclause.h"
#include "../convergencetest.h"
#include "../gelmanconvergencetest.h"
#include "../gibbsparams.h"
#include <iostream>
using namespace std;

class GibbsGist;

class GibbsScheduler {
public:
   size_t numAtoms;
   size_t chainId;
   GibbsGist *st;

   GibbsScheduler(size_t chainId_, size_t inNumAtoms, GibbsGist *inState);

   ~GibbsScheduler();

   static void  *callMemberFunction(void *arg);

   static bool genTruthValueForProb(const double &p);

   double getProbabilityOfPred(const size_t &predIdx);

   void *performGibbsStep();

   void updateWtsForGndPreds(vector<size_t> &gndPredIndices);

   void gndPredFlippedUpdates(const size_t &gndPredIdx);

};
#endif
