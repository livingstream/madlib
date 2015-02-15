#ifndef _H_GIBBS_SCHEDULER_H
#define _H_GIBBS_SCHEDULER_H
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>     
#include <sstream>
#include "../../common/timer.h"
#include "gibbs_cgla.h"
#include "../../state/groundclause.h"
#include "../convergencetest.h"
#include "../gelmanconvergencetest.h"
#include "../gibbsparams.h"
#include <iostream>
#include <fstream>
#include<pthread.h>
using namespace std;

class GibbsGist;

class GibbsScheduler {
public:
   static pthread_mutex_t lock;
   size_t numAtoms;
   size_t numClauses;
   const size_t chainId;
   GibbsGist *st;

   vector<bool> loc_truthValues;
   vector<double> loc_wtsWhenFalse;
   vector<double> loc_wtsWhenTrue;
   vector<double> loc_numTrue;
   vector<double> loc_numTrueTemp;
   vector<bool> loc_affectedGndPredFlag;
   vector<int> loc_numTrueLits;
   vector<size_t> loc_affectedGndPredIndices;
   
   GibbsScheduler(size_t chainId_, size_t inNumAtoms, size_t inNumClauses, GibbsGist *inState);

   ~GibbsScheduler();
 
   void init();

   void initTruthValuesAndWts();

   void randomInitGndPredsTruthValues();

   static void  *callMemberFunction(void *arg);

   bool genTruthValueForProb(const double &p);

   double getProbabilityOfPred(const size_t &predIdx);

   void *performGibbsStep();

   void updateWtsForGndPreds(vector<size_t> &gndPredIndices);

   void gndPredFlippedUpdates(const size_t &gndPredIdx);
   static inline void printLog(int threadId, char* msg);

};
#endif
