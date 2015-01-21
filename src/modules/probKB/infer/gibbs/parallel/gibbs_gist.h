#ifndef _H_GIBBS_GIST_H
#define _H_GIBBS_GIST_H
#include<vector>
#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include "gibbs_scheduler.h"
#include "gibbs_cgla.h"
#include "../../state/variablestate.h"
using namespace std;
#define NUM_THREADS 10

class GibbsGist {
public:
   pthread_t threads[NUM_THREADS];
   vector<GibbsScheduler *> localSchedulers;
   GibbsCGLA *gGla;
   vector<GibbsCGLA *> glaVec;
   vector<GibbsScheduler *> gibbsVec;
   vector<vector< bool> > truthValues;
   vector<vector<double> > wtsWhenFalse;
   vector<vector<double> > wtsWhenTrue;
   vector<vector<double> > numTrue;
   vector<vector<double> > numTrueTemp;
   vector<vector<bool> > affectedGndPredFlag;
   vector<vector<int> > numTrueLits;
   vector<vector<size_t> > affectedGndPredIndices;
   VariableState *varst;
   vector<double> probs;

   size_t numAtoms;
   size_t numClauses;
   size_t numChains;

public:
   GibbsGist(size_t numAtoms_, size_t numClauses_, VariableState *varst_);

   void initTruthValues();
   
   void infer();

   void finalize();
};
#endif
