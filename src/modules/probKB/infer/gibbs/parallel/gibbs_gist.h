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
   GibbsCGLA *gGla;
   vector<GibbsScheduler *> gibbsVec;
   vector<vector<double> > numTrueTemp;
   VariableState *varst;
   vector<double> probs;

   vector<vector<double> > numTrue;
   size_t numAtoms;
   size_t numClauses;
   size_t numChains;

public:
   GibbsGist(size_t numAtoms_, size_t numClauses_, VariableState *varst_);
  
   ~GibbsGist();
 
   void initTruthValues();
   
   void infer();

   void finalize();
};
#endif
