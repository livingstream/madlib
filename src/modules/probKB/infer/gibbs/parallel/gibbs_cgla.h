#ifndef H_GIBBS_CGLA_H
#define H_GIBBS_CGLA_H
#include <iostream>
#include <vector>
#include "../gibbsparams.h"
#include "../convergencetest.h"
#include "../gelmanconvergencetest.h"
using namespace std;

class GibbsGist;

class GibbsCGLA {
public:
   bool done;
   int sample;
   double secondsElapsed;
   double startTimeSec;
   double currentTimeSec;
   bool burnIn;
   size_t numAtoms;
   size_t numChains;
   Timer timer;
   GibbsParams *para;
   GibbsGist *state;
   GelmanConvergenceTest **burnConvergenceTests;
   ConvergenceTest **gibbsConvergenceTests;
   GibbsCGLA(size_t numAtoms_, size_t numChains_, GibbsGist *state_);

   ~GibbsCGLA();

   void merge(GibbsCGLA *cgla);

   bool shouldIterate();

   void initConvergenceTests();

   void deleteConvergenceTests();
};
#endif
