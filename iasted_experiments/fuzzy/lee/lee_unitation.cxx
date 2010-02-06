/* ----------------------------------------------------------------------------
  lee_unitation.cxx

  Neal Richter
  Montana State University
    
  Copyright (c) 2005 Neal Richter

 DESCRIPTION:
   1+3 unitation GA with lee fuzzy adaptive rules

   
   Uses GAlib library and example source as a guide/template
   Copyright (c) 1995-199 Massachusetts Institute of Technology and Matthew Wall.
   http://lancet.mit.edu/ga/
   Thanks MATT!
---------------------------------------------------------------------------- */
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <math.h>
#include <sys/timeb.h>
#include <ga/GA1DBinStrGenome.h>

#include "GASimpleGA_1p3.h"

#include "fuzzy.h"

#define GENOME_SIZE   20
#define POPULATION_SIZE   10
#define NUM_GEN 1000

static int obj_case = 0;
static int global_population_size = POPULATION_SIZE;

#define NUM_SAMPLES_MAX	       2500
#define NUM_SAMPLES_MIN	       250
#define NUM_SAMPLES_DEFAULT	   500

#define FUZZY_HIGH           3
#define FUZZY_MED            2
#define FUZZY_LOW            1

#define MUT_LOW             (float)(1/(float)GENOME_SIZE)
#define MUT_MED             (float)(2/(float)GENOME_SIZE)
#define MUT_HIGH            (float)(4/(float)GENOME_SIZE)

#define MUT_DELTA_LOW       0.5
#define MUT_DELTA_MED       1.5
#define MUT_DELTA_HIGH      2.5

#define MUT_UPPER           0.5
#define MUT_LOWER           (float)(1/((float)2*GENOME_SIZE))

#define CROSS_LOW           0.95
#define CROSS_MED           0.97
#define CROSS_HIGH          1.00

#define CROSS_DELTA_LOW     0.5
#define CROSS_DELTA_MED     1.0
#define CROSS_DELTA_HIGH    1.75

#define CROSS_UPPER           1.0
#define CROSS_LOWER           0.75

//global variables -- for convienence
int current_generation;
float current_score_mean;
float current_score_stddev;

int   supreme_fit_count;
int   vhigh_fit_count;
int   high_fit_count;
int   med_fit_count;
int   low_fit_count;
int   vlow_fit_count;

int evals;


int max_objfunc_fitness = 0;

#include "unitation_header.h"


//functions
float Objective(GAGenome &);
int  lee_fuzzy_adapt(float *, float *, float, float, float, float);

int
main(int argc, char** argv) 
{
  time_t time_now;
  struct timeb mytime;
  unsigned int seed = 0;
  int stop_evolving = FALSE;

  
  int    unchanged_best_fitness_count = 0;
  float  new_pcrossover = 0;
  float  new_pmutate = 0;
  float  old_best_fitness = 0;
  float  current_best_fitness = 0;
  float  current_worst_fitness = 0;
  float  current_score_variance = 0;
  float  best_fitness_save = 0;
  float  fitness_delta_since_save = 0;
  int    fuzzy_adapt_action =  FALSE;

  if((argc != 2) && (argc != 3))
  {
    cout << "usage: "<< argv[0] << "objfunc#  [popsize]" << endl;
    cout << "\t\tobjfuncs:\n";
    cout << "\t\t\t:1 = ONEMAX\n";
    cout << "\t\t\t:2 = NEEDLE\n";
    cout << "\t\t\t:3 = DECTRAP\n";
    cout << "\t\t\t:4 = DECTWOTRAP\n";
    cout << "\t\t\t:5 = TWOTRAP\n";
    exit(1);
  }
  
  obj_case = atoi(argv[1]);
  if((argc == 3) && (strlen(argv[2]) > 0))
    global_population_size = atoi(argv[2]);
  
  cout << "Unitation" << GENOME_SIZE << " bit mutation only\n\n";
  cout << "GASimpleGA_1p3 (non-overlapping populations)\n";
  cout << "Lee & Takagi Fuzzy Adaptive Rules for mutation rate control\n";
  cout << endl;
  
  //default random seed
  ftime(&mytime);
  seed = (unsigned int)mytime.millitm;

  GA1DBinaryStringGenome worst(GENOME_SIZE, Objective);

  switch(obj_case)
  {

      case 1:// ONEMAX
          cout << "ONEMAX Unitation Obj" << endl;
          max_objfunc_fitness = GENOME_SIZE;
          worst.set(0,GENOME_SIZE);
          worst.unset(0,GENOME_SIZE);
          break;
      case 2:// NEEDLE
          cout << "NEEDLE Unitation Obj" << endl;
          max_objfunc_fitness = GENOME_SIZE;
          worst.set(0,GENOME_SIZE);
          worst.unset(0,GENOME_SIZE);
          worst.set(0,1);
          break;
      case 3:// DECTRAP
          cout << "DECTRAP Unitation Obj" << endl;
          max_objfunc_fitness = GENOME_SIZE;
          worst.set(0,GENOME_SIZE);
          worst.unset(0,GENOME_SIZE);
          worst.set(0,2);
          break;
      case 4:// DECTWOTRAP
          cout << "DECTWOTRAP Unitation Obj" << endl;
          max_objfunc_fitness = GENOME_SIZE;
          worst.set(0,GENOME_SIZE);
          worst.unset(0,GENOME_SIZE);
          worst.set(0,9);
          break;
      case 5:// TWOTRAP
          cout << "TWOTRAP Unitation Obj" << endl;
          max_objfunc_fitness = GENOME_SIZE;
          worst.set(0,GENOME_SIZE);
          worst.unset(0,GENOME_SIZE);
          worst.set(0,9);
          break;

      default:
          cout << "NO Obj func defined" << endl;
          exit(1);

  }


  GA1DBinaryStringGenome genome(GENOME_SIZE, Objective);
  GASimpleGA_1p3 ga(genome);

  ga.elitist(gaFalse);
  cout << "No Elitism." << endl;

  ga.populationSize(global_population_size);
  ga.nGenerations(NUM_GEN);
  cout << "PopulationSize = " <<  ga.populationSize(); 
  cout << " Num Generations = " << ga.nGenerations() << endl;
  
  ga.recordDiversity(gaTrue);
  ga.pMutation((float)((float)1/(float)GENOME_SIZE));
  ga.pMutation(0.10);
  ga.pCrossover(0.0);

  ga.scoreFilename("bog.dat");
  ga.selectScores(GAStatistics::AllScores);
  ga.flushFrequency(1);
  ga.scoreFrequency(1);

  GANoScaling noscale;
  ga.scaling(noscale);

  GARouletteWheelSelector selector;
  ga.selector(selector);

  ga.parameters(argc, argv);

  time_now = time(NULL);
  cout << "Init Time: " << ctime(&time_now) << endl;
  cout << "Initializing" << endl;;
  ga.initialize(seed);
  cout << "Evolving" << endl;
  time_now = time(NULL);
  cout << "Evolve Start Time: " << ctime(&time_now) << endl;
  
  while((!ga.done()) && (stop_evolving == FALSE))
  {
    supreme_fit_count = 0;
    vhigh_fit_count = 0;
    high_fit_count = 0;
    med_fit_count = 0;
    low_fit_count = 0;
    vlow_fit_count = 0;

    ga.step();
    time_now = time(NULL);
    cout << "Time: " << ctime(&time_now) << endl;
    current_generation = ga.generation();
    cout << "Generation:" << current_generation << endl;
    cout << "best-so-far individual is: \n" << ga.statistics().bestIndividual() << "\n";
    cout << "best-so-far individual's score is: \n" << ga.statistics().bestIndividual().score() << "\n";
    cout << endl << ga.statistics() << endl;
    cout << endl;
    cout << "#supreme_fit_count  " << supreme_fit_count << endl;
    cout << "#vhigh_fit_count  " << vhigh_fit_count << endl;
    cout << "#high_fit_count  " << high_fit_count << endl;
    cout << "#med_fit_count  " << med_fit_count  << endl;
    cout << "#low_fit_count  " << low_fit_count << endl;
    cout << "#vlow_fit_count  " << vlow_fit_count << endl;

    cout.flush();
    
    if(ga.statistics().bestIndividual().score() >= max_objfunc_fitness)
        stop_evolving = TRUE;
    
    current_score_mean = ga.statistics().current(GAStatistics::Mean);
    current_best_fitness =  ga.statistics().current(GAStatistics::Maximum); 
    current_worst_fitness =  ga.statistics().current(GAStatistics::Minimum); 

    if (current_generation == 1)
        best_fitness_save = current_best_fitness;

    fitness_delta_since_save = current_best_fitness - best_fitness_save;
    if (fitness_delta_since_save < 0)
        fitness_delta_since_save = 0;       //fitness shouldn't decrease
    
    new_pcrossover = ga.pCrossover();
    new_pmutate = ga.pMutation();
    
    //adaptive step
    fuzzy_adapt_action = lee_fuzzy_adapt(&new_pcrossover, &new_pmutate, current_score_mean, 
                                         current_worst_fitness, current_best_fitness, fitness_delta_since_save);
    if(fuzzy_adapt_action == TRUE)
    {
        best_fitness_save = current_best_fitness;
        fitness_delta_since_save = 0;
    }

    //adjust the output parameters
#if 0    
    //ga.pCrossover(new_pcrossover);
#error NO CORSSOVER!
#endif
    ga.pMutation(new_pmutate);
  
    cout << "pMutation " << ga.pMutation()  << endl;
  
  }
  cout << endl;

  cout << "---" << endl;
  cout << "Done with generations.";
  cout << "best individual is: \n" << ga.statistics().bestIndividual() << "\n";
  cout << "\n" << ga.statistics() << "\n";
  cout << endl << endl << "EvCA done." << endl;
  cout << "--------------------------------------------------------------------" << endl;

  return 0;
}

//---------


//note the Rice fuzzy library returns zero membership values for calls with 
//vals above or below limit like this
//   ftriangle(13, 6, 12, 12);
//so make sure that the upper & lower membership functions 'span' the range of
//valid values!

//Based on Lee & Takagi's "Dynamic Control of Genetic Algorithms using Fuzzy Logic Techniques,"
//http://www.cs.berkeley.edu/~leem/
//
// Note:
//    1) Population size is kept constant, unlike their paper.
//    2) A sugeno-style output is used rather than a mandami output (center of gravity defuzzify)
//       unlike the paper

//Parameters:
//    pcross, pmut [output]
//    af: average fitness
//    wf: worst fitness
//    bf: best fitness
//    cf:  change in fitness since last control action
//
//   return:  TRUE if  pcross or pmut changed, FALSE otherwise
int  lee_fuzzy_adapt(float * pcross, float * pmut, float af, float wf, float bf, float cf)
{
    int ret = FALSE;
    int  fz_ab_ratio = 0;       //  average fit/best fit
    int  fz_wa_ratio = 0;       //  worst fit/average fit
    int  fz_changefit = 0;      //  change in fitness since last control action
    float ab_ratio = 0;
    float wa_ratio = 0;
    float temp1 = 0;
    float temp2 = 0;
    float cross_delta = 0;
    float mut_delta = 0;

    //get fuzzy values


    //these membership functions are skewed for empirical fit to EvCA numbers
    
    //average fit/best fit
    ab_ratio = (af/bf);
    temp1 = ftriangle(ab_ratio, 0.0*max_objfunc_fitness, 0.6*max_objfunc_fitness, 0.7*max_objfunc_fitness);
    fz_ab_ratio = FUZZY_LOW;
    temp2 = ftriangle(ab_ratio, 0.68*max_objfunc_fitness, 0.77*max_objfunc_fitness, 0.83*max_objfunc_fitness);
    if(temp2 > temp1)
    {
        fz_ab_ratio = FUZZY_MED;
        temp1 = temp2;
    }
    temp2 = ftriangle(ab_ratio, 0.82*max_objfunc_fitness, 0.90*max_objfunc_fitness, 1.0*max_objfunc_fitness);
    if(temp2 > temp1)
        fz_ab_ratio = FUZZY_HIGH;
    

    //  worst fit/average fit
    wa_ratio = (wf/af);
    temp1 = ftriangle(wa_ratio, 0.0*max_objfunc_fitness, 0.8*max_objfunc_fitness, 0.945*max_objfunc_fitness);
    fz_wa_ratio = FUZZY_LOW;
    temp2 = ftriangle(wa_ratio, 0.94*max_objfunc_fitness, 0.96*max_objfunc_fitness, 0.975*max_objfunc_fitness);
    if(temp2 > temp1)
    {
        fz_wa_ratio = FUZZY_MED;
        temp1 = temp2;
    }
    temp2 = ftriangle(wa_ratio, 0.97*max_objfunc_fitness, 1.0*max_objfunc_fitness, 1.0*max_objfunc_fitness);
    if(temp2 > temp1)
        fz_wa_ratio = FUZZY_HIGH;

    //  change in fitness since last control action
    //  note:  fitness is in [0,1]
    temp1 = ftriangle(cf, 0.0*max_objfunc_fitness, 0.0*max_objfunc_fitness, 0.03*max_objfunc_fitness);
    fz_changefit = FUZZY_LOW;
    temp2 = ftriangle(cf, 0.028*max_objfunc_fitness, 0.05*max_objfunc_fitness, 0.08*max_objfunc_fitness);
    if(temp2 > temp1)
    {
        fz_changefit = FUZZY_MED;
        temp1 = temp2;
    }
    temp2 = ftriangle(cf, 0.075*max_objfunc_fitness, 0.14*max_objfunc_fitness, 0.14*max_objfunc_fitness);
    if(temp2 > temp1)
        fz_changefit = FUZZY_HIGH;

        
    //fuzzy rule base

    //crossover rules
/* 1*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_LOW))
            cross_delta = CROSS_DELTA_LOW;
/* 2*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_MED;
/* 3*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_HIGH))
            cross_delta = CROSS_DELTA_LOW;
/* 4*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_HIGH;
/* 5*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_LOW))
            cross_delta = CROSS_DELTA_MED;
/* 6*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_LOW;
/* 7*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_HIGH))
            cross_delta = CROSS_DELTA_LOW;
/* 8*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_LOW))
            cross_delta = CROSS_DELTA_MED;
/* 9*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_MED;
/*10*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_HIGH;
/*11*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_HIGH))
            cross_delta = CROSS_DELTA_LOW;
/*12*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_LOW))
            cross_delta = CROSS_DELTA_LOW;
/*13*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_LOW;
/*14*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_HIGH))
            cross_delta = CROSS_DELTA_HIGH;
/*15*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_MED))
            cross_delta = CROSS_DELTA_LOW;
/*16*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_HIGH))
            cross_delta = CROSS_DELTA_LOW;
/*17*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_LOW))
            cross_delta = CROSS_DELTA_MED;
/*18*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_HIGH))
            cross_delta = CROSS_DELTA_MED;

    //mut rules
/* 1*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_LOW))
            mut_delta = MUT_DELTA_LOW;
/* 2*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_MED;
/* 3*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_MED))
            mut_delta = MUT_DELTA_MED;
/* 4*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_HIGH;
/* 5*/  if ((fz_ab_ratio == FUZZY_LOW) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_LOW))
            mut_delta = MUT_DELTA_MED;
/* 6*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_LOW))
            mut_delta = MUT_DELTA_MED;
/* 7*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_MED))
            mut_delta = MUT_DELTA_HIGH;
/* 8*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_MED;
/* 9*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_LOW))
            mut_delta = MUT_DELTA_MED;
/*10*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_LOW))
            mut_delta = MUT_DELTA_HIGH;
/*11*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_MED))
            mut_delta = MUT_DELTA_MED;
/*12*/  if ((fz_ab_ratio == FUZZY_MED) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_MED;
/*13*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_LOW))
            mut_delta = MUT_DELTA_LOW;
/*14*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_LOW) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_HIGH;
/*15*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_MED) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_HIGH;
/*16*/  if ((fz_ab_ratio == FUZZY_HIGH) && (fz_wa_ratio == FUZZY_HIGH) && (fz_changefit == FUZZY_HIGH))
            mut_delta = MUT_DELTA_LOW;
            
            
    //defuzzify output parameters (sugeno style)
    // use a quadratic 'beta' function on the xxx_deltas

    if((mut_delta == 0) && (cross_delta == 0))
    {
        ret = FALSE;
    }
    else
    {
        if(mut_delta != 0)
        {
            *pmut = *pmut * sqrt(sqrt((double)mut_delta));

            if(*pmut > MUT_UPPER)
                *pmut = MUT_UPPER;
            else if(*pmut < MUT_LOWER)
                *pmut = MUT_LOWER;
        }
            
        if(cross_delta != 0)
        {
            *pcross = *pcross * sqrt(sqrt((double)cross_delta));

            if(*pcross > CROSS_UPPER)
                *pcross = CROSS_UPPER;
            else if(*pcross < CROSS_LOWER)
                *pcross = CROSS_LOWER;
        }

         ret = TRUE;
    }


    cout << "ab_ratio:" << ab_ratio << "  [fz_ab_ratio:" << fz_ab_ratio << "]" << endl;
    cout << "wa_ratio:" << wa_ratio << "  [fz_wa_ratio:" << fz_wa_ratio << "]" << endl;
    cout << "cf:" << cf << "  [fz_changefit:" << fz_changefit << "]" << endl;
    cout << "mut_delta:" <<  mut_delta << "  pmut:" << *pmut << endl;
    cout << "cross_delta:" <<  cross_delta << "  pcross:" << *pcross << endl;

    return(ret);
}
