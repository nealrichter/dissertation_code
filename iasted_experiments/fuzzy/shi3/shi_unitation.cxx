/* ----------------------------------------------------------------------------
  shi_unitation.cxx

  Neal Richter
  Montana State University
    
  Copyright (c) 2005 Neal Richter

 DESCRIPTION:
   1+3 unitation GA with shi fuzzy adaptive rules

   
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

#include "royalroad_header.h"

#define GENOME_SIZE   20
#define POPULATION_SIZE   10
#define NUM_GEN 1000

int evals;
static int obj_case = 0;
static int global_population_size = POPULATION_SIZE;
static int global_genome_size = GENOME_SIZE;
static int global_num_gen = NUM_GEN;

float max_objfunc_fitness = 0;

#include "unitation_header.h"
#include "longpath_header.h"

#define LATE_GEN_THRESH     .6

#define NUM_SAMPLES_MAX	       2500
#define NUM_SAMPLES_MIN	       250
#define NUM_SAMPLES_DEFAULT	   500

#define FUZZY_HIGH           3
#define FUZZY_MED            2
#define FUZZY_LOW            1

#define MUT_LOW             (float)(1/(float)global_genome_size)
#define MUT_MED             (float)(2/(float)global_genome_size)
#define MUT_HIGH            (float)(4/(float)global_genome_size)

#define CROSS_LOW           0.95
#define CROSS_MED           0.97
#define CROSS_HIGH          1.00


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
  
unsigned int seed = 0;

//functions
void  shi_fuzzy_adapt(float *, float *, float, float, int);

int
main(int argc, char** argv) 
{
  time_t time_now;
  struct timeb mytime;
  int stop_evolving = FALSE;

  int    unchanged_best_fitness_count = 0;
  float  new_pcrossover = 0;
  float  new_pmutate = 0;
  float  old_best_fitness = 0;
  float  current_best_fitness = 0;
  float  current_score_variance = 0;
  
  if((argc != 2) && (argc != 3))
  {
    cout << "usage: "<< argv[0] << "objfunc#  [popsize]" << endl;
    cout << "\t\tobjfuncs:\n";
    cout << "\t\t\t:1 = ONEMAX\n";
    cout << "\t\t\t:2 = NEEDLE\n";
    cout << "\t\t\t:3 = DECTRAP\n";
    cout << "\t\t\t:4 = DECTWOTRAP\n";
    cout << "\t\t\t:5 = TWOTRAP\n";
    cout << "\t\t\t:6 = Royal Road\n";
    cout << "\t\t\t:7 = LongPath\n";
    exit(1);
  }
  
  obj_case = atoi(argv[1]);
  if((argc == 3) && (strlen(argv[2]) > 0))
    global_population_size = atoi(argv[2]);
  
  if((obj_case >= 1) && (obj_case <= 5))
  {
      cout << "Unitation" << global_genome_size << " bits, mutation only\n\n";
      cout << "GASimpleGA_1p3 (non-overlapping populations)";
  }
  else if(obj_case == 6)
  {
      global_genome_size = royalroad_nbits;
      cout << "Royal Road " << global_genome_size << " bits, mutation only\n\n";
      cout << "GASimpleGA_1p3 (non-overlapping populations)";
  }
  else if(obj_case == 7)
  {
      global_genome_size = LONGPATH_BITS;
      cout << "LongPath " << global_genome_size << " bits, mutation only\n\n";
      cout << "GASimpleGA_1p3 (non-overlapping populations) - MAXIMIZER";
      cout << endl;
  }
  
  cout << "Shi Fuzzy Adaptive Rules for mutation rate control\n";
  cout << endl;
  
  cout << " MUT_LOW = " << MUT_LOW;
  cout << " MUT_MED = " << MUT_MED;
  cout << " MUT_HIGH = " << MUT_HIGH << endl;
  cout << " CROSS_LOW = " << CROSS_LOW;
  cout << " CROSS_MED = " << CROSS_MED;
  cout << " CROSS_HIGH = " << CROSS_HIGH << endl;
  cout << endl;
  
  //default random seed
  ftime(&mytime);
  seed = (unsigned int)mytime.millitm;
  srand(seed);
  
  GA1DBinaryStringGenome worst(global_genome_size, Objective);

  switch(obj_case)
  {

      case 1:// ONEMAX
          cout << "ONEMAX Unitation Obj" << endl;
          max_objfunc_fitness = global_genome_size;
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          break;
      case 2:// NEEDLE
          cout << "NEEDLE Unitation Obj" << endl;
          max_objfunc_fitness = global_genome_size;
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          worst.set(0,1);
          break;
      case 3:// DECTRAP
          cout << "DECTRAP Unitation Obj" << endl;
          max_objfunc_fitness = global_genome_size;
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          worst.set(0,2);
          break;
      case 4:// DECTWOTRAP
          cout << "DECTWOTRAP Unitation Obj" << endl;
          max_objfunc_fitness = global_genome_size;
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          worst.set(0,4);
          break;
      case 5:// TWOTRAP
          cout << "TWOTRAP Unitation Obj" << endl;
          max_objfunc_fitness = global_genome_size;
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          worst.set(0,8);
          break;
      case 6:// Royal Road
          cout << "Royal Road Obj" << endl;
            
          { 
            GA1DBinaryStringGenome best(global_genome_size, Objective);
            best.set(0,global_genome_size);
            max_objfunc_fitness = Objective(best);
            cout << "Royal Road max is " << max_objfunc_fitness << endl;
            royalroad_highestLevel=0;
          }
          
          global_num_gen = RR_NUM_GEN;
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          break;

      case 7:// LongPath
          cout << "LongPath" << endl;
            
          max_objfunc_fitness = longpath_best_obj_val;
          cout << "LongPath best (min) is " << max_objfunc_fitness << endl;
          
          worst.set(0,global_genome_size);
          worst.unset(0,global_genome_size);
          break;

      default:
          cout << "NO Obj func defined" << endl;
          exit(1);

  }

  GA1DBinaryStringGenome genome(global_genome_size, Objective);
  GASimpleGA_1p3 ga(genome);

  ga.elitist(gaFalse);
  cout << "No Elitism." << endl;

  ga.populationSize(global_population_size);
  ga.nGenerations(global_num_gen);
  cout << "PopulationSize = " <<  ga.populationSize(); 
  cout << " Num Generations = " << ga.nGenerations() << endl;
  
  ga.recordDiversity(gaTrue);
  ga.pMutation((float)((float)1/(float)global_genome_size));
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
  
  GAPopulation temp = ga.population();
  temp.remove(0);
  temp.add(worst);
  ga.population(temp);

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
    
    if(ga.generation() == 0)
    {
        temp.remove(0);
        temp.add(worst);
        ga.pMutation(0.0);
        temp.evaluate(gaTrue);
        ga.population(temp);
        ga.step();
        cout << "NEAL worst is: \n" << ga.population().best() << "\n";
        ga.pMutation((float)((float)1/(float)global_genome_size));
    }
    else
        ga.step();
  
    time_now = time(NULL);
    cout << "Time: " << ctime(&time_now) << endl;
    current_generation = ga.generation();
    cout << "Generation:" << current_generation << endl;
    cout << "best-so-far individual is: \n" << ga.statistics().bestIndividual() << "\n";
    cout << "best-so-far individual's score is: \n" << ga.statistics().bestIndividual().score() << "\n";
    
    if(obj_case == 6)
    {
        if((ga.generation() % 50) == 0)
            cout << endl << ga.statistics() << endl;
    }
    else
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
    current_score_stddev = ga.statistics().current(GAStatistics::Deviation);
    current_best_fitness =  ga.statistics().current(GAStatistics::Maximum); 

    //get adaptive step parameters ready

    if (current_best_fitness == old_best_fitness)
        unchanged_best_fitness_count++;
    else
        unchanged_best_fitness_count = 0;

    old_best_fitness = current_best_fitness;
    
    new_pcrossover = ga.pCrossover();
    new_pmutate = ga.pMutation();

    //scale & calculate variance
    current_score_variance =  current_score_stddev*100;
    current_score_variance *=  current_score_variance;
    current_score_variance *= (float)1/100;
	    

    
    //adaptive step
    shi_fuzzy_adapt(&new_pcrossover, &new_pmutate, current_best_fitness, 
                     current_score_variance, unchanged_best_fitness_count);

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
  
  if(obj_case == 6)
  {
    cout << "the highest level achieved was " << royalroad_highestLevel << "\n";
  }

  cout << "\n" << ga.statistics() << "\n";
  cout << endl << endl << "EvCA done." << endl;
  cout << "--------------------------------------------------------------------" << endl;

  return 0;
}

float Objective(GAGenome &genome)
{

  if(obj_case == 6)
  {
      return( RoyalRoad_Objective(genome) );
  }
  else if(obj_case == 7)
  {
      return( LongPath_Objective(genome) );
  }
  else
  {
      return( Unitation_Objective(genome) );
  }

}

//---------

//note the Rice fuzzy library returns zero membership values for calls with 
//vals above or below limit like this
//   ftriangle(13, 6, 12, 12);
//so make sure that the upper & lower membership functions 'span' the range of
//valid values!
//
//  NOTE: 5/9/2002 Added 'late gen = high mutation rule"
//
void  shi_fuzzy_adapt(float * pcross, float * pmut, float bf, float vf, int uf)
{
    int  fz_bf = 0;
    int  fz_vf = 0;
    int  fz_uf = 0;
    int  fz_mr = 0;
    int  fz_cr = 0;
    float temp1 = 0;
    float temp2 = 0;

    //get fuzzy values

    //best fitness
    //adjusted/scaled to .5 <-> 1.0
    temp1 = ftriangle(bf, 0.0*max_objfunc_fitness, 0.5*max_objfunc_fitness, 0.7*max_objfunc_fitness);
    fz_bf = FUZZY_LOW;
    temp2 = ftriangle(bf, 0.65*max_objfunc_fitness, 0.775*max_objfunc_fitness, 0.9*max_objfunc_fitness);
    if(temp2 > temp1)
    {
        fz_bf = FUZZY_MED;
        temp1 = temp2;
    }
    temp2 = ftriangle(bf, 0.85*max_objfunc_fitness, 1.0*max_objfunc_fitness, 1.0*max_objfunc_fitness);
    if(temp2 > temp1)
        fz_bf = FUZZY_HIGH;
    
    
    //variance of fitness - 
    //adjusted for empirical values observed
    temp1 = ftriangle(vf, 0.0*max_objfunc_fitness, 0.0*max_objfunc_fitness, 0.12*max_objfunc_fitness);
    fz_vf = FUZZY_LOW;
    temp2 = ftriangle(vf, 0.1*max_objfunc_fitness, 0.16*max_objfunc_fitness, 0.22*max_objfunc_fitness);
    if(temp2 > temp1)
    {
        fz_vf = FUZZY_MED;
        temp1 = temp2;
    }
    temp2 = ftriangle(vf, 0.20*max_objfunc_fitness, .32*max_objfunc_fitness, 1*max_objfunc_fitness);
    if(temp2 > temp1)
        fz_vf = FUZZY_HIGH;

    //unchanged fitness
    temp1 = ftriangle(uf, 0.0, 0.0, 6);
    fz_uf = FUZZY_LOW;
    temp2 = ftriangle(uf, 3, 6, 9);
    if(temp2 > temp1)
    {
        fz_uf = FUZZY_MED;
        temp1 = temp2;
    }
    temp2 = ftriangle(uf, 6, 12, global_num_gen);
    if(temp2 > temp1)
        fz_uf = FUZZY_HIGH;

        
    //fuzzy rule base
    if (fz_bf == FUZZY_LOW)
    {
        fz_mr = FUZZY_LOW;
        fz_cr = FUZZY_HIGH;
    }
    if ((fz_bf == FUZZY_MED) && (fz_uf == FUZZY_LOW))
    {
        fz_mr = FUZZY_LOW;
        fz_cr = FUZZY_HIGH;
    }
    if ((fz_bf == FUZZY_MED) && (fz_uf == FUZZY_MED))
    {
        fz_mr = FUZZY_MED;
        fz_cr = FUZZY_MED;
    }
    if ((fz_bf == FUZZY_HIGH) && (fz_uf == FUZZY_LOW))
    {
        fz_mr = FUZZY_LOW;
        fz_cr = FUZZY_HIGH;
    }
    if ((fz_bf == FUZZY_HIGH) && (fz_uf == FUZZY_MED))
    {
        fz_mr = FUZZY_MED;
        fz_cr = FUZZY_MED;
    }

    if((vf != 0) && (global_population_size > 10))
    {
        if ((fz_uf == FUZZY_HIGH) && (fz_vf == FUZZY_LOW))
        {
            fz_mr = FUZZY_HIGH;
            fz_cr = FUZZY_LOW;
        }
        if ((fz_uf == FUZZY_HIGH) && (fz_vf == FUZZY_MED))
        {
            fz_mr = FUZZY_HIGH;
            fz_cr = FUZZY_LOW;
        }
        if ((fz_uf == FUZZY_HIGH) && (fz_vf == FUZZY_HIGH))
        {
            fz_mr = FUZZY_LOW;
            fz_cr = FUZZY_LOW;
        }
    }
    else
    {
        fz_cr = FUZZY_LOW;

        //choose mut randomly (1,2,3)
        int rnum = 1+(int)(3.0*rand()/(RAND_MAX+1.0));

        fz_mr = rnum;
        
    }

    //late generation rule  to reset the MUT rate to be higher
    //if((current_generation >= (int)(LATE_GEN_THRESH*global_population_size)) && (fz_mr == FUZZY_LOW))
    //    fz_mr = FUZZY_HIGH;

    //defuzzify output parameters

    //mutation rate
    //altered to pick MED as 0.03 (original EvCA value)
    //not symetric!
    //jan/25/  aggressive settings
    if (fz_mr == FUZZY_LOW)
        *pmut = MUT_LOW;
    else if (fz_mr == FUZZY_MED)
        *pmut = MUT_MED;
    else if (fz_mr = FUZZY_HIGH)
        *pmut = MUT_HIGH;


    //crossover rate
    //altered to pick HIGH as 1.00 (original EvCA value)
    //not symetric! and more conservative than SHI original values
    //1/25  less corssover variation.
    if (fz_cr == FUZZY_LOW)
        *pcross = CROSS_LOW;
    else if (fz_cr == FUZZY_MED)
        *pcross = CROSS_MED;
    else if (fz_cr == FUZZY_HIGH)
        *pcross = CROSS_HIGH;


    cout << "bf:" << bf << "  [fz_bf:" << fz_bf << "]" << endl;
    cout << "vf:" << vf << "  [fz_vf:" << fz_vf << "]" << endl;
    cout << "uf:" << uf << "  [fz_uf:" << fz_uf << "]" << endl;
    cout << "fz_mr:" <<  fz_mr << "  pmut:" << *pmut << endl;
    cout << "fz_cr:" <<  fz_cr << "  pcross:" << *pcross << endl;

    
}
