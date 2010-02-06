/* ----------------------------------------------------------------------------
  adaptive_mu_unitation.cxx

  Neal Richter
  Montana State University
    
  Copyright (c) 2005 Neal Richter

 DESCRIPTION:
   unitation, adaptive mutrate - mutation only

   droste mu update rule/cycle
   back mu update
   
   Uses GAlib library and example source as a guide/template
   Copyright (c) 1995-199 Massachusetts Institute of Technology and Matthew Wall.
   http://lancet.mit.edu/ga/
   Thanks MATT!
---------------------------------------------------------------------------- */
#include <stdio.h>
#include <sys/timeb.h>
#include <iostream.h>
#include <math.h>
#include <ga/GA1DBinStrGenome.h>

#include "GASimpleGA_1p3_adaptive.h"

#include "royalroad_header.h"

#define GENOME_SIZE   20
#define POPULATION_SIZE   1
#define NUM_GEN 1000

#define MUT_RATE_LOWER_LIMIT    (float)((float)1/(float)global_genome_size)
#define MUT_RATE_UPPER_LIMIT    0.5


int evals;

static int obj_case = 0;
static int global_population_size = POPULATION_SIZE;
static int global_genome_size = GENOME_SIZE;
static int global_num_gen = NUM_GEN;

float max_objfunc_fitness = 0;

#include "unitation_header.h"
#include "longpath_header.h"

int
main(int argc, char** argv) 
{
  struct timeb tp;
  unsigned int seed = 0;
  int stop_evolving = FALSE;
  
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
      cout << endl;
  }
  else if(obj_case == 6)
  {
      global_genome_size = royalroad_nbits;
      cout << "Royal Road " << global_genome_size << " bits, mutation only\n\n";
      cout << "GASimpleGA_1p3 (non-overlapping populations)";
      cout << endl;
  }
  else if(obj_case == 7)
  {
      global_genome_size = LONGPATH_BITS;
      cout << "LongPath " << global_genome_size << " bits, mutation only\n\n";
      cout << "GASimpleGA_1p3 (non-overlapping populations) - MAXIMIZER";
      cout << endl;
  }


#ifdef RECHENBERG_ADAPTIVE
  cout << "Rechenberg 1/5 Rule Adaptive Update of Mutation Rate" << endl;
#elif THIERENS_CONSTGAIN_ADAPTIVE
  cout << "Thierens Constant Gain Adpative Update of Mutation Rate" << endl;
#elif THIERENS_DECLINING_ADAPTIVE
  cout << "Thierens Declining Adpative Update of Mutation Rate" << endl;
#else
#error ACK!
  cout << "Static Mutation Rate" << endl;
#endif

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


  //default random seed

  ftime(&tp);
  seed = tp.millitm;

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
  ga.pMutation(MUT_RATE_LOWER_LIMIT);
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
  
  cout << "initializing..."; cout.flush();
  ga.initialize(seed);
  
  //chose adaptive method
#ifdef RECHENBERG_ADAPTIVE
  ga.adaptiveMethod(RECHENBERG_ADAPTIVE_METHOD);
#elif THIERENS_CONSTGAIN_ADAPTIVE
  ga.adaptiveMethod(THIERENS_CONSTGAIN_ADAPTIVE_METHOD);
#elif THIERENS_DECLINING_ADAPTIVE
  ga.adaptiveMethod(THIERENS_DECLINING_ADAPTIVE_METHOD);
#endif
  
  reset_distribution();
  cout << "evolving"; cout.flush();
  while((!ga.done()) && (stop_evolving == FALSE))
  {

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
  
    cout << endl;
    cout << "Generation:" << ga.generation() << endl;
    output_distribution(ga);
    reset_distribution();
    cout << "best-so-far individual is: \n" << ga.statistics().bestIndividual() << "\n";
    cout << "best-so-far individual's score is: \n" << ga.statistics().bestIndividual().score() << "\n";

    if(obj_case == 6)
    {
        if((ga.generation() % 50) == 0)
            cout << endl << ga.statistics() << endl;
    }
    else
        cout << endl << ga.statistics() << endl;
    
    cout.flush();

    if(ga.statistics().bestIndividual().score() >= max_objfunc_fitness)
        stop_evolving = TRUE;

    //update mutation rate via rule
    {
#ifdef RECHENBERG_ADAPTIVE
        if((ga.generation() % 5) == 0)
            ga.updateMutation(MUT_RATE_LOWER_LIMIT, MUT_RATE_UPPER_LIMIT);
#else
        ga.updateMutation(MUT_RATE_LOWER_LIMIT, MUT_RATE_UPPER_LIMIT);
#endif

        cout << endl << "New mutation rate is " <<  ga.pMutation() << endl;
    }
  }
  cout << endl;

  cout << "Done with generations.";
  cout << "best individual is: \n" << ga.statistics().bestIndividual() << "\n";
  
  if(obj_case == 6)
  {
    cout << "the highest level achieved was " << royalroad_highestLevel << "\n";
  }

  cout << "\n" << ga.statistics() << "\n";
  cout << endl << endl << "EvCA done." << endl;
  cout << "----------------------------------------------------------------" << endl;

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
