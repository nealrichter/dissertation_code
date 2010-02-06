// $Header: /usr/local/cvsroot/dissertation/code/galib/iasted/expmnts/common/GASimpleGA_1p3.cxx,v 1.3 2005/05/15 20:17:34 nealr Exp $
/* ----------------------------------------------------------------------------
  gasimple.C
  mbwall 28jul94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

  Source file for the simple genetic algorithm object.
---------------------------------------------------------------------------- */
#include <GASimpleGA_1p3.h>
#include <ga/garandom.h>


GAParameterList&
GASimpleGA_1p3::registerDefaultParameters(GAParameterList& p) {
  GAGeneticAlgorithm::registerDefaultParameters(p);

  p.add(gaNelitism, gaSNelitism,
	GAParameter::BOOLEAN, &gaDefElitism);

  return p;
}

GASimpleGA_1p3::GASimpleGA_1p3(const GAGenome& c) : GAGeneticAlgorithm(c){
  oldPop = pop->clone();

  el = gaTrue;
  params.add(gaNelitism, gaSNelitism, GAParameter::BOOLEAN, &el);
}
GASimpleGA_1p3::GASimpleGA_1p3(const GAPopulation& p) : GAGeneticAlgorithm(p){
  oldPop = pop->clone();

  el = gaTrue;
  params.add(gaNelitism, gaSNelitism, GAParameter::BOOLEAN, &el);
}
GASimpleGA_1p3::GASimpleGA_1p3(const GASimpleGA_1p3& ga) : GAGeneticAlgorithm(ga){
  oldPop = (GAPopulation *)0;
  copy(ga);
}
GASimpleGA_1p3::~GASimpleGA_1p3(){
  delete oldPop;
}
GASimpleGA_1p3&
GASimpleGA_1p3::operator=(const GASimpleGA_1p3& ga){
  if(&ga != this) copy(ga); 
  return *this;
}
void 
GASimpleGA_1p3::copy(const GAGeneticAlgorithm & g){
  GAGeneticAlgorithm::copy(g);
  const GASimpleGA_1p3& ga = DYN_CAST(const GASimpleGA_1p3&,g);
  el = ga.el;
  if(oldPop) oldPop->copy(*(ga.oldPop));
  else oldPop = ga.oldPop->clone();
  oldPop->geneticAlgorithm(*this);
}


int
GASimpleGA_1p3::setptr(const char* name, const void* value){
  int status = GAGeneticAlgorithm::setptr(name, value);

  if(strcmp(name, gaNelitism) == 0 ||
     strcmp(name, gaSNelitism) == 0){
    el = (*((int*)value) != 0 ? gaTrue : gaFalse);
    status = 0;
  }
  return status;
}

int
GASimpleGA_1p3::get(const char* name, void* value) const {
  int status = GAGeneticAlgorithm::get(name, value);

  if(strcmp(name, gaNelitism) == 0 || 
     strcmp(name, gaSNelitism) == 0){
    *((int*)value) = (el == gaTrue ? 1 : 0);
    status = 0;
  }
  return status;
}

void 
GASimpleGA_1p3::objectiveFunction(GAGenome::Evaluator f){
  GAGeneticAlgorithm::objectiveFunction(f);
  for(int i=0; i<pop->size(); i++)
    oldPop->individual(i).evaluator(f);
}

void 
GASimpleGA_1p3::objectiveData(const GAEvalData& v){
  GAGeneticAlgorithm::objectiveData(v);
  for(int i=0; i<pop->size(); i++)
    pop->individual(i).evalData(v);
}

const GAPopulation&
GASimpleGA_1p3::population(const GAPopulation& p) {
  if(p.size() < 1) {
    GAErr(GA_LOC, className(), "population", gaErrNoIndividuals);
    return *pop;
  }

  GAGeneticAlgorithm::population(p);
  oldPop->copy(*pop->clone());
  oldPop->geneticAlgorithm(*this);

  return *pop;
}

int 
GASimpleGA_1p3::populationSize(unsigned int n) {
  GAGeneticAlgorithm::populationSize(n);
  oldPop->size(n);
  return n;
}

int 
GASimpleGA_1p3::minimaxi(int m) { 
  GAGeneticAlgorithm::minimaxi(m);
  if(m == MINIMIZE)
    oldPop->order(GAPopulation::LOW_IS_BEST);
  else
    oldPop->order(GAPopulation::HIGH_IS_BEST);
  return minmax;
}










// Initialize the population, set the random seed as needed, do a few stupidity
// checks, reset the stats.  We must initialize the old pop because there is no
// guarantee that each individual will get initialized during the course of our
// operator++ operations.  We do not evaluate the old pop because that will 
// happen as-needed later on.
void
GASimpleGA_1p3::initialize(unsigned int seed)
{
  GARandomSeed(seed);

  pop->initialize();
  pop->evaluate(gaTrue);	// the old pop will get it when the pops switch
//  oldPop->initialize();

  stats.reset(*pop);

  if(!scross) 
    GAErr(GA_LOC, className(), "initialize", gaErrNoSexualMating);
}


//   Evolve a new generation of genomes.  When we start this routine, pop
// contains the current generation.  When we finish, pop contains the new 
// generation and oldPop contains the (no longer) current generation.  The 
// previous old generation is lost.  We don't deallocate any memory, we just
// reset the contents of the genomes.
//   The selection routine must return a pointer to a genome from the old
// population.
void
GASimpleGA_1p3::step()
{
  int i, mut;
  float child_score, mom_score;
  GAGenome *mom;          // tmp holders for selected genomes

  GAPopulation *tmppop;		// Swap the old population with the new pop.
  tmppop = oldPop;		// When we finish the ++ we want the newly 
  oldPop = pop;			// generated population to be current (for
  pop = tmppop;			// references to it from member functions).
  int max_index = 0;

// Generate the individuals in the temporary population from individuals in 
// the main population.

  for(i=0; i<pop->size(); i++)
  {
      mom = &(oldPop->select());  
      stats.numsel += 1;		// keep track of number of selections

      //get mom score
      mom_score = mom->score();

      GAPopulation * sub_pop = new GAPopulation(*mom,4);
      int max_index = 0;

      sub_pop->size(4);
      sub_pop->initialize();
      sub_pop->evaluate(gaFalse);

      //copy mom
      sub_pop->individual( 0 ).copy(*mom);
      sub_pop->individual( 0 ).evaluator(pop->individual( 0 ).evaluator());

      //produce 3 children
      for(int j = 1; j < 4; j++)
      {
          sub_pop->individual( j ).copy(*mom);
          sub_pop->individual( j ).evaluator(pop->individual( 0 ).evaluator());
          stats.nummut += (mut = sub_pop->individual( j ).mutate(pMutation()));

          //get child score
          child_score = sub_pop->individual( j ).score();
          stats.numeval++;

          //find max
          if(max_index != j)
          {
              if(child_score >= sub_pop->individual( max_index ).score())
                  max_index = j;
          }

      }//end of production of 3 children


      //put max in new population
      pop->individual( i ).copy(sub_pop->individual( max_index ));
  }
  
  stats.numrep += pop->size();

  pop->evaluate(gaTrue);    // get info about current pop for next time

// If we are supposed to be elitist, carry the best individual from the old
// population into the current population.  Be sure to check whether we are
// supposed to minimize or maximize.

  if(minimaxi() == GAGeneticAlgorithm::MAXIMIZE) {
    if(el && oldPop->best().score() > pop->best().score())
      oldPop->replace(pop->replace(&(oldPop->best()), GAPopulation::WORST), 
		      GAPopulation::BEST);
  }
  else {
    if(el && oldPop->best().score() < pop->best().score())
      oldPop->replace(pop->replace(&(oldPop->best()), GAPopulation::WORST), 
		      GAPopulation::BEST);
  }

  stats.update(*pop);		// update the statistics by one generation
}
