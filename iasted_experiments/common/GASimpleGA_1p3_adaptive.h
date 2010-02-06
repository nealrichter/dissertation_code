// $Header: /usr/local/cvsroot/dissertation/code/galib/iasted/expmnts/common/GASimpleGA_1p3_adaptive.h,v 1.1 2005/03/28 07:00:20 nealr Exp $
/* ----------------------------------------------------------------------------
  gasimple.h
  mbwall 28jul94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

  Header file for the simple genetic algorithm class.
---------------------------------------------------------------------------- */
#ifndef _ga_gasimple_h_
#define _ga_gasimple_h_

#include <ga/GABaseGA.h>

#define RECHENBERG_ADAPTIVE_METHOD           1
#define THIERENS_CONSTGAIN_ADAPTIVE_METHOD   2
#define THIERENS_DECLINING_ADAPTIVE_METHOD   3

//#define RECHENBERG_ALPHA                    (float) 1.1
#define RECHENBERG_ALPHA                    (float) 2.0
#define THIERENS_CONSTGAIN_LAMBDA           (float) 1.1
#define THIERENS_CONSTGAIN_OMEGA            (float) 1.5

#define THIERENS_DECADAPT_LAMBDA            (float) 2.0
#define THIERENS_DECADAPT_OMEGA             (float) 2.0
#define THIERENS_DECADAPT_GAMMA             (float) 0.95

class GASimpleGA_1p3 : public GAGeneticAlgorithm {
public:
  GADefineIdentity("GASimpleGA_1p3", GAID::SimpleGA);

  static GAParameterList& registerDefaultParameters(GAParameterList&);

public:
  GASimpleGA_1p3(const GAGenome&);
  GASimpleGA_1p3(const GAPopulation&);
  GASimpleGA_1p3(const GASimpleGA_1p3&);
  GASimpleGA_1p3& operator=(const GASimpleGA_1p3&);
  virtual ~GASimpleGA_1p3();
  virtual void copy(const GAGeneticAlgorithm&);

  virtual void initialize(unsigned int seed=0);
  virtual void step();
  GASimpleGA_1p3 & operator++() { step(); return *this; }

  virtual int setptr(const char* name, const void* value);
  virtual int get(const char* name, void* value) const;

  GABoolean elitist() const {return el;}
  GABoolean elitist(GABoolean flag)
    {params.set(gaNelitism, (int)flag); return el=flag;}

  virtual int minimaxi() const {return minmax;}
  virtual int minimaxi(int m);

  virtual const GAPopulation& population() const {return *pop;}
  virtual const GAPopulation& population(const GAPopulation&);
  virtual int populationSize() const {return pop->size();}
  virtual int populationSize(unsigned int n);
  virtual GAScalingScheme& scaling() const {return pop->scaling();}
  virtual GAScalingScheme& scaling(const GAScalingScheme & s)
    {oldPop->scaling(s); return GAGeneticAlgorithm::scaling(s);}
  virtual GASelectionScheme& selector() const {return pop->selector(); }
  virtual GASelectionScheme& selector(const GASelectionScheme& s)
    {oldPop->selector(s); return GAGeneticAlgorithm::selector(s);}
  virtual void objectiveFunction(GAGenome::Evaluator f);
  virtual void objectiveData(const GAEvalData& v);
  
  virtual float updateMutation(float, float);

  virtual int adaptiveMethod(int);
  virtual int adaptiveMethod(void);

protected:
  GAPopulation *oldPop;		// current and old populations
  GABoolean el;			// are we elitist?

//private:
  int adaptive_method;
  int rechenberg_success;
  int rechenberg_failure;
  int thierens_constgain_case0_success;
  int thierens_constgain_case1_success;
  int thierens_constgain_case2_success;
  int thierens_constgain_case3_success;
  int thierens_decadapt_case0_success;
  int thierens_decadapt_case1_success;
  int thierens_decadapt_case2_success;
  int thierens_decadapt_case3_success;
  float new_mu;

};



#ifndef NO_STREAMS
inline ostream& operator<< (ostream& os, GASimpleGA_1p3 & arg)
{ arg.write(os); return(os); }
inline istream& operator>> (istream& is, GASimpleGA_1p3 & arg)
{ arg.read(is); return(is); }
#endif

#endif
