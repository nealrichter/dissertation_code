#ifndef _ROYALROAD_HEADER_H_
#define _ROYALROAD_HEADER_H_

//**************************************************

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

float Objective(GAGenome &);

/* ----------------------------------------------------------------------------
  ex20.C
  mbwall 5sep95
  Copyright (c) 1995-1996  Massachusetts Institute of Technology

 DESCRIPTION:
   This example runs the royal road problem.  See the comments near the 
objective functions for details about the function itself.
   Some of this was copied (at least partially) from the galopps genetic
algorithm library and from the pga package.  I used a bunch of globals in this
example - not good programming style, but it gets the job done.
---------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <ga/ga.h>


// This is the objective function for computing Holland's 1993 ICGA version
// of the Royal Road problem.  It has been corrected per GAList volume 7
// number 23, 8/26/93.  No bonus points are awarded for a given level until 
// it has been achieved (this fixes Holland's coding error in GAList).
//   Holland posed this problem as a challenge to test the 
// performance of genetic algorithms.  He indicated that, with the parameter 
// settings of
//
//     schemata size = 8
//     bits between schemata = 7
//     m* = 4
//     U* = 1.0
//     u = 0.3
//     v = 0.02
//
// he could attain royal_road_level 3 most of the time within
// 10,000 function evaluations.  He challenged other GA users to match or beat
// that performance.  He indicated that he used a population size of 512 to
// obtain his solutions, and did NOT use a "simple genetic algorithm."
//   The genome for this problem is a single-dimension bit string with length
// defined by the block size and gap size as:
//
//     length = (blocksize+gapsize) * (2^K)
//
// where K= 1,2,3, or 4.  Holland used K = 4.

#define RR_NUM_GEN  20000

#define RR_NBLOCKS 16		// this number is 2^K

const int RR_BLOCKSIZE= 4; //8;	// block size - length of target schemata
const int RR_GAPSIZE=3; //7;	// gap size - number of bits between target schemata
const int RR_MSTAR=2; //4;	// Holland's m* - up to this many bits in low level
			// block gets reward
const float RR_USTAR=1.0;     // Holland's U* - first block earns this
const float RR_U=0.3;      // Holland's u - increment for lowest level match
const float RR_V=0.02;     // Holland's v - reward/penalty per bit

int royalroad_nbits = (RR_BLOCKSIZE+RR_GAPSIZE)*RR_NBLOCKS;
int royalroad_blockarray[RR_NBLOCKS];
int royalroad_highestLevel=0;

float
RoyalRoad_Objective(GAGenome & c){
  GA1DBinaryStringGenome & genome = (GA1DBinaryStringGenome &)c;

  float score = 0.0;
  int total, i, j, index, n;

  int stop = 0;
// do the lowest level blocks first

  n = 0;
  for(i=0; i<RR_NBLOCKS; i++) {
      total = 0;
      for(j=i*(RR_BLOCKSIZE + RR_GAPSIZE); j<i*(RR_BLOCKSIZE+RR_GAPSIZE)+RR_BLOCKSIZE; j++)
          if(genome.gene(j) == 1) total++;  // count the bits in the block
      if(total > RR_MSTAR && total < RR_BLOCKSIZE)
          score -= (total-RR_MSTAR)*RR_V;
      else if(total <= RR_MSTAR)
          score += total * RR_V; 
      if(total == RR_BLOCKSIZE) {
          royalroad_blockarray[i] = 1;
          n++;
      }
      else{
          royalroad_blockarray[i] = 0;
      }
  }

// bonus for filled low-level blocks

  if(n > 0) score += RR_USTAR + (n-1)*RR_U;

// now do the higher-level blocks

  n = RR_NBLOCKS;		// n is now number of filled low level blocks
  int proceed = 1;		// should we look at the next higher level?
  int level = 0;
  while ((n > 1) && proceed) {
      proceed = 0;
      total = 0;
      /* there are n valid blocks in the royalroad_blockarray each time */
      /* round, so n=2 is the last.                           */
      for(i=0,index=0; i<(n/2)*2; i+=2,index++) {
          if(royalroad_blockarray[i] == 1 && royalroad_blockarray[i+1] == 1) {
              total++;
              proceed = 1;
              royalroad_blockarray[index] = 1;
          }
          else{
              royalroad_blockarray[index] = 0;
          }
      }
      if(total > 0){
          score += RR_USTAR + (total-1)*RR_U;
          level++;
      }
      n /= 2;
  }

  if(royalroad_highestLevel < level) royalroad_highestLevel = level;
 
  return(score*10);
}








//**************************************************


#endif /* _ROYALROAD_HEADER_H_ */
