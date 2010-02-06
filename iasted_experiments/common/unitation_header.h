#ifndef _UNITATION_HEADER_H_
#define _UNITATION_HEADER_H_

//**************************************************

#define TRUE  1
#define FALSE 0

float Unitation_Objective(GAGenome &);
void output_distribution(GASimpleGA_1p3 &ga);
void reset_distribution();
int unitation_distribution[GENOME_SIZE+1];

#if (GENOME_SIZE == 10)
            int twotrap_fitness[GENOME_SIZE+1] = {10, 5, 0, 2, 4, 6, 4, 2, 0, 5, 10};
            int dectwotrap_fitness[GENOME_SIZE+1] = {8, 6, 4, 2, 0, 10, 0, 2, 4, 6, 8};
#elif (GENOME_SIZE == 20)
            int twotrap_fitness[GENOME_SIZE+1] = {20, 15, 10, 5, 0, 2, 4, 6, 8, 10, 12, 10, 8, 6, 4, 2, 0, 5, 10, 15, 20};
            int dectwotrap_fitness[GENOME_SIZE+1] = {16, 14, 12, 10, 8, 6, 4, 2, 1, 0, 20, 0, 1, 2, 4, 6, 8,10, 12, 14, 16};
#else
#error NO 2TRAP obj func for this genome size
#error NO DEC2TRAP obj func for this genome size
#endif


//Unitation objective func

float
Unitation_Objective(GAGenome& g) 
{
    GA1DBinaryStringGenome & genome = (GA1DBinaryStringGenome &)g;
    float ones_count=0.0001;
    float final_score = 0.0;
    
    //count 1s
    for(int i=0; i<genome.length(); i++)
      ones_count += genome.gene(i);

    switch(obj_case)
    {
        case 1: // ONEMAX

            final_score = ones_count;
            break;

        case 2: // NEEDLE

            if(ones_count == genome.length())
                final_score = ones_count;
            else
                final_score = 1;
            break;

        case 3:// DECTRAP

            if(ones_count == genome.length())
                final_score = ones_count;
            else
                final_score = genome.length() * (1 - ((1+ones_count)/( genome.length() )));
            break;

        case 4: // DECTWOTRAP

            final_score = (double) dectwotrap_fitness[(int)ones_count];
            break;

        case 5: // TWOTRAP

            final_score = (double) twotrap_fitness[(int)ones_count];
            break;


        default:
            cout << "ERROR NO Objective func defined" << endl;
            exit(1);
    }

    evals++;

    return(final_score);
}


//void output_distribution(int popsize)
void output_distribution(GASimpleGA_1p3 &ga)
{
    int popsize = ga.populationSize();
    int idx = 0;
    GAPopulation pop = ga.population();
    
    cout << "Evals = " << evals << endl;

    cout << "Scaled Population distribution (unitation distribution)" << endl;

    for(int i=0; i< popsize; i++)
    {
        idx = (int) pop.individual(i).score();
        unitation_distribution[idx]++;
    }

    for (int i=0; i<(GENOME_SIZE+1); i++)
        cout << (float)((float)unitation_distribution[i]/popsize)  << ", ";

    cout << endl;

}


void reset_distribution()
{ 
    for (int i=0; i<(GENOME_SIZE+1); i++)
       unitation_distribution[i] = 0;

    evals = 0;
}


//**************************************************


#endif /* _UNITATION_HEADER_H_ */
