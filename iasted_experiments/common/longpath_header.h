#ifndef _LONGPATH_HEADER_H_
#define _LONGPATH_HEADER_H_

//**************************************************

#define HEX__(n) 0x##n##LU

/* 8-bit conversion function */
#define B8__(x) ((x&0x0000000FLU)?1:0)  \
               +((x&0x000000F0LU)?2:0)  \
               +((x&0x00000F00LU)?4:0)  \
               +((x&0x0000F000LU)?8:0)  \
               +((x&0x000F0000LU)?16:0) \
               +((x&0x00F00000LU)?32:0) \
               +((x&0x0F000000LU)?64:0) \
               +((x&0xF0000000LU)?128:0)

/*
 * *** USER MACROS ***
 */

/* for upto 8-bit binary constants */
#define B8(d) ((unsigned char)B8__(HEX__(d)))

/* for upto 16-bit binary constants, MSB first */
#define B16(dmsb,dlsb) (((unsigned short)B8(dmsb)<<8) + B8(dlsb))

/*
 * Sample usage:
 *    B8(01010101) = 85
 *    B16(10101010,01010101) = 43605
 */
    
//**************************************************

#define TRUE  1
#define FALSE 0

#define LONGPATH_BITS    9

float LongPath_Objective(GAGenome &);

static int OnPath(GA1DBinaryStringGenome&);

static float longpath_max_val = (float) 3 * pow(2,(float)((float)LONGPATH_BITS-1)/(float)2) - 2;

float longpath_best_obj_val = longpath_max_val * 2;

static float genome_pos = -1;

//from pg 107 & 108 of Gunter Rudolph's PhD Thesis
//l = 9 , 9 bit genome!!
//
//  USE MINIMIZING GA!!!!!
//
//Code is a bit wierd.. but I coded it to follow 
//Rudolph's psuedo-code extactly.  The OnPath
//func could return genome_pos... but hey, it's
//research code.  Like that is an excuse.

int longpath_pos[47] = {  
                          B16(00000000, 00000000), \ //#0
                          B16(00000000, 00000001), \ //#1
                          B16(00000000, 00000011), \ //#2
                          B16(00000000, 00000111), \ //#3
                          B16(00000000, 00000110), \ //#4
                          B16(00000000, 00001110), \ //#5
                          B16(00000000, 00011110), \ //#6
                          B16(00000000, 00011111), \ //#7
                          B16(00000000, 00011011), \ //#8
                          B16(00000000, 00011001), \ //#9
                          B16(00000000, 00011000), \ //#10
                          B16(00000000, 00111000), \ //#11
                          B16(00000000, 01111000), \ //#12
                          B16(00000000, 01111001), \ //#13
                          B16(00000000, 01111011), \ //#14
                          B16(00000000, 01111111), \ //#15
                          B16(00000000, 01111110), \ //#16
                          B16(00000000, 01101110), \ //#17
                          B16(00000000, 01100110), \ //#18
                          B16(00000000, 01100111), \ //#19
                          B16(00000000, 01100011), \ //#20
                          B16(00000000, 01100001), \ //#21
                          B16(00000000, 01100000), \ //#22
                          B16(00000000, 11100000), \ //#23
                          B16(00000001, 11100000), \ //#24
                          B16(00000001, 11100001), \ //#25
                          B16(00000001, 11100011), \ //#26
                          B16(00000001, 11100111), \ //#27
                          B16(00000001, 11100110), \ //#28
                          B16(00000001, 11101110), \ //#29
                          B16(00000001, 11111110), \ //#30
                          B16(00000001, 11111111), \ //#31
                          B16(00000001, 11111011), \ //#32
                          B16(00000001, 11111001), \ //#33
                          B16(00000001, 11111000), \ //#34
                          B16(00000001, 10111000), \ //#35
                          B16(00000001, 10011000), \ //#36
                          B16(00000001, 10011001), \ //#37
                          B16(00000001, 10011011), \ //#38
                          B16(00000001, 10011111), \ //#39
                          B16(00000001, 10011110), \ //#40
                          B16(00000001, 10001110), \ //#41
                          B16(00000001, 10000110), \ //#42
                          B16(00000001, 10000111), \ //#43
                          B16(00000001, 10000011), \ //#44
                          B16(00000001, 10000001), \ //#45
                          B16(00000001, 10000000)  \ //#46
                       };

//longpath objective func

float
LongPath_Objective(GAGenome& g) 
{
    GA1DBinaryStringGenome & genome = (GA1DBinaryStringGenome &)g;
    int ones_count=0;
    float final_score = 0.0;

    if(genome.length() != LONGPATH_BITS)
    {
        cout << "Genome length must be" << LONGPATH_BITS << "!!" << endl;
        exit(1);
    }

    genome_pos = -1;

    //count 1s
    for(int i=0; i<genome.length(); i++)
        ones_count += genome.gene(i);

    evals++;

    /* minimizng func
    if(OnPath(genome))
        final_score = longpath_max_val - genome_pos;
    else
        final_score = longpath_max_val + ones_count;
    */
    
    //maximiztion version
    if(OnPath(genome))
        final_score = longpath_max_val + genome_pos;
    else
        final_score = longpath_max_val - ones_count;
    
    return(final_score);  
}

static int OnPath(GA1DBinaryStringGenome& genome)
{
    //convert string to number
    int genome_numvalue = 0;
    genome_pos = -1;

    int place = 8;
    for(int i=0; i<genome.length(); i++)
    {
        genome_numvalue += genome.gene(i) * (int) pow(2,place);
        place--;
    }

    for(int i = 0; i < 47; i++)
    {
        if(longpath_pos[i] == genome_numvalue)
        {
            genome_pos = i;  
            break; //could short circuit!
        }
    }

    if(genome_pos != -1)
        return(true);
    else
        return(false);
}


//**************************************************


#endif /* _LONGPATH_HEADER_H_ */
