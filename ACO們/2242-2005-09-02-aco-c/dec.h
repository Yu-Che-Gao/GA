#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"timer.h"
#include"ls.h"
#include"parameters.h"

extern FILE *best, *stats, *plot ;
extern FILE* runfile;
extern long int seed;
extern int ex; // the example number
extern long int nevaluations ,nevalbest[runs];
extern 	float midht;
extern clock_t  start_time;
extern double runtime[runs], besttime, passedtime;
extern int tourneypos, tourneysize, tourneylist[ANTS];
extern char *inpfname;

struct variable 
{
	float ulimit, llimit,increment;
	int nvalues;
	double *trail , *prob;
	char name[20];
} ;
extern struct variable var[nov];

struct ant_struct 
{
	int variable[nov];
	float ofn;
	#if(constrained ==1)
	float penalty;
	#endif
	float bestofn;
	int iter;
};

extern struct ant_struct  ant[ANTS],bestant,itbestant,itworstant,runbestant[runs];


struct values 
{
	float var[nov], ofn;
	int iter;
	#if(constrained ==1)
	float penalty;
	#endif
};
extern struct values  itbestval,itworstval,bestval, runbestval[runs];

#if (lbt==1)
extern struct ant_struct lbest[ANTS];
#endif
