#include"dec.h"
#include"func.h"

FILE* runfile;
void init_out(char * fin )
{
int varno;
char finp[strlen(fin)];
char finp1[strlen(fin)];
char finp2[strlen(fin)];

char finprun[strlen(fin)];
strcpy(finprun, fin);
runfile= fopen( strcat(finprun, ".run.out" ), "w");
strcpy(finp, fin);
strcpy(finp1, fin);
strcpy(finp2, fin);
printf("%s\n",fin); 
printf("%s\n",finp); 
printf("%s\n",finp1); 
printf("finp2 is %s\n",finp2); 
if(finp !=finp1)
printf("different names \n\n\n\n"); 

best= fopen( strcat(finp, ".best.tex" ), "w");
stats = fopen( strcat(finp1, ".stat.out" ), "w");
plot = fopen( strcat(finp2, ".plot.out" ), "w");
printf("strcat gives %s\n\n", strcat(finp2, ".plot.out" ));
//getchar();
if(plot == NULL)
{
	plot = fopen( strcat(finp, ".plot.out" ), "w");
	printf("error in opening plot \n");
}

if(best == NULL)
printf("error in opening best \n");

if(stats == NULL)
printf("error in opening stats \n");

int ranksinras,noofelitist,noofvars,noants,   nr;
int maxcycles,nantsls, nvaluesls, maxcycls; 
float evap , lsevap ,localevap,alp, initialtrail, que; 
noants = ANTS;	
nr = runs;
noofvars = nov;
ranksinras = ras_ranks;
noofelitist = elitist_ants;
initialtrail = tau0;
maxcycles = ncmax;
evap = rho;
localevap = gamma;
que = q0;
alp = alpha;



fprintf(best,"\\begin{table}[!hp] \\label{algodetails}\n  \\begin{center}\n \\caption{Algorithm details  with %s }\n \\begin{tabular}{|l|r|} \\hline \n", inpfname);
if(lbt)
{
fprintf(best,"%15s&","Algorithm used:");
fprintf(best, "Local Best Tour \\\\\n & Ant Algorithm\\\\\n");
}

if(elitist)
{
fprintf(best,"%15s&","Algorithm used:");
fprintf(best, "Elitist Ant System\\\\\n");
fprintf(best, "Elitist Ants & %d \\\\\n",noofelitist );
}
if(ras)
{
fprintf(best,"%15s&","Algorithm used:");
fprintf(best, "Rank Based \\\\\n &Ant System\\\\\n");
fprintf(best, "\\hline\n");
fprintf(best, "Ranks Used & %d \\\\ \n Elitist Ants & %d \\\\\n", ranksinras,noofelitist );
}
if(taco)
{
fprintf(best,"%15s&","Algorithm used:");
fprintf(best, "Tournament Based Ant \\\\\n &Colony Optimization with \\\\\n Tournament size & %d \\\\\n", tourneysize);
}
if(as)
{
fprintf(best,"%15s&","Algorithm used:");
fprintf(best, "Basic Ant System\\\\\n");
}
if(acs)
{
fprintf(best,"%15s&","Algorithm used:");
fprintf(best, "Ant Colony System \\\\\n q0 & %f\\\\\n",que);
}
fprintf(best, "\\hline\n");

fprintf(best, "No of variables &  %d,\\\\\nNo of ants & %d \\\\\n", noofvars,noants );
fprintf(best, "Max no of cycles &  %d\\\\\n Max no of runs & %d\\\\\nEvaporation (rho) &  %f \\\\\n",maxcycles, nr ,evap );
fprintf(best,"Initial trail & %f \\\\\n",initialtrail );
fprintf(best, "\\hline\n");
if(localupdate)
{
fprintf(best, "Local update & used\\\\\n");
}
fprintf(best, "$\\alpha$ &  %f,\\\\\nLocal update evap( $\\gamma$) &  %f \\\\\n \\hline ",alp ,localevap);
if(usels)
{
fprintf(best, "Local search & \\\\\n used with& ");
if(lswithbest)
fprintf(best, "Global Best Ant\\\\\n");
fprintf(best, "\\hline\n");



if(lswithitbest)
fprintf(best, "Iteration Best Ant\\\\\n");
nvaluesls= lsnvalues;
maxcycls = lsncmax;
lsevap = lsrho;
nantsls= LSANTS;
fprintf(best, "In local search & \\\\ \\hline \n Max no of Cycles &  %d\\\\\nNo of Ants  & %d\\\\\nEvaporation ($ls \\rho$ ) &  %f \\\\\n",maxcycls, nantsls ,lsevap );
}
fprintf(best,"\\hline \n \\end{tabular}\n \\end{center}\n \\end{table}\n");

fprintf(best, "\n\n\n");


//fprintf(stats ,"#iteration no #iteration best #iteration worst  #meanval;ues # std dev for iter #time elapsed\n");
fprintf(stats,"%15s","#Iteration no.");
fprintf(stats,"%15s","Iteration Best");
fprintf(stats,"%15s","Iteration worst");
fprintf(stats,"%15s","Mean ofn");
fprintf(stats,"%15s","Std-dev");
fprintf(stats,"%15s","Time elapsed");
fprintf(stats,"\n");

//fprintf(plot ,"#iteration no #iteration best #global best \n");
fprintf(plot,"%15s","#Iteration no.");
fprintf(plot,"%15s","Iteration Best");
fprintf(plot,"%15s","Global Best");
fprintf(plot,"%15s","Time elapsed");
for(varno=0;varno<nov;varno++)
{fprintf(plot,"%12s.gb",var[varno].name);
}
for(varno=0;varno<nov;varno++)
{fprintf(plot,"%12s.ib",var[varno].name);
}
fprintf(plot,"\n");

/*
fprintf(runfile,"%15s","\\#run ");
fprintf(runfile,"%15s","#obj-function");
fprintf(runfile,"%15s","#Iterations");
fprintf(runfile,"%15s","#Time");
fprintf(runfile,"%15s","#Evaluations");
for(varno=0;varno<nov;varno++)
	{
		fprintf(runfile,"%15s",var[varno].name);
	}
fprintf(runfile,"\n");

*/


}//end of init_out 


void initialize_ants_variables(int runno)
{
int antno;
nevaluations = 0;
nevalbest[runno] = 0;
/*fprintf(plot,"run no is %4d\n",runno);*/
for(antno=0; antno < ANTS; antno++)
{
init_ant(&ant[antno],0,0.0);
#if (lbt ==1)
init_ant(&lbest[antno],0,1e25);
#endif
}
init_ant(&bestant,0,1e25);

bestval.ofn = 1e25;
}//ants initialized

void iteration_init(int iterationno)
{
int antno;
for(antno=0;antno<ANTS;antno++)
{
ant[antno].iter = iterationno;
}

init_ant(&itbestant, iterationno, 1e25);
itbestval.ofn = 1e25;
init_ant(&itworstant, iterationno,0.0);
}//end of iteration init 

 void init_ant(struct ant_struct *antptr, int iterationno,float iofn)
{
antptr->ofn = iofn;
antptr->iter= iterationno;
antptr->bestofn= 1e25;
#if (constrained ==1)
antptr->penalty = 0.0;
#endif

}//endof init ant

void freememory()
{
int vno;

for(vno=0;vno<nov;vno++)
{
free(var[vno].trail);
free(var[vno].prob);
}

}

