#include"dec.h"
#include"func.h"

void print_ant(struct ant_struct * antptr)
{
int vno;
for(vno=0;vno<nov;vno++)
{
printf("%20s = %3d \n",var[vno].name, antptr->variable[vno]);
}
printf("objective funcction = %10.6f \n", antptr->ofn);
printf("iteration no = %4d \n",antptr->iter);
}//end of print ant

void copy_ant(struct ant_struct *from, struct ant_struct *to)
{
int vno;
to->ofn  = from->ofn;
to->iter = from->iter;
to->bestofn= from->bestofn;
for(vno= 0;vno<nov;vno++)
{
to->variable[vno] = from->variable[vno];
}
#if (constrained ==1)
to->penalty = from->penalty;
#endif

}//end of copy ant

void copy_val(struct values  *from ,struct values *to)
{
int vno ;
to->ofn = from->ofn;
to->iter = from->iter;
for(vno=0;vno<nov;vno++)
{
to->var[vno] = from->var[vno];
}
to->iter = from->iter;
}//end of copy values

void lsstats(int iterationno,int runno)
{

if(bestant.ofn>lsbestant->ofn)
{
if(bestval.ofn<lsbestval->ofn)
printf("error in updating in lsstats values\n");

find_lsbestant_valueno(lsbestval, lsbestant);
printf("ls worked in it no=%6d where bestantofn= %10f and lsbestantofn = %10f \n",iterationno, bestant.ofn , lsbestant->ofn);
copy_ant(lsbestant,&bestant);
copy_val(lsbestval,&bestval);
bestant.iter= iterationno;
bestval.iter = iterationno;
nevalbest[runno] = nevaluations;
}

}//end of lsstats

void find_lsbestant_valueno(struct values *lsvalptr, struct ant_struct *lsantptr)
{
int vno,valueno;
for(vno=0;vno<nov;vno++)
{
if(lsvalptr->var[vno] == var[vno].ulimit )
{
lsantptr->variable[vno] = (var[vno]).nvalues -1;
//printf("lsvalptr->var[%d] = %f , var[%d].ulimit = %f  lsantptr->variable[%d] = %d (var[%d]).nvalues-1 = %d\n",vno,lsvalptr->var[vno],vno, var[vno].ulimit, vno,lsantptr->variable[vno] ,vno, (var[vno]).nvalues -1); 
//printf("lsantptr->variable[%d] = %d \n\n", vno, lsantptr->variable[vno]);
}
else 
{
valueno = (int ) ((lsvalptr->var[vno] - var[vno].llimit)/ var[vno].increment);

//printf("valueno = %d lsvalptr->var[%d] = %f , var[%d].llimit = %f  \nlsantptr->variable[%d] = %d (var[%d]).nvalues-1 = %d\n increment = %f",valueno,vno,lsvalptr->var[vno],vno, var[vno].llimit, vno,lsantptr->variable[vno] ,vno, (var[vno]).nvalues -1, var[vno].increment); 

lsantptr->variable[vno] = valueno;
//printf("lsantptr->variable[%d] = %d \n\n", vno, lsantptr->variable[vno]);
}
}
}//end of find lsbestantvalueno

void update_stats(int iterationno,int runno)
{
int vno;
float   mean_value;
double stddev;


if(itbestant.ofn<bestant.ofn)
{
if(itbestval.ofn>bestval.ofn)
printf("in iteration no %d error in updating values\n",iterationno);
copy_ant(&itbestant,&bestant);
copy_val(&itbestval,&bestval);
besttime = elapsed_time( VIRTUAL );
nevalbest[runno] = nevaluations;
/*
ls(&bestval,iterationno);
lsstats(iterationno , runno);
*/
}

passedtime = elapsed_time(VIRTUAL);  
//fprintf(plot,"%4d %15.6f %15.6f  %3.15f \n ",iterationno,itbestant.ofn,bestant.ofn,passedtime);
fprintf(plot,"%15d",iterationno );
fprintf(plot,"%15f",itbestant.ofn );
fprintf(plot,"%15f",bestant.ofn );
fprintf(plot,"%15f",passedtime );

for(vno=0;vno<nov;vno++)
{
fprintf(plot,"%15.3f",  bestval.var[vno]);
}
for(vno=0;vno<nov;vno++)
{
fprintf(plot,"%15.3f",  itbestval.var[vno]);
}

fprintf(plot,"\n");
mean_value  = find_mean();
stddev = std_deviation_ants( (double) mean_value  );
//fprintf( stats , " %6d   %15.3f  %15.3f  %15.3f  %15.3f %3.15f \n", iterationno , itbestant.ofn, itworstant.ofn ,mean_value , stddev , passedtime);
fprintf(stats,"%15d",iterationno );
fprintf(stats,"%15f",itbestant.ofn );
fprintf(stats,"%15f",itworstant.ofn );
fprintf(stats,"%15f",mean_value );
fprintf(stats,"%15f",stddev );
fprintf(stats,"%15f",passedtime );
fprintf(stats,"\n");

}// end of update stats 



double std_deviation_ants(  double d_mean ) 
{
  long int j;
  double   dev;

  dev = 0.;
  for ( j = 0 ; j < ANTS; j++ ) {
    dev += ((double)(ant[j]).ofn - d_mean) * ((double)(ant[j]).ofn - d_mean);
  }
  return sqrt(dev/(double)(ANTS- 1));
}//end fo std deviation ants



float find_mean()
{
float  sum=0.0;
int i;
for(i=0;i<ANTS; i++)
{
 sum += (ant[i]).ofn;
}

return (sum / ANTS);
}//end of mean


void print_itstats(int iterationno, int runno)
{

printf("runno is %d iteration no %d \n bestofn = %15.3f@ iter %4d itbestofn = %15.3f \n ", runno ,iterationno ,bestant.ofn,bestant.iter, itbestval.ofn);

}

	FILE *fpath;
void report_run (int runno)
{

	char finpfile[strlen(inpfname)];
	if(runno == 0)
	{
	strcpy(finpfile, inpfname);
	fpath= fopen( strcat(finpfile, ".path.out" ), "w");
	}
int vno;
printf("Reporting run no = %d \n",runno+1);
fprintf(plot ,"#Runno is %d \n",runno+1);
fprintf(stats,"#Runno is %d \n",runno+1);

passedtime = elapsed_time(VIRTUAL);  

copy_ant(&bestant,&runbestant[runno]);
copy_val(&bestval,&runbestval[runno]);
runtime[runno] = besttime;
fprintf(best,"\\begin{table}[!hp] \\label{result%d}\n  \\begin{center}\n \\caption{Results obtained for run no %d with %s }\n \\begin{tabular}{|r|l|} \\hline \n", 1+ runno, 1+ runno, inpfname);

fprintf(best,"%15s&","Run-No");

fprintf(best,"%15d\n",1+ runno);
fprintf(best,"\\\\\\hline\n");
fprintf(best,"%15s&","Time-taken");
fprintf(best,"%15f\\\\\n", passedtime);

#if (constrained ==1)
fprintf(best,"%15s& ","Best-ofn");
fprintf(best,"%15f\\\\\n",bestval.ofn);
fprintf(best,"%15s& ","Penalty");
fprintf(best,"%15f\\\\\n",bestant.penalty);
fprintf(best,"%15s& ","Iteration");
fprintf(best,"%15d\\\\\n",1+ bestval.iter);
fprintf(best,"%15s& ","Evaluations");
fprintf(best,"%15d\\\\\n",nevalbest[runno]);
fprintf(best,"%15s& ","After-time");
fprintf(best,"%15f\\\\\n",besttime);
//fprintf(best,"best ofn = %15.6f with penaalty %15.15f\\\\ found at iter %d after %ld evaluations ,after time   =  %3.15f\\\\\n\n",bestval.ofn,bestant.penalty,1+bestval.iter,nevalbest[runno],besttime);
#endif
#if (constrained ==0)
fprintf(best,"%15s& ","Best-Ofn");
fprintf(best,"%15f\\\\\n",bestval.ofn);
fprintf(best,"%15s& ","Iteration");
fprintf(best,"%15f\\\\\n",1+ bestval.iter);
fprintf(best,"%15s& ","Evaluations");
fprintf(best,"%15f\\\\\n",nevalbest[runno]);
fprintf(best,"%15s& ","After-time");
fprintf(best,"%15f\\\\\n",besttime);
//fprintf(best,"best ofn = %15.6f  found at iter %d after %ld evaluations ,after time   =  %3.15f\n\n",bestval.ofn,1+bestval.iter,nevalbest[runno],besttime);
#endif

fprintf(best,"$\\mathcal{X}_i$& value \n");
fprintf(best,"\\\\\\hline\n");
//fprintf(fpath, "%15d %15f\n",0,midht  );

fprintf(fpath,"set arrow from %d,%f to %d,%f size 0.3,30 \n", 0 , midht ,  1 ,  midht + ( var[0].nvalues - 1.0 ) / 2.0  -  bestant.variable[0]);
for(vno=0;vno<nov;vno++)
	{
	fprintf(best,"%15s& %15.3f\\\\\n", var[vno].name, bestval.var[vno]);
	
	if(vno < nov -1 )
	{
	fprintf(fpath,"set arrow from %d,%f to %d,%f size 0.3,30 \n", vno+1 , midht + ( var[vno].nvalues - 1.0 ) / 2.0  -  bestant.variable[vno],  vno+2 ,  midht + ( var[1+ vno].nvalues - 1.0 ) / 2.0  -  bestant.variable[1+ vno]);
	}
	
	//fprintf(fpath,"%15d %15f\n", vno+1 , midht + ( var[vno].nvalues - 1.0 ) / 2.0  -  bestant.variable[vno]);
	}

	vno-- ;
	fprintf(fpath,"set arrow from %d,%f to %d,%f  size 0.3,30 \n", vno+1 , midht + ( var[vno].nvalues - 1.0 ) / 2.0  -  bestant.variable[vno],  vno+2 ,  midht ); 
	//fprintf(fpath, "%15d %15f\n",vno+1,midht  );
	fprintf(fpath, "\n");



fprintf(best,"\\hline \n \\end{tabular}\n \\end{center}\n \\end{table}\n");

	if(runno == runs - 1 )
	{
	fclose(fpath);
	}

	fprintf(plot,"\n");
}//run reported


void final_stats()

{
char finp[strlen(inpfname)];
strcpy(finp, inpfname);
FILE *table= fopen( strcat(finp, ".table.tex" ), "w");

int varno;

fprintf(table,"\n");
fprintf(table,"\n");
fprintf(table,"\\documentclass{article} \n \\usepackage[dvips]{graphics}\\usepackage{graphicx} \n \\begin{document} \nInput file name is %15s\n\n",inpfname);

fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.constructiongraph.eps}} \n \\caption{Convergence of objective function for %s.} \\label{conv-ofn} \\end{figure}\n",inpfname , inpfname);
	fprintf(table," \\input{%s.best.tex}", inpfname);
	fprintf(table," \\input{%s.var.tex}", inpfname);

fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.obj.eps}} \n \\caption{Convergence of objective function for %s.} \\label{conv-ofn} \\end{figure}\n",inpfname , inpfname);
fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.var.eps}} \n \\caption{Convergence of variable values for %s.} \\label{conv-var} \\end{figure}\n",inpfname , inpfname);
fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.dev.eps}} \n \\caption{Convergence history of run time statistics for %s.} \\label{conv-run} \\end{figure}\n",inpfname , inpfname);
fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.time.eps}} \n \\caption{Time required to find best solution for %d runs with  %s.} \\label{time} \\end{figure}\n",inpfname , runs,inpfname);

/*
\begin{figure}[!htp]
\centering
\scalebox{0.7}{\includegraphics{bwants.eps}}
\caption{Biological ants find shortest path by using pheromone based communication.}
\label{antspath}
\end{figure}
*/

fprintf(table,"\\begin{table}[!hp] \\label{allresults}\n  \\begin{center}\n \\caption{Results obtained for  %d runs with %s }\n \\begin{tabular}{|", runs, inpfname);
for(varno=0;varno<nov+5;varno++)
{fprintf(table,"|c|");
}
#if(constrained==1)
	fprintf(table,"|c|");
#endif
fprintf(table,"|} \\hline \\hline \n");


fprintf(table,"%15s&","\\#Run");
fprintf(table,"%15s&","Obj-function");
#if(constrained==1)
fprintf(table,"%15s&","penalty");
#endif
fprintf(table,"%15s&"," Iterations");
fprintf(table,"%15s&"," After time");
fprintf(table,"%15s&"," Reqd. Eval's");
for(varno=0;varno<nov;varno++)
{
if(varno!=(nov-1))
fprintf(table,"%15s&",var[varno].name);
else
fprintf(table,"%15s",var[varno].name);
}
fprintf(table,"\\\\ \\hline \\hline\n");


fprintf(runfile,"%15s","#Run");
fprintf(runfile,"%15s","Obj-function");
#if(constrained==1)
fprintf(runfile,"%15s","penalty");
#endif
fprintf(runfile,"%15s"," Iterations");
fprintf(runfile,"%15s"," After time");
fprintf(runfile,"%15s"," Reqd. Eval's");
for(varno=0;varno<nov;varno++)
{
if(varno!=(nov-1))
fprintf(runfile,"%15s",var[varno].name);
else
fprintf(runfile,"%15s\n",var[varno].name);
}


printf("Writing  final stats \n");
int i,vno , runno ;
double meanofn=0.0, meaniter = 0.0,meaneval=0.0, meantime=0.0, bestrunofn, worstrunofn;
double worstruntime,bestruntime;
double ofnarray[runs];
double stddev,stditer,stdeval, stdtime;
long int iterarray[runs], bestruniter, worstruniter,worstruneval=0 , bestruneval= ANTS*ncmax;
worstruntime = runtime[0];
bestruntime = runtime[0];
bestrunofn = 1e25;
worstrunofn = 0.0;
worstruniter = 1;
bestruniter = ncmax;


for(i=0;i<runs;i++)
{
	 meanofn +=  (runbestval[i]).ofn;
  	meaniter += (runbestval[i]).iter;
  	meaneval += nevalbest[i];
  	meantime += runtime[i];
  	ofnarray[i] =  (double) (runbestval[i]).ofn;
  	iterarray[i] = (long int )(1+ (runbestval[i]).iter);
	
	fprintf(table,"%15d&",1+i);
	fprintf(runfile,"%15d",1+i);
	
	#if (constrained ==1)
	fprintf(table,"%15f&", (runbestval[i]).ofn);
	fprintf(table,"%15f&", (runbestval[i]).penalty);
	fprintf(table,"%15d&",(runbestval[i]).iter );
	fprintf(table,"%15f&",runtime[i] );
	fprintf(table,"%15ld&",nevalbest[i] );
	#endif
	#if (constrained ==0)
		//fprintf(table,"%10.6f  %4d  %10f %10ld  &",(runbestval[i]).ofn,(runbestval[i]).iter,runtime[i],nevalbest[i]);
	fprintf(table,"%15f&", (runbestval[i]).ofn);
	fprintf(table,"%15d&",(runbestval[i]).iter );
	fprintf(table,"%15f&",runtime[i] );
	fprintf(table,"%15ld&",nevalbest[i] );
	#endif



	#if (constrained ==1)
	fprintf(runfile,"%15f", (runbestval[i]).ofn);
	fprintf(runfile,"%15f", (runbestval[i]).penalty);
	fprintf(runfile,"%15d",(runbestval[i]).iter );
	fprintf(runfile,"%15f",runtime[i] );
	fprintf(runfile,"%15ld",nevalbest[i] );
	#endif
	#if (constrained ==0)
	fprintf(runfile,"%15f", (runbestval[i]).ofn);
	fprintf(runfile,"%15d",(runbestval[i]).iter );
	fprintf(runfile,"%15f",runtime[i] );
	fprintf(runfile,"%15ld",nevalbest[i] );
	#endif




	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"%15.3f&", (runbestval[i]).var[vno]);
		else
		fprintf(table,"%15.3f", (runbestval[i]).var[vno]);

		if(vno!= (nov-1))
		fprintf(runfile,"%15.3f", (runbestval[i]).var[vno]);
		else
		fprintf(runfile,"%15.3f", (runbestval[i]).var[vno]);
	}
	
  	if(bestruntime > runtime[i])
  	{
  		bestruntime = runtime[i];
  	}
  	if (worstruntime < runtime[i])
        {
	 	worstruntime = runtime[i];
	}


	if(bestrunofn > ofnarray[i])
  	{
  	bestrunofn = ofnarray[i];
  	}
  	if (worstrunofn < ofnarray[i])
         {
	 	worstrunofn = ofnarray[i];
	 }


	if(bestruniter > iterarray[i])
  	{
  		bestruniter = iterarray[i];
  	}
  	if (worstruniter < iterarray[i])
         {
	 	worstruniter = iterarray[i];
	 }

	if (bestruneval > nevalbest[i])
  	{
	       	bestruneval = nevalbest[i];
	}
	if(worstruneval < nevalbest[i])
  	{
	       	worstruneval = nevalbest[i];
	}
	
	

	fprintf(table,"\\\\ \\hline \n");
	fprintf(runfile,"\n");
	
 	}
	
	meanofn /= runs;
 	meaniter /= runs;
 	meaneval /= runs;
 	meantime /= runs;
	stddev = std_deviationr( ofnarray , runs , (double)meanofn);
	stdtime = std_deviationr( runtime , runs , (double)meantime);
	stditer = std_deviation( iterarray, runs, (double)meaniter);
	stdeval = std_deviation(nevalbest, runs, (double)meaneval);



	fprintf(table,"%15s&","mean");
	fprintf(table,"%15f&",meanofn );
	#if(constrained==1)
		fprintf(table,"& ");
	#endif
	fprintf(table,"%15f&",meaniter );
	fprintf(table,"%15f&",meantime );
	fprintf(table,"%15f&",meaneval );
	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"& ");
		else
		fprintf(table,"\\\\ \\hline \n");

	}


	fprintf(table,"%15s&","best");
	fprintf(table,"%15f&",bestrunofn );
	#if(constrained==1)
		fprintf(table,"& ");
	#endif
	fprintf(table,"%15d&",bestruniter );
	fprintf(table,"%15f&",bestruntime );
	fprintf(table,"%15d&",bestruneval );
	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"& ");
		else
		fprintf(table,"\\\\ \\hline \n");

	}



	fprintf(table,"%15s&"," worst");
	fprintf(table,"%15f&", worstrunofn );
	#if(constrained==1)
		fprintf(table,"& ");
	#endif
	fprintf(table,"%15d&", worstruniter );
	fprintf(table,"%15f&", worstruntime );
	fprintf(table,"%15d&", worstruneval );
	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"& ");
		else
		fprintf(table,"\\\\ \\hline \n");

	}





	fprintf(table,"%15s&"," std-dev");
	fprintf(table,"%15f&", stddev);
	#if(constrained==1)
		fprintf(table,"& ");
	#endif
	fprintf(table,"%15f&", stditer);
	fprintf(table,"%15f&", stdtime);
	fprintf(table,"%15f&", stdeval);
	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"& ");
		else
		fprintf(table,"\\\\ \\hline \n");

	}





	fprintf(table,"%15s&"," ave/best");
	fprintf(table,"%15f&",((meanofn - bestrunofn)/(1e-10 + bestrunofn)) );
	#if(constrained==1)
		fprintf(table,"& ");
	#endif
	fprintf(table,"%15f&", ((meaniter- bestruniter)/(bestruniter+1e-10)));   
	fprintf(table,"%15f&", ((meantime- bestruntime)/(bestruntime+1e-10))); 
	fprintf(table,"%15f&", ((meaneval- bestruneval)/(bestruneval+1e-10))); 
	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"& ");
		else
		fprintf(table,"\\\\ \\hline \n");

	}

	fprintf(table,"%15s&","worst/best");
	fprintf(table,"%15f&", ((worstrunofn-bestrunofn)/(1e-10+ bestrunofn))); 
	#if(constrained==1)
		fprintf(table,"& ");
	#endif
	fprintf(table,"%15f&", ((worstruniter-bestruniter)/(bestruniter + 1e-10)));   
	fprintf(table,"%15f&",  ((worstruntime-bestruntime)/(1e-10+ bestruntime))); 
	fprintf(table,"%15f&",  ((worstruneval-bestruneval)/(1e-10+ bestruneval))); 
	for(vno=0;vno<nov;vno++)
	{
		if(vno!= (nov-1))
		fprintf(table,"& ");
		else
		fprintf(table,"\\\\ \\hline \n");

	}

	
	fprintf(table,"\\hline \n \\end{tabular}\n \\end{center}\n \\end{table}\n");
	fprintf(table,"\n");
	fprintf(table,"\n");


	

	
	fprintf(runfile,"%15s","#");
	fprintf(runfile,"%15s","#obj-function");
	fprintf(runfile,"%15s","#Iterations");
	fprintf(runfile,"%15s","#Time");
	fprintf(runfile,"%15s","#Evaluations");
	fprintf(runfile,"\n");
	fprintf(runfile,"%15s","# mean");
	fprintf(runfile,"%15s","best");
	fprintf(runfile,"%15s"," worst");
	fprintf(runfile,"%15s"," std-dev");
	fprintf(runfile,"%15s"," ave/best");
	fprintf(runfile,"%15s"," worst/best");
	fprintf(runfile,"\n#");
	fprintf(runfile,"%15f",meanofn );
	fprintf(runfile,"%15f",bestrunofn );
	fprintf(runfile,"%15f", worstrunofn );
	fprintf(runfile,"%15f", stddev);
	fprintf(runfile,"%15f",((meanofn - bestrunofn)/(1e-10 + bestrunofn)) );
	fprintf(runfile,"%15f", ((worstrunofn-bestrunofn)/(1e-10+ bestrunofn))); 
	fprintf(runfile,"\n#");
	
	
	//fprintf(best, "\% final statistics  \n");
	//fprintf(best, "\%meanofn = %15.15f   best ofn = %15.15f , worstofn = %15.15f standard-deviation ofn = %15.15f\n\n \n", meanofn,bestrunofn, worstrunofn, stddev  );
		//fprintf(best, "\%average ofn  = %3.15f more than bestofn , worstofn  = %3.15f more than bestofn\n\n \n",((meanofn- bestrunofn)/(1e-10+bestrunofn)), ((worstrunofn-bestrunofn)/(1e-10+ bestrunofn)));  

	//fprintf(best, "\%mean iterations = %6.2f bestrun iter = %6d worstruniter = %6d standard dev iterations = %6.2f \n\n\n", meaniter,bestruniter,worstruniter,  stditer);
	
	//fprintf(best, "\%average iter = %3.15f more than best iter , worstiter = %3.15f more than bestiter \n\n\n",((meaniter- bestruniter)/(bestruniter+1e-10)), ((worstruniter-bestruniter)/(bestruniter + 1e-10)));    
	







	
	fprintf(table,"\\begin{table}[!hp] \\label{stats}\n  \\begin{center}\n \\caption{Statistics for %d runs with %s}\n \\begin{tabular}{||c||c||c||c||c||c||c||}\n \\hline \\hline \n",runs , inpfname);
	fprintf(table,"%15s&"," ");
	fprintf(table,"%15s&","mean");
	fprintf(table,"%15s&","best");
	fprintf(table,"%15s&"," worst");
	fprintf(table,"%15s&"," std-dev");
	fprintf(table,"%15s&"," ave/best");
	fprintf(table,"%15s"," worst/best");
	fprintf(table,"\\\\ \\hline \n");
	fprintf(table,"%15s&","obj-function");
	fprintf(table,"%15f&",meanofn );
	fprintf(table,"%15f&",bestrunofn );
	fprintf(table,"%15f&", worstrunofn );
	fprintf(table,"%15f&", stddev);
	fprintf(table,"%15f&",((meanofn - bestrunofn)/(1e-10 + bestrunofn)) );
	fprintf(table,"%15f", ((worstrunofn-bestrunofn)/(1e-10+ bestrunofn))); 
	fprintf(table,"\\\\ \\hline \n");
	

	
	fprintf(runfile,"%15s","#");
	fprintf(runfile,"%15s","#obj-function");
	fprintf(runfile,"%15s","#Iterations");
	fprintf(runfile,"%15s","#Time");
	fprintf(runfile,"%15s","#Evaluations");
	fprintf(runfile,"\n");
	fprintf(runfile,"%15s","# mean");
	fprintf(runfile,"%15s","best");
	fprintf(runfile,"%15s"," worst");
	fprintf(runfile,"%15s"," std-dev");
	fprintf(runfile,"%15s"," ave/best");
	fprintf(runfile,"%15s"," worst/best");
	fprintf(runfile,"\n#");
	fprintf(runfile,"%15f",meanofn );
	fprintf(runfile,"%15f",bestrunofn );
	fprintf(runfile,"%15f", worstrunofn );
	fprintf(runfile,"%15f", stddev);
	fprintf(runfile,"%15f",((meanofn - bestrunofn)/(1e-10 + bestrunofn)) );
	fprintf(runfile,"%15f", ((worstrunofn-bestrunofn)/(1e-10+ bestrunofn))); 
	fprintf(runfile,"\n#");
	
	
	//fprintf(best, " \%final statistics  \n");
	//fprintf(best, "\%meanofn = %15.15f   best ofn = %15.15f , worstofn = %15.15f standard-deviation ofn = %15.15f\n\n \n", meanofn,bestrunofn, worstrunofn, stddev  );
		//fprintf(best, "\%average ofn  = %3.15f more than bestofn , worstofn  = %3.15f more than bestofn\n\n \n",((meanofn- bestrunofn)/(1e-10+bestrunofn)), ((worstrunofn-bestrunofn)/(1e-10+ bestrunofn)));  

	//fprintf(best, "\%mean iterations = %6.2f bestrun iter = %6d worstruniter = %6d standard dev iterations = %6.2f \n\n\n", meaniter,bestruniter,worstruniter,  stditer);
	
	//fprintf(best, "\%average iter = %3.15f more than best iter , worstiter = %3.15f more than bestiter \n\n\n",((meaniter- bestruniter)/(bestruniter+1e-10)), ((worstruniter-bestruniter)/(bestruniter + 1e-10)));    
	
	fprintf(table,"%15s&","Iterations");
	fprintf(table,"%15f&",meaniter );
	fprintf(table,"%15d&",bestruniter );
	fprintf(table,"%15d&", worstruniter );
	fprintf(table,"%15f&", stditer);
	fprintf(table,"%15f&", ((meaniter- bestruniter)/(bestruniter+1e-10)));   
	fprintf(table,"%15f", ((worstruniter-bestruniter)/(bestruniter + 1e-10)));   
	fprintf(table,"\\\\ \\hline \n");
	
	
	
	fprintf(runfile,"%15f",meaniter );
	fprintf(runfile,"%15d",bestruniter );
	fprintf(runfile,"%15d", worstruniter );
	fprintf(runfile,"%15f", stditer);
	fprintf(runfile,"%15f", ((meaniter- bestruniter)/(bestruniter+1e-10)));   
	fprintf(runfile,"%15f", ((worstruniter-bestruniter)/(bestruniter + 1e-10)));   
	fprintf(runfile,"\n#");
	
	
//fprintf(best, "\%mean best time = %6.15f besttime = %6.15f worsttime = %6.15f standard dev besttime  = %6.6f \n\n\n", meantime, bestruntime, worstruntime ,stdtime);
	
	//fprintf(best, "\%average runtime = %3.15f more than best runtime , worstruntime = %3.15f more than bestruntime \n\n\n",((meantime- bestruntime)/(bestruntime+1e-10)), ((worstruntime-bestruntime)/(1e-10+ bestruntime)));  
	
	fprintf(table,"%15s&","Time");
	fprintf(table,"%15f&",meantime );
	fprintf(table,"%15f&",bestruntime );
	fprintf(table,"%15f&", worstruntime );
	fprintf(table,"%15f&", stdtime);
	fprintf(table,"%15f&", ((meantime- bestruntime)/(bestruntime+1e-10))); 
	fprintf(table,"%15f",  ((worstruntime-bestruntime)/(1e-10+ bestruntime))); 
	fprintf(table,"\\\\ \\hline \n");
	

	fprintf(runfile,"%15f",meantime );
	fprintf(runfile,"%15f",bestruntime );
	fprintf(runfile,"%15f", worstruntime );
	fprintf(runfile,"%15f", stdtime);
	fprintf(runfile,"%15f", ((meantime- bestruntime)/(bestruntime+1e-10))); 
	fprintf(runfile,"%15f",  ((worstruntime-bestruntime)/(1e-10+ bestruntime))); 
	fprintf(runfile,"\n#");

	
	
	//fprintf(best, "\%mean evaluations  = %6.15f besteval = %6d worsteval = %6d standard dev besteval  = %6.6f \n\n\n", meaneval, bestruneval, worstruneval ,stdeval);

	
	//fprintf(best, "\%average evaluations = %3.15f more than best evaluations , worst evaluations = %3.15f more than bestevaluations \n\n\n\n",((meaneval- bestruneval)/(bestruneval+ 1e-10)), ((worstruneval-bestruneval)/(bestruneval+1e-10)));  
	
	fprintf(table,"%15s&","Evaluations");
	fprintf(table,"%15f&",meaneval );
	fprintf(table,"%15d&",bestruneval );
	fprintf(table,"%15d&", worstruneval );
	fprintf(table,"%15f&", stdeval);
	fprintf(table,"%15f&", ((meaneval- bestruneval)/(bestruneval+1e-10))); 
	fprintf(table,"%15f ",  ((worstruneval-bestruneval)/(1e-10+ bestruneval))); 
	fprintf(table,"\\\\ \\hline \n\n\n\n");
	fprintf(table,"\\hline \n \\end{tabular}\n \\end{center}\n \\end{table}\n ");
	fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.run.eps}} \n \\caption{Values of Best Objective function found for %d runs} \\label{obj-run} \\end{figure}\n", inpfname, runs);
	fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.eval.eps}} \n \\caption{Evaluations requred to find the best solution for each of  %d runs.} \\label{eval-run}\\end{figure}\n", inpfname, runs );
	fprintf(table,"\\begin{figure}[!htp]\n \\centering \\scalebox{1.0}{\\includegraphics{%s.iter.eps}} \n \\caption{Iterations requred to find the best solution for each of %d runs.} \\label{iter-run} \\end{figure}\n", inpfname , runs );

	fprintf(table," \\end{document}");
	fclose(table);
	
	
	fprintf(runfile,"%15f",meaneval );
	fprintf(runfile,"%15d",bestruneval );
	fprintf(runfile,"%15d", worstruneval );
	fprintf(runfile,"%15f", stdeval);
	fprintf(runfile,"%15f", ((meaneval- bestruneval)/(bestruneval+1e-10))); 
	fprintf(runfile,"%15f",  ((worstruneval-bestruneval)/(1e-10+ bestruneval))); 
	fprintf(runfile,"\n#");
	fclose(runfile);
	
	

}//end of final stats


double std_deviationr( double *values, long int max, double d_mean) 
{
  long int j;
  double   dev;

  if (max <= 1)
    return 0.;
  dev = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    dev += ((double)values[j] - d_mean) * ((double)values[j] - d_mean);
  }
  return sqrt(dev/(double)(max - 1));
}



double std_deviation( long int *values, long int max, double i_mean ) 
/*    
      FUNCTION:       compute the standard deviation of an integer array  
      INPUT:          pointer to array, length of array, mean 
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   dev = 0.;

  if (max <= 1)
    return 0.;
  for ( j = 0 ; j < max; j++ ) {
    dev += ((double)values[j] - i_mean) * ((double)values[j] - i_mean);
  }
  return sqrt(dev/(double)(max - 1));
}
