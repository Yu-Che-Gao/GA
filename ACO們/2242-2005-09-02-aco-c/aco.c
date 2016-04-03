
#include"dec.h"
#include"func.h"

main(int argc , char *argv[])
{
inpfname = argv[1];
int runno, itno, antno;
double time_used;


read_data(inpfname);
print_data();

init_out(inpfname);


time_used = elapsed_time( VIRTUAL );
printf("Initialization took %.10f seconds\n",time_used);

for(runno=0;runno<runs;runno++)
{
seed = (long int ) time(NULL);
 start_timers();

 initialize_ants_variables(runno);


 initialize_trail();

 for(itno=0;  bestant.ofn != 436.00 &&  itno<ncmax ;itno++)
 {
	 
 iteration_init(itno);
  
 find_values();
  
 analysis( itno);
 

update_stats(itno,runno);

#if (usels == 1)

if(lswithitbest ==1)
{
ls(&itbestval,itno);
lsstats(itno,runno);
}

if(lswithbest ==1)
{
ls(&bestval,itno);
lsstats(itno,runno);
}

#endif

 trail();

if(!(itno%statsafterit))
{ print_itstats(itno ,runno);}

}//end of nc max iterations

update_stats(itno,runno);
report_run(runno);

}//end of runs

final_stats();

freememory();


}//end of main

 
