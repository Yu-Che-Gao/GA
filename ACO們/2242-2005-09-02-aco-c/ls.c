#include"dec.h"
#include"func.h"
#include"ls.h"

extern struct variable *lsvar;
extern struct ant_struct *lsant, *lsbestant, *lsitbestant, *lsitworstant;
extern struct values *lsitbestval, *lsitworstval,*lsbestval;

void ls(struct values *gvalptr ,int globalitno)
{
int  itno, antno;

lsget_data(gvalptr);

lsprint_data();

 lsinitialize_trail();

 for(itno=0;itno<lsncmax ;itno++)
 {
	 
 lsiteration_init(itno);
  
 lsfind_values( );
  
 lsanalysis( itno);
 
 lsupdate_stats(itno);

 lstrail();

 //lsprint_itstats(itno );

 }//end of nc max iterations

lsfreememory();

}//end of main


void lsfreememory()
{
int vno,antno;

for(vno=0;vno<nov;vno++)
{
	free((lsvar+vno)->trail);
	free((lsvar+vno)->prob);
}
/*
 * struct variable *lsvar;
struct ant_struct *lsant, *lsbestant, *lsitbestant, *lsitworstant;
struct values ,*lsitbestval, *lsitworstval;*/

free(lsvar);
free(lsant);
free(lsitbestant);
free(lsitworstant);
free(lsitbestval);
free(lsitworstval);
}//end of free ls mem

void lsprint_itstats(int iterationno)
{

printf("in local search  iteration no %d \n lsbestofn = %15.3f@ iter %4d lsitbestofn = %15.3f \n ", iterationno ,lsbestant->ofn,lsbestant->iter, lsitbestant->ofn);

}


void lsupdate_stats(int iterationno)
{
int vno;

if(lsitbestant->ofn<lsbestant->ofn)
{
copy_ant(lsitbestant,lsbestant);
copy_val(lsitbestval,lsbestval);
}

}// end of update stats 



void lsanalysis(int iterationno)
{
int antno,vno,valueno;
int i;
struct values *tempval;
tempval = (struct values *) malloc(sizeof(struct values));
tempval->iter = iterationno;

for(antno=0;antno<LSANTS;antno++)
{
(lsant+antno)->iter = iterationno;

lsdecode_varvalues ( tempval->var, antno);

#if (constrained ==1)
(lsant+antno)->ofn = objective_function(tempval->var, &(lsant+antno)->penalty );
#endif


#if (constrained ==0)
(lsant+antno)->ofn = objective_function(tempval->var );
#endif


tempval->ofn = (lsant+antno)->ofn;


if((lsant+antno)->ofn < (lsant+antno)->bestofn )
{
 (lsant+antno)->bestofn = (lsant+antno )->ofn;
 (lsant+antno)->iter = iterationno;
}


if((lsant+antno)->ofn < lsitbestant->ofn)
{
copy_val(tempval , lsitbestval);
copy_ant((lsant+antno),lsitbestant);
}

if((lsant+antno)->ofn >lsitworstant->ofn)
{
copy_val(tempval, lsitworstval);
copy_ant((lsant+antno), lsitworstant);
}
}//endof antno loop

free(tempval);
}//end of analysis


void lsdecode_varvalues(float * varvalues,int antno)
{
int vno , valueno;
for(vno=0;vno<nov;vno++)
{
 if((lsant+antno)->variable[vno] == ((lsvar+vno)->nvalues - 1 ))
  {
  varvalues[vno] = (lsvar+vno)->ulimit;
  }
  else
  {
    varvalues[vno] = (  (lsvar+vno)->llimit + (((lsant+antno)->variable[vno]) * (lsvar+vno)->increment ));
  }
}//end of vno loop

}//end of decode values

void lsiteration_init(int iterationno)
{
int antno;
for(antno=0;antno<LSANTS;antno++)
{
(lsant+antno)->iter = iterationno;
}

init_ant(lsitbestant, iterationno, 1e10);
init_ant(lsitworstant, iterationno,0.0);
}//end of iteration init 
/*
void lsinitialize_ants_vars( )
{
int antno;
lsant = (struct ant_struct*)calloc(LSANTS , sizeof(struct ant_struct));
lsitbestant =(struct ant_struct *) malloc(sizeof(struct ant_struct));
lsbestant =(struct ant_struct*) malloc(sizeof(struct ant_struct));
lsitworstant = (struct ant_struct *) malloc(sizeof(struct ant_struct));

lsitworstval =(struct values *) malloc(sizeof(struct values));
lsitbestval = (struct values *)malloc(sizeof(struct values));
lsbestval = (struct values *)malloc(sizeof(struct values));

for(antno=0; antno < LSANTS; antno++)
{
init_ant((lsant+antno),0,0.0);
}

init_ant(lsbestant,0,1e10);

}//ants initialized
*/
void lsget_data(struct values *gvalptr)
{
int vno;

lsvar = (struct variable *) calloc(nov,sizeof(struct variable));

for(vno=0;vno<nov;vno++)
{
	
strcpy((lsvar+vno)->name, (var[vno]).name);
	
(lsvar+vno)->ulimit = gvalptr->var[vno] + (lsnvalues * (var[vno]).increment); 	
if((lsvar+vno)->ulimit > (var[vno]).ulimit )
{
(lsvar+vno)->ulimit = (var[vno]).ulimit;
}

(lsvar+vno)->llimit = gvalptr->var[vno] - (lsnvalues * (var[vno]).increment); 	
if((lsvar+vno)->llimit<(var[vno]).llimit)
{
(lsvar+vno)->llimit = (var[vno]).llimit;
}

(lsvar+vno)->increment = (var[vno]).increment;

lsfind_nvalue(vno);


(lsvar+vno)->trail = (double *) calloc((lsvar+vno)->nvalues ,sizeof(double));
(lsvar+vno)->prob = (double *) calloc((lsvar+vno)->nvalues ,sizeof(double));

}//end of nov loop


int antno;
lsant = (struct ant_struct*)calloc(LSANTS , sizeof(struct ant_struct));
lsitbestant =(struct ant_struct *) malloc(sizeof(struct ant_struct));
lsbestant =(struct ant_struct*) malloc(sizeof(struct ant_struct));
lsitworstant = (struct ant_struct *) malloc(sizeof(struct ant_struct));

lsitworstval =(struct values *) malloc(sizeof(struct values));
lsitbestval = (struct values *)malloc(sizeof(struct values));
lsbestval = (struct values *)malloc(sizeof(struct values));

for(antno=0; antno < LSANTS; antno++)
{
init_ant((lsant+antno),0,0.0);
}

init_ant(lsbestant,0,1e10);


}//end of get data



void lsprint_data()
{
int vno;

char finp[strlen(inpfname)];
strcpy(finp, inpfname);

FILE *input= fopen( strcat(finp, ".lsvar.out" ), "w");
for(vno=0;vno<nov;vno++)
{
fprintf(input," the name for variable no %s \n",(lsvar+vno)->name);
fprintf(input,"the lower limit for %s is %f \n",(lsvar+vno)->name,(lsvar+vno)->llimit);
fprintf(input,"the upper limit for %s is %f  \n",(lsvar+vno)->name,(lsvar+vno)->ulimit);
fprintf(input," the increment for %s is %f \n",(lsvar+vno)->name,(lsvar+vno)->increment);
fprintf(input," the no of values for %s is %d\n",(lsvar+vno)->name,(lsvar+vno)->nvalues);
}
fclose(input);
}//end of lsprint data



void lsfind_nvalue(int vno)
{
float nv,quotient;
int divisor, went= 1 ;


nv =  ( (((lsvar+vno)->ulimit) - ((lsvar+vno)->llimit ))) / ((lsvar+vno)->increment );

divisor = (int) nv;

(lsvar+vno)->nvalues =1 + divisor; 

quotient = nv /(( float)divisor);

if((quotient !=  1.00000000000000));
{
	if ( (fabs((quotient - 1.000000000000000) ) > 0.01))
 	{
 	 ((lsvar+vno)->nvalues) ++;
 	} 
}


}//end of find lsnvalues
	
	
	
 
void lsfind_values()
{
int vno,antno;
for(vno=0;vno<nov;vno++)
{
for(antno=0;antno<LSANTS;antno++)
{ 
(lsant+antno)->variable[vno] =  lsfind_var_valueno(vno);
}
}
}//end of lsfind values


int lsfind_var_valueno(int vno)
{
double q , sum ,ranno;

int valueno, lsnmax;
lsfind_prob(vno);
q=ran01(&seed);
ranno=ran01(&seed);
valueno = lsfind_no_of_max( vno, &lsnmax );
/*printf(" for %d vno, q= %f ,ranno= %f , lsnmax =%d valueno= %d \n",vno,q,ranno,lsnmax,valueno);*/

if((q<lsq0) && (lsnmax==1))
 {
 #if (localupdate==1)
 (lsvar+vno)->trail[valueno] *= (1-lsrho);
 #endif
 return valueno;
 }
else
{
 sum=0.0;
 for(valueno=0;valueno< (lsvar + vno)-> nvalues && sum < ranno ; valueno++)
 {
 sum += (lsvar+vno)->prob[valueno] ;
 }//endfo valueno loop

  #if (localupdate==1)
	if(valueno> 0 && valueno -1 < (lsvar + vno)-> nvalues )
	{
  	(lsvar+vno)->trail[valueno-1] *= (1-rho);
	}
  #endif
 return (valueno-1);
}//end of else

}//end of find valueno



int lsfind_no_of_max(int vno,int *noofmax)
/* returns the value no for which the max value is obtained
 * and changes the value in no of max (lsnmax) */
{
int valueno, chooseno , wentinside = 0, inloop = 0;
*noofmax= (lsvar+vno)->nvalues;
double max=0;
for(valueno=0;valueno< (lsvar+vno)->nvalues; valueno++)
{
if(max < (lsvar+vno)->prob[valueno])
 {
  max =  (lsvar+vno)->prob[valueno];
  chooseno = valueno;
 }
}//endfo valueno loop

for(valueno=0;valueno < (lsvar+vno)->nvalues; valueno++)
{
		
if(max > (lsvar+vno)->prob[valueno])
 {
  wentinside++;
  (*noofmax)= (*noofmax) -1;
 }
}//endfo valueno loop
/*printf("no of wentinside= %d , noofmax= %d\n" , wentinside, *noofmax);*/

return chooseno;
}//end of find no of max
  

void lsfind_prob(int vno)
{
int valueno;
float denominator=0.0;

for(valueno=0;valueno<((lsvar+vno)->nvalues) ; valueno++)
{
denominator += pow((double) (lsvar+vno)->trail[valueno] , (double) alpha);
}

for(valueno=0;valueno<((lsvar+vno)->nvalues) ; valueno++)
{
((lsvar+vno)->prob[valueno]) = pow((double) ((lsvar+vno)->trail[valueno]) , (double) alpha)/denominator;
}

}//end of find prob



void lsevaporate_trail()
{
int vno , valueno;

for(vno=0;vno<nov;vno++)
{
 for(valueno=0;valueno<(lsvar+vno)->nvalues;valueno++)
  {
  ((lsvar+vno)->trail[valueno]) *= (1-lsrho);
  }
}//end of vno loop
}//end of evaporate trail

void lstrail()
{
int antno;

lsevaporate_trail();

for(antno=0;antno<LSANTS;antno++)
{
lsupdate_trail(lsant+antno);
}

#if (elitist ==1)
lsupdate_trail(lsbestant);
#endif

lsupdate_trail(lsitbestant);
}//end of trail


void lsupdate_trail(struct ant_struct * antptr)
{
int vno ,valueno;

for(vno=0;vno<nov;vno++)
{
  valueno= antptr->variable[vno];
 (lsvar+vno)->trail[antptr->variable[vno]] +=  bestant.ofn /*antptr->bestofn*/  / antptr->ofn;
}
}//end of update trail


void lsinitialize_trail()
{
int vno , valueno;

for(vno=0;vno<nov;vno++)
{
for(valueno=0;valueno<((lsvar+vno)->nvalues);valueno++)
 {
 ((lsvar+vno)->trail[valueno]) = lstau0;
  }
}
}//end of trail initialization

