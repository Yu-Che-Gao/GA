#include"dec.h"
#include"func.h"
#include"rand.h"

void find_values()
{
int vno,antno;
for(vno=0;vno<nov;vno++)
{
for(antno=0;antno<ANTS;antno++)
{ 
(ant[antno]).variable[vno] =  find_var_valueno(vno);
/*printf("value chosen for %d ant at %d variable is %d \n",antno,vno, (ant[antno]).variable[vno]);*/
}
}
}//end of find values


int find_var_valueno(int vno)
{
double q , sum ,ranno;

int valueno, nmax;
find_prob(vno);
q=ran01(&seed);
ranno=ran01(&seed);
valueno = find_no_of_max( vno, &nmax );
/*printf(" for %d vno, q= %f ,ranno= %f , nmax =%d valueno= %d \n",vno,q,ranno,nmax,valueno);*/

if((q<q0) && (nmax==1))
 {
 #if (localupdate==1)
 (var[vno]).trail[valueno] *= (1-gamma);
 #endif
 return valueno;
 }
else
{
 sum=0.0;
 for(valueno=0; valueno -1 < var[vno].nvalues && sum < ranno ; valueno++)
 {
 sum += (var[vno]).prob[valueno] ;
 }//endfo valueno loop

  #if (localupdate==1)
  if(valueno > 0 && valueno -1 < var[vno].nvalues)
	{
  	(var[vno]).trail[valueno-1] *= (1-rho);
	}
  #endif
 return (valueno-1);
}//end of else

}//end of find valueno



int find_no_of_max(int vno,int *noofmax)
/* returns the value no for which the max value is obtained
 * and changes the value in no of max (nmax) */
{
int valueno, chooseno , wentinside = 0, inloop = 0;
*noofmax= var[vno].nvalues;
double max=0;
for(valueno=0;valueno< (var[vno]).nvalues; valueno++)
{
if(max < (var[vno]).prob[valueno])
 {
  max =  (var[vno]).prob[valueno];
  chooseno = valueno;
 }
}//endfo valueno loop

/*printf("max prob for %d vno is %f for valueno %d \n\n ", vno,max, chooseno);*/
for(valueno=0;valueno < (var[vno]).nvalues; valueno++)
{
		
if(max > (var[vno]).prob[valueno])
 {
  wentinside++;
  (*noofmax)= (*noofmax) -1;
 }
}//endfo valueno loop
/*printf("no of wentinside= %d , noofmax= %d\n" , wentinside, *noofmax);*/

return chooseno;
}//end of find no of max
  

void find_prob(int vno)
{
int valueno;
float denominator=0.0;

for(valueno=0;valueno<(var[vno]).nvalues ; valueno++)
{
denominator += pow((double) var[vno].trail[valueno] , (double) alpha);
}

for(valueno=0;valueno<(var[vno]).nvalues ; valueno++)
{
var[vno].prob[valueno] = pow((double) var[vno].trail[valueno] , (double) alpha)/denominator;
}

}//end of find prob
