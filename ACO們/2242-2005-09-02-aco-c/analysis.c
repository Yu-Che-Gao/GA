
#include"dec.h"
#include"func.h"

void analysis(int iterationno)
{
int antno,vno,valueno;
struct values tempval;
tempval.iter = iterationno;
for(antno=0;antno<ANTS;antno++)
{
ant[antno].iter = iterationno;
decode_varvalues ( tempval.var, antno);
/*
for(vno=0;vno<nov;vno++)
{
printf("variable no %d = %f for antno=%d \n",vno, tempval.var[vno], antno);
}
*/

#if (constrained ==0)
(ant[antno]).ofn = objective_function(tempval.var );
#endif

#if (constrained ==1)
(ant[antno]).ofn = objective_function(tempval.var ,&((ant[antno]).penalty));
#endif


tempval.ofn = (ant[antno]).ofn;

/*printf("tempval.ofn = %f ant[%d].ofn = %f \n\n",tempval.ofn,antno, (ant[antno]).ofn);*/

if((ant[antno]).ofn < (ant[antno]).bestofn )
{
 ant[antno].bestofn = (ant[antno]).ofn;
 ant[antno].iter = iterationno;
}


if((ant[antno]).ofn < itbestant.ofn)
{
copy_val(&tempval , &itbestval);
copy_ant(&ant[antno],&itbestant);
}

if((ant[antno]).ofn >itworstant.ofn)
{
copy_val(&tempval, &itworstval);
copy_ant(&ant[antno], &itworstant);
}
}//endof antno loop

}//end of analysis

void decode_varvalues(float * varvalues,int antno)
{
int vno , valueno;
for(vno=0;vno<nov;vno++)
{
 if((ant[antno]).variable[vno] == ((var[vno]).nvalues - 1 ))
  {
  varvalues[vno] = var[vno].ulimit;
  }
  else
  {
    varvalues[vno] = (  (var[vno]).llimit + (((ant[antno]).variable[vno])* (var[vno]).increment ));
  }
}//end of vno loop

}//end of decode varvaules
