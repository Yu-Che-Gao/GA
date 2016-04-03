#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#include"parameters.h"
#define pi 3.14
#define ks 10
#define ke 500


#if (constrained ==1)
float objective_function(float *varvalues, float * antpenalty);
#endif
#if (constrained ==0)
float objective_function(float *varvalues);
#endif


extern long int nevaluations;

#define var10 0
#define ex1 1

#if(constrained==1)
float objective_function( float *varvalues, float *antpenalty)
/* takes array of varvalues of length nov and
 * returns objective function which is a float values*/
{
	extern long int nevaluations;
			nevaluations++;
float ofn;

float objfn, x1,x2 ;
x1 = varvalues[0];
x2 = varvalues[1];

# if(ex1)
objfn = 1 +( (x1*x1 +x2 -11) *(x1*x1 +x2 -11) ) + ( (x1 +x2*x2 - 7) * (x1
			+ x2*x2 -7));

float g1,g2 ; 
g1 = 1.0 - (( pow((x1 -5.0),2) + pow((x2),2) ) / 26.0 ) ;
g2 = 1.0 - ( ( 4 * x1 + x2 ) / 20.0) ;

		float sumconstr = 0.0, penaltyconst = 1;
		if(g1 >=0)
		{
		sumconstr += g1; 
		}
		if(g2 >=0)
		{
		sumconstr += g1; 
		}
	*antpenalty = penaltyconst * sumconstr;

	return objfn * ( 1 + *antpenalty) ;

#endif
		//printf("x1= %f x2  is %g \n",x1,x2);
		//printf("sumconstr= %f g1 is %g \n",sumconstr,g1);
		//return objfn ;

#if(var10)
float constr[7]= {0.0,0.0,0.0,0.0,0.0,0.0}, sumconstr=0.0;
float e,a,b,c1,c2,d1,d2,g1,g2,h,i;
//float x[]  =   { 1.000 , 6.000 , 0.000 , 0.000 , 0.000 , 8.000 , 1.000 , 5.000 , 0.000 , 6.000 , 5.000}  ;
// float x[]  =    { 1, 6, 0, 0,0, 8, 1,5,0,6,5};


a=x[0];
b=x[1];
c1=x[2];
c2=x[3];
d1=x[4];
d2=x[5];
g1=x[6];
g2=x[7];
h=x[8];
i=x[9];
e= 20- (a+b+c2+d2);


constr[0] = ((b-c1)/7 )-1;
constr[1] = ((b+c2-d1)/6) -1;
constr[2] = ((b+c2 - g1)/5) -1;
constr[3] = ((d2 -g2 )/3) -1;
constr[4] = ((c2+d2 -h -i )/2) -1;
constr[5] = (e/5.0) -1  ;
constr[6] =  -e   ;
//printf("e = %f\n",e);
	int cno;
	for(cno=0;cno<7;cno++)
	{
		if(constr[cno] >=0)
		{
		sumconstr += constr[cno];
		}
		//printf("constraint [%d] = %f \n",cno,constr[cno]);
	}
sumconstr *= ks;

ofn = (60*a +10*b + 20*c1 +20*c2 + 5*d1+ 5*d2 + 30*e + 15*g1 +15*g2 + 12*h +6*i)* ( 1 + sumconstr ); 


*antpenalty = sumconstr;

//printf(" ofn = %15.15f , sumconstr = %6.8f \n", ofn,sumconstr);   


return ofn;

#endif
}




#endif





#define ex2 1
#define f2 0
#define f3 0
#define f1 0


#if(constrained==0)
float objective_function( float *varvalues )
/* takes array of varvalues of length nov and
 * returns objective function which is a float values*/
{
nevaluations++;

float x1,x2,x3 , x4,  objfn ;
int i;
x1 = varvalues[0];
x2 = varvalues[1];
//x3 = varvalues[2];
//x4 = varvalues[3];



if(ex2)
	{
		objfn = 1 +( (x1*x1 +x2 -11) *(x1*x1 +x2 -11) ) + ( (x1 +x2*x2 - 7) * (x1
					+ x2*x2 -7));
//	objfn = 1 +( (x1*x1 +x2 -11) *(x1*x1 +x2 -11) ) + ( (x1 +x2*x2  - 7) * (x1 + x2*x2  -7));
	}

else if(f2)
	objfn = 1 + 100* (x1*x1 -x2)*(x1*x1 -x2) + (1-x2)*(1-x2) ;

else if(f3)
	{
	objfn = 50  ;
	for(i=0;i<5;i++)
		{
			objfn += (varvalues[i] * varvalues[i] - 10*cos(2 * pi *varvalues[i]));
		}
	}
else if(f1)
	{
	objfn = 0  ;
	for(i=0;i<3;i++)
		{
		objfn += (varvalues[i] * varvalues[i]) ;
		}
	}




return objfn ;
}
#endif

