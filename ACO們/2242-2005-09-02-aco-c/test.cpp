#include<stdio.h>
#include<math.h>
#include"parameters.h"




#if (constrained ==1)
float objective_function(float *varvalues, float * antpenalty);
#endif
#if (constrained ==0)
float objective_function(float *varvalues);
#endif

long nevaluations;



main()
{
//	nevaluations = 0;

		float values[nov] ;
		float cost , penalty;
		values[0] = 1  ;
		values[1] = 3  ;
		penalty = 0.0;
		cost = objective_function( values, &penalty);
		printf(" ofn = %15.15f , sumconstr = %6.8f ", cost,penalty);

}
