#include"rand.h"
long int seed;

double ran01( long *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

int rnd(int low,int high)
/* Pick a random integer between low and high */
{
    int i;

    if(low >= high)
        i = low;
    else
    {
        i = (int)((ran01(&seed) * (high - low + 1)) + low);
        if(i > high) i = high;
    }
    return(i);
}

double rndreal(double lo ,double hi)
/* real random number between specified limits */
{
    return((ran01(&seed) * (hi - lo)) + lo);
}

