
#include"dec.h"
#include"func.h"
/*----------------------------------------------------------*/
/* Files for tournament selection :                         */
/* Source : sga.c (c) E.Goldberg                            */
/*----------------------------------------------------------*/

void preselect_tour()
{
    reset1();
    tourneypos = 0;
}


int tour_select()
{
    int pick, winner, i;

    /* If remaining members not enough for a tournament, then reset list */
start_select :
    if((ANTS - tourneypos) < tourneysize)
    {
        reset1();
        tourneypos = 0;
    }

    /* Select tourneysize structures at random and conduct a tournament */
    winner=tourneylist[tourneypos];
/* Added by RBA */
    if( winner < 0 || winner > ANTS-1) {
                                             printf("\n Warning !! ERROR1");
                                             printf(" tourpos = %d",tourneypos);
                                             printf(" winner = %d",winner);
                                             preselect_tour();
                                             goto start_select; }
    for(i=1; i<tourneysize; i++)
    {
        pick=tourneylist[i+tourneypos];
/* Added by RBA */
        if (pick < 0 || pick > ANTS-1) { preselect_tour();
                                             printf("\n Warning !! ERROR2");
                                             goto start_select; }
	// case 1:
       #if (constrained == 1)
	if (ant[winner].penalty > ant[pick].penalty) winner = pick;
	else if ((ant[winner].penalty <= 0.0) && (ant[pick].penalty <= 0.0))
	  {
	      if( ant[pick].ofn <ant[winner].ofn) winner=pick;
	  }
       #endif
      
      #if (costrained ==0)
	      if( ant[pick].ofn <  ant[winner].ofn) winner=pick;
	#endif
    }
    
    /* Update tourneypos */
    tourneypos += tourneysize;
    return(winner);
}


#if (sharing==1)
double distanc(one,two)
     int one,two;
{
  int k;
  double sum;
  
  sum = 0.0;
  for (k=0; k<nvar_bin; k++) 
    sum += square((oldpop[one].xbin[k]-oldpop[two].xbin[k])/(xbin_upper[k]-xbin_lower[k]));
  for (k=0; k<nvar_real; k++) 
    sum += square((oldpop[one].xreal[k]-oldpop[two].xreal[k])/(xreal_upper[k]-xreal_lower[k]));

  return (sqrt(sum/(nvar_bin + nvar_real)));
};

int tour_select_constr() 
{ 
  int pick, winner, i, minus=0, rand_pick, rand_indv, flag, indv; 
  
  /* If remaining members not enough for a tournament, then reset list */ 
  start_select : 
    if((ANTS - tourneypos) < tourneysize) 
      { 
        reset1(); 
        tourneypos = 0; 
      } 
  
  /* Select tourneysize structures at random and conduct a tournament */ 
  winner = tourneylist[tourneypos]; 
  /* Added by RBA */ 
  if( winner < 0 || winner > ANTS-1) 
    {
      printf("\n Warning !! ERROR1"); 
      printf(" tourpos = %d",tourneypos); 
      printf(" winner = %d",winner); 
      preselect_tour(); 
      goto start_select; 
    } 
  for(i=1; i<tourneysize; i++) 
    { 
      pick = tourneylist[i+tourneypos];

      if((oldpop[winner].penalty>0.0) && (oldpop[pick].penalty<=0.0))
	winner = pick;
      else if((oldpop[winner].penalty>0.0)&&(oldpop[pick].penalty>0.0))
	{	
	  if(oldpop[pick].penalty < oldpop[winner].penalty) 
	    winner=pick;
	}
      else if((oldpop[winner].penalty<=0.0)&&(oldpop[pick].penalty<=0.0))
	{
	  if(distanc(winner,pick) < sigma_share) 
	    {
	      if (MINM * oldpop[pick].ofn < MINM * oldpop[winner].ofn) 
		winner=pick;
	    }
	  else
	    {
	      minus = -1;
	      for (indv = flag = 0; indv<critical_size && flag==0; indv++) 
		{
		  rand_indv = rnd(0,ANTS-1);
		  rand_pick = tourneylist[rand_indv];
		  if(oldpop[rand_pick].penalty <= 0.0) 
		    {
		      if(distanc(winner,rand_pick) < sigma_share) 
			{
			  flag = 1;
			  if (MINM * oldpop[rand_pick].ofn < MINM * oldpop[winner].ofn) 
			    winner=rand_pick;
			}
		    }
		}
	    }
	}
      if (pick < 0 || pick > ANTS-1) 
	{ 
	  preselect_tour(); 
	  printf("\n Warning !! ERROR2"); 
	  goto start_select; 
	} 
    } 
    /* Update tourneypos */ 
  tourneypos += tourneysize + minus; 
  return(winner); 
}

#endif

void reset1()
/* Name changed from reset because of clash with lib. function - RBA */
/* Shuffles the tourneylist at random */
{
    int i, rand1, rand2, temp_site;

    for(i=0; i<ANTS; i++) 

    for(i=0; i < ANTS; i++)
    {
        rand1= rnd(0,ANTS-1);
        rand2=  rnd(0,ANTS-1);
        temp_site = tourneylist[rand1];
        tourneylist[rand1]=tourneylist[rand2];
        tourneylist[rand2]=temp_site;
    }
}

