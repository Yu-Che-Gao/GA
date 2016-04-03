#include"dec.h"
#include"func.h"

void initialize_trail()
{
int vno , valueno;

for(vno=0;vno<nov;vno++)
{
for(valueno=0;valueno<(var[vno]).nvalues;valueno++)
 {
 (var[vno]).trail[valueno] = tau0;
  }
}
}//end of trail initialization

void evaporate_trail()
{
int vno , valueno;

for(vno=0;vno<nov;vno++)
{
for(valueno=0;valueno<(var[vno]).nvalues;valueno++)
 {
 (var[vno]).trail[valueno] *= (1-rho);
 
}

}
}//end of evaporate trail

void trail()
{

evaporate_trail();

#if (ras ==1)
ras_update();
#endif

#if (as ==1)
as_update();
#endif

#if(acs ==1)
acs_update();
#endif

#if (lbt==1)
update_trail(&lbest[antno]);
#endif


#if (elitist ==1)
eas_update();
#endif

#if (taco == 1)
taco_update();
#endif

}//end of trail


void taco_update()
{
int tour_select();
int antno, updateantno;
for(antno=0;antno<ANTS;antno++)
{
updateantno = tour_select();
update_trail(&ant[updateantno]); 
}

update_trail_weighted(&bestant, elitist_ants);
}


void acs_update()
{

	update_trail_weighted( &bestant, 1 );

}


void as_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
*/
{
int antno;

    for ( antno = 0 ; antno < ANTS ; antno++ )
	update_trail( &ant[antno] );
}



void eas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Elitist Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants plus elitist ant deposit pheromones on matrix "pheromone"
*/
{
    int   antno;

    for ( antno = 0 ; antno < ANTS ; antno++ )
	update_trail( &ant[antno] );
    update_trail_weighted(&bestant , elitist_ants );
}




void ras_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is 
                      anyway not critical w.r.t. CPU time given that ras_ranks is 
		      typically very small.
*/
{
    int i,antno , target;
    float help_b[ANTS],b;

    for ( antno = 0 ; antno < ANTS ; antno++ )
	help_b[antno] = ant[antno].ofn;

    for ( i = 0 ; i < ras_ranks-1 ; i++ ) 
       {
	b = help_b[0]; target = 0;
	for ( antno = 0 ; antno <ANTS ; antno++ ) 
	{
	    if ( help_b[antno] < b ) 
	    {
		b = help_b[antno]; target = antno;
	    }
	}
	help_b[target] = 1e25;
	update_trail_weighted( &ant[target], ras_ranks-i-1 );
      }
    update_trail_weighted( &bestant, ras_ranks );
}



void update_trail(struct ant_struct * antptr)
{
int vno ,valueno;

for(vno=0;vno<nov;vno++)
{
  valueno= antptr->variable[vno];
if(antptr->ofn != 0.0)
{
 //(var[vno]).trail[valueno] +=  /*bestant.ofn*/( 1e-25+antptr->bestofn  / antptr->ofn);
 (var[vno]).trail[valueno] += 1; 
}
else
 (var[vno]).trail[valueno] +=  1.0;
}
}//end of update trail


void update_trail_weighted(struct ant_struct * antptr, int weight)
{
int vno ,valueno;

for(vno=0;vno<nov;vno++)
{
  valueno= antptr->variable[vno];
if(antptr->ofn != 0.0)
{
 //(var[vno]).trail[valueno] += (weight) *(1e-25+ antptr->bestofn)/( antptr->ofn);
 (var[vno]).trail[valueno] += (weight);
}
else
 (var[vno]).trail[valueno] += (weight) ;

}
}//end of update trail
