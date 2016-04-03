void as_update( void );
void acs_update();
void eas_update( void );
void taco_update();
void ras_update();
void swap(long int v[], long int i, long int j);
void sort(long int v[], long int left, long int right);
void preselect_tour();
int tour_select();
double distanc(int one,int two);
int tour_select_constr(); 
void reset1();
void ls(struct values *gvalptr, int );
void find_lsbestant_valueno(struct values *lsvalptr, struct ant_struct *lsantptr);
void init_out(char *);
void lsstats(int,int  );
void read_data(char *);
void update_trail_weighted(struct ant_struct * antptr, int weight);
void initialize_trail();
void initialize_ants_variables(int runno);
//void lsinitialize_ants_vars(int itno  );
int find_var_valueno(int vno);
int find_no_of_max(int vno,int *noofmax);
void find_prob(int vno);
void decode_varvalues( float *, int );

#if (constrained ==1)
float objective_function(float *varvalues, float * antpenalty);
#endif
#if (constrained ==0)
float objective_function(float *varvalues);
#endif
void copy_ant(struct ant_struct *from , struct ant_struct *to);
void iteration_init(int iterationno);
void find_values();
void copy_val(struct values  *from ,struct values *to);
void analysis(int iterationno);
void update_stats(int iterationno,int runno);
void print_itstats(int iterationno, int runno);
void trail();
void report_run(int runno);
void final_stats();
void find_nvalue(int vno);
void init_ant(struct ant_struct *antptr, int iterationno, float iofn);
void evaporate_trail();
void update_trail(struct ant_struct * antptr);
void print_ant(struct ant_struct * antptr);
void print_data();
double std_deviation_ants(  double d_mean ) ;
float find_mean();
double std_deviationr( double *values, long int max, double d_mean );
double std_deviation( long int *values, long int max, double i_mean ) ;

#include"rand.h"
