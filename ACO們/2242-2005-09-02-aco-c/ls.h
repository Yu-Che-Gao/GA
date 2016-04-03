/*#define LSANTS 5	
#define lsnvalues 2
#define lsncmax 5
#define lsrho 0.2
#define lsalpha 1
#define lstau0 1
#define lsq0 0.7*/
#include<string.h>
struct variable *lsvar;
struct ant_struct *lsant, *lsbestant, *lsitbestant, *lsitworstant;
struct values *lsitbestval, *lsitworstval,*lsbestval;


void lsfreememory();
void lsprint_itstats(int iterationno);
void lsupdate_stats(int iterationno);
void lsanalysis(int iterationno);
void lsdecode_varvalues(float * varvalues,int antno);
void lsiteration_init(int iterationno);
void lsget_data(struct values *gvalptr);
void lsprint_data();
void lsfind_nvalue(int vno);
void lsfind_values();
int lsfind_var_valueno(int vno);
int lsfind_no_of_max(int vno,int *noofmax);
void lsfind_prob(int vno);
void lsevaporate_trail();
void lstrail();
void lsupdate_trail(struct ant_struct * antptr);
void lsinitialize_trail();

