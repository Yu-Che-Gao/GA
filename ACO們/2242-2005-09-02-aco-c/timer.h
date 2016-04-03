
extern clock_t start_time;

extern double elapsed;

int time_expired();

void start_timers();


typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;

extern double elapsed_time(TIMER_TYPE type);

extern double elapsed;

