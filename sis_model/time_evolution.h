#include <gsl/gsl_rng.h>

#ifndef time_evolution_h
#define time_evolution_h

double ninfection_ave_value();

double nrecover_ave_value();

void ninfection_count_plus_one();

void nrecover_count_plus_one();

void print_ninfection();

void print_nrecover();

void kernel_init(int nnode);

void kernel_linear(double alpha, double gamma, double dt, int nnode, int nedge, int* sigma, int* edges, gsl_rng* rng);

/** alpha : transmition rate
*** gamma : recover rate
***    dt : time interval
*** nnode : number of node
*** nedge : number of edge
*** sigma : status of node
*** edges : list of edge
***   rng : random number generator
**/

void kernel_close(int nnode, int* sigma);

double kernel_total_infected_time(double dt);

#endif
