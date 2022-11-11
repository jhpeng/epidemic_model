#include <gsl/gsl_rng.h>

#ifndef time_evolution_h
#define time_evolution_h

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

double total_infected_time_value();

double total_ninfection_value();

double total_nrecover_value();

#endif
