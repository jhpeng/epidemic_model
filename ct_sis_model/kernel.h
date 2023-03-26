#include <gsl/gsl_rng.h>

#ifndef kernel_h
#define kernel_h

int random_integer_with_distribution(const gsl_rng *rng, double *weight_dist, double total_weight, int len);

void kernel_uniform(double alpha, double gamma, double* dt, int nnode, int nedge, int* sigma, int* edges, gsl_rng* rng);

#endif
