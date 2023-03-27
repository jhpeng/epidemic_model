#include <gsl/gsl_rng.h>

#ifndef kernel_h
#define kernel_h

int random_integer_with_distribution(const gsl_rng *rng, double *weight_dist, double total_weight, int len);

int kernel_uniform_scan(double alpha, double gamma, double* dt, int nnode, int nedge, const int* sigma, const int* edges, int* propose, gsl_rng* rng);

void kernel_update(int np, const int* propose, int* sigma);

void kernel_free();

#endif
