#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "graphs.h"

int main() {
    int N=100;
    double rho=1.0;
    unsigned long int seed = 398479375;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int nnode,nedge;

    for(int i=0;i<1000;i++){
        int* edges = kinship_graphs_generator(&nnode, &nedge, N, rho, rng);

        printf("nnode=%d, nedge=%d \n",nnode,nedge);
        free(edges);
    }
    
    gsl_rng_free(rng);
    return 0;
}
