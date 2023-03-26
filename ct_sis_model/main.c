#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "graphs.h"
#include "kernel.h"

int main() {
    int N = 1000;
    double rho = 1.0;
    double alpha = 1.0;
    double gamma = 1.0;
    double initial_ratio = 0.5; 
    double T = 1000;
    unsigned long int seed = 398479375;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int nnode,nedge;

    int* edges = kinship_graphs_generator(&nnode, &nedge, N, rho, rng);

    int* sigma = (int*)malloc(sizeof(int)*nnode);
    while(1){
    for(int i=0; i<nnode; i++) {
        if(initial_ratio>gsl_rng_uniform_pos(rng)) {
            sigma[i]=1;
        } else {
            sigma[i]=0;
        }
    }

    double current_time=0;
    while(current_time<T) {
        double dt=0;
        kernel_uniform(alpha, gamma, &dt, nnode, nedge, sigma, edges, rng);

        current_time += dt;
    }

    double infected_ratio=0;
    for(int i=0; i<nnode; i++) {
        infected_ratio += sigma[i];
    }
    infected_ratio = infected_ratio/nnode;
    printf("%.10f \n", infected_ratio);
    }

    free(sigma);
    free(edges);
    gsl_rng_free(rng);
    return 0;
}
