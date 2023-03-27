#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>

#include "graphs.h"
#include "kernel.h"

int main(int argc, char* argv[]) {
    int N = 100;
    double rho = 1.0;
    double alpha = 1.0;
    double gamma = 1.0;
    double initial_ratio = 0.5; 
    double T = 100;
    unsigned long int seed = 0;


    int opt;
    while ((opt = getopt(argc, argv, "N:r:a:g:i:T:s:")) != -1) { // Fixed option string
        switch (opt) {
            case 'N':
                N = atoi(optarg);
                break;
            case 'r':
                rho = atof(optarg);
                break;
            case 'a':
                alpha = atof(optarg);
                break;
            case 'g':
                gamma = atof(optarg);
                break;
            case 'i':
                initial_ratio = atof(optarg);
                break;
            case 'T':
                T = atof(optarg);
                break;
            case 's':
                seed = strtoul(optarg, NULL, 0);
                break;
            default:
                fprintf(stderr, "Usage: %s [-N n] [-r rho] [-a alpha] [-g gamma] [-i initial_ratio] [-T T] [-s seed]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    while(1) {
    int nnode,nedge;

    int* edges = kinship_graphs_generator(&nnode, &nedge, N, rho, rng);
    double beta = ((double)(nedge*2))/nnode;

    int* sigma   = (int*)malloc(sizeof(int)*nnode);
    int* propose = (int*)malloc(sizeof(int)*nnode);


    clock_t start, end;
    double cpu_time_used;

    start = clock();
    for(int i=0; i<nnode; i++) {
        if(initial_ratio>gsl_rng_uniform_pos(rng)) {
            sigma[i]=1;
        } else {
            sigma[i]=0;
        }
    }

    int np=0;
    double current_time=0;
    while(current_time<T) {
        double dt=0;
        np = kernel_uniform_scan(alpha, gamma, &dt, nnode, nedge, sigma, edges, propose, rng);

        current_time += dt;

        if(current_time<T) kernel_update(np, propose, sigma);
    }

    double infected_ratio=0;
    for(int i=0; i<nnode; i++) {
        infected_ratio += sigma[i];
    }
    infected_ratio = infected_ratio/nnode;

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("N=%d rho=%.4f beta=%.6f alpha=%.4f if=%.10f t=%f\n", N,rho,beta,alpha,infected_ratio,cpu_time_used);

    free(sigma);
    free(propose);
    free(edges);
    kernel_free();
    }

    gsl_rng_free(rng);
    return 0;
}
