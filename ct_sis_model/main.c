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
    int nsample = 1000;
    unsigned long int seed = 0;


    int opt;
    while ((opt = getopt(argc, argv, "N:r:a:g:i:T:n:s:")) != -1) { // Fixed option string
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
            case 'n':
                nsample = atoi(optarg);
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

    for(int i_sample=0; i_sample<nsample; i_sample++) {
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

    int infected_number=0;
    for(int i=0; i<nnode; i++) {
        infected_number += sigma[i];
    }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("N=%d rho=%.4f beta=%.6f alpha=%.4f if=%d t=%f\n", N,rho,beta,alpha,infected_number,cpu_time_used);

    FILE* beta_file = fopen("beta.dat","a");
    FILE* if_file = fopen("infected_number.dat","a");
    FILE* cput_file = fopen("cpu_time_used.dat","a");

    fprintf(beta_file,"%.12f \n", beta);
    fprintf(if_file,"%d \n", infected_number);
    fprintf(cput_file,"%.4f \n", cpu_time_used);

    fclose(beta_file);
    fclose(if_file);
    fclose(cput_file);

    free(sigma);
    free(propose);
    free(edges);
    kernel_free();
    }

    gsl_rng_free(rng);
    return 0;
}
