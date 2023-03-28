#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>

#include "graphs.h"
#include "kernel.h"

void quasistationary_state_sampling(const int* edges, int nnode, int nedge, double alpha, double gamma, double t, int nsample, gsl_rng* rng) {
    int* sigma     = (int*)malloc(sizeof(int)*nnode);
    int* sigma_p   = (int*)malloc(sizeof(int)*nnode);
    int* propose   = (int*)malloc(sizeof(int)*nnode);

    for(int i=0; i<nnode; i++) {
        sigma[i]=1;
        sigma_p[i]=1;
    }

    FILE* if_file = fopen("infected_number.dat","a");
    FILE* if2_file = fopen("infected_number_2.dat","a");
    FILE* cput_file = fopen("cpu_time_used.dat","a");
    FILE* nt_file = fopen("number_trial_run.dat","a");


    clock_t start, end;
    double cpu_time_used;
    int ntrial=0;

    start = clock();
    for(int i_sample=0; i_sample<nsample; i_sample++) {
        int np=0;
        double dt=0;
        double current_time=0;
        while(current_time<t) {
            np = kernel_uniform_scan(alpha, gamma, &dt, nnode, nedge, sigma, edges, propose, rng);

            current_time += dt;

            if(current_time<t) kernel_update(np, propose, sigma);
        }

        int infected_number=0;
        for(int i=0; i<nnode; i++) {
            infected_number += sigma[i];
        }

        ntrial++;
        if(infected_number==0) {
            for(int i=0; i<nnode; i++) {
                sigma[i] = sigma_p[i];
            }
            i_sample--;
        } else {
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            double inf2 = ((double)infected_number)*((double)infected_number);

            printf("if=%d if2=%.12e t=%f ntrial=%d\n", infected_number, inf2, cpu_time_used, ntrial);

            fprintf(if_file,"%d \n", infected_number);
            fprintf(if2_file,"%.12e \n", inf2);
            fprintf(cput_file,"%.4f \n", cpu_time_used);
            fprintf(nt_file,"%d \n", ntrial);

            start = clock();
            ntrial=0;
        }
    }

    fclose(if_file);
    fclose(if2_file);
    fclose(cput_file);
    fclose(nt_file);


    free(sigma);
    free(sigma_p);
    free(propose);
    kernel_free();
}

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
                fprintf(stderr, "Usage: %s [-N n] [-r rho] [-a alpha] [-g gamma] [-i initial_ratio] [-T T] [-n nsample] [-s seed]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int nnode,nedge;

    int* edges = kinship_graphs_generator(&nnode, &nedge, N, rho, rng);
    double beta = ((double)(nedge*2))/nnode;

    quasistationary_state_sampling(edges, nnode, nedge, alpha, gamma, T, nsample, rng) ;

    FILE* beta_file = (FILE*)fopen("ave_number_edge_per_node.dat","a");
    fprintf(beta_file, "%.12f \n", beta);
    fclose(beta_file);

    FILE* edge_file = (FILE*)fopen("edges.dat","a");
    for(int i=0; i<nedge; i++) {
        fprintf(edge_file, "%d %d \n", edges[i*2],edges[i*2+1]);
    }
    fclose(edge_file);

    free(edges);
    gsl_rng_free(rng);
    return 0;
}
