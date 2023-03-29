#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>

#include "graphs.h"
#include "kernel.h"

/**
 * @brief Calculates the mean of an array of doubles (static function).
 *
 * This static function computes the mean (average) value of a given sequence of doubles. If the length is less than or equal to zero,
 * the function returns zero.
 *
 * @note This function is static and can only be used within the translation unit (source file) it is defined in.
 *
 * @param sequence Pointer to an array of doubles representing the input sequence.
 * @param length The length of the input sequence.
 * @return The mean value of the input sequence, or 0.0 if the length is less than or equal to zero.
 */
static double mean(double* sequence, int length) {
    if (length <= 0) {
        return 0.0; 
    }

    double mean_value=0;
    for(int i=0; i<length; i++) {
        mean_value += sequence[i];
    }

    return mean_value/length;
}

/**
 * @brief Performs quasistationary state sampling for an epidemic process.
 *
 * This function performs quasistationary state sampling for an epidemic process, calculating the infected ratio,
 * infected ratio squared, and number of trials. It stores these values in separate output files.
 *
 * @param edges Pointer to an array of edge pairs representing the network connections.
 * @param nnode The total number of nodes in the network.
 * @param nedge The total number of edges in the network.
 * @param alpha The alpha parameter for the epidemic model.
 * @param gamma The gamma parameter for the epidemic model.
 * @param t The time threshold for the sampling process.
 * @param nsample The number of samples to be generated.
 * @param rng Pointer to a GSL random number generator.
 */

void quasistationary_state_sampling(const int* edges, int nnode, int nedge, double alpha, double gamma, double t, int nsample, gsl_rng* rng) {
    int* sigma     = (int*)malloc(sizeof(int)*nnode);
    int* sigma_p   = (int*)malloc(sizeof(int)*nnode);
    int* propose   = (int*)malloc(sizeof(int)*nnode);

    for(int i=0; i<nnode; i++) {
        sigma[i]=1;
        sigma_p[i]=1;
    }

    int length = 1000;
    double* infected_ratio_sequence = (double*)malloc(sizeof(double)*length);
    double* infected_ratio2_sequence = (double*)malloc(sizeof(double)*length);
    double* number_trial_sequence = (double*)malloc(sizeof(double)*length);


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
            double inf2 = ((double)infected_number)*((double)infected_number);

            infected_ratio_sequence[i_sample%length] = ((double)infected_number);
            infected_ratio2_sequence[i_sample%length] = inf2;
            number_trial_sequence[i_sample%length] = ntrial;

            if((i_sample+1)%length==0) {
                end = clock();
                cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

                double infected_ratio_mean = mean(infected_ratio_sequence,length)/nnode;
                double infected_ratio2_mean = mean(infected_ratio2_sequence,length)/nnode/nnode;
                double number_trial_mean = mean(number_trial_sequence,length);

                printf("inf=%.12e inf2=%.12e ntrial=%.12e cpu_time=%.4f\n", infected_ratio_mean, infected_ratio2_mean, number_trial_mean, cpu_time_used);

                FILE* if_file = fopen("infected_ratio.dat","a");
                FILE* if2_file = fopen("infected_ratio2.dat","a");
                FILE* cput_file = fopen("cpu_time_used.dat","a");
                FILE* nt_file = fopen("number_trial_run.dat","a");

                fprintf(if_file,"%.12e \n", infected_ratio_mean);
                fprintf(if2_file,"%.12e \n", infected_ratio2_mean);
                fprintf(cput_file,"%.4f \n", cpu_time_used);
                fprintf(nt_file,"%.12e \n", number_trial_mean);

                fclose(if_file);
                fclose(if2_file);
                fclose(cput_file);
                fclose(nt_file);

                start = clock();
            }

            ntrial=0;
        }
    }


    free(sigma);
    free(sigma_p);
    free(propose);
    free(infected_ratio_sequence);
    free(infected_ratio2_sequence);
    free(number_trial_sequence);
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
