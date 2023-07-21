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

void qs_sampling_memory(const int* edges, int nnode, int nedge, double alpha, double gamma, int nsample, gsl_rng* rng) {
    int* sigma = (int*)malloc(sizeof(int)*nnode);
    int* propose = (int*)malloc(sizeof(int)*nnode);
    int* distrib = (int*)malloc(sizeof(int)*nnode);

    for(int i=0;i<nnode;i++) {
        if(gsl_rng_uniform_pos(rng)<0.5) {
            sigma[i] = 1;
        } else {
            sigma[i] = 0;
        }
        distrib[i] = 0;
    }

    int M = 100;
    int* sigma_history = (int*)malloc(sizeof(int)*nnode*M);

    for(int i=0;i<(nnode*M);i++) {
        if(gsl_rng_uniform_pos(rng)<0.5) {
            sigma_history[i] = 1;
        } else {
            sigma_history[i] = 0;
        }
    }

    clock_t start, end;
    double cpu_time_used;

    start = clock();
    printf("start of thermalization...\n");

    double pr = 0.01;
    double dt = 0.0;
    int np = 0;
    int nreject = 0;
    int nthermal = M*1000;
    for(int i_thermal=0; i_thermal<nthermal;) {
        dt = 0;
        while(dt*pr<gsl_rng_uniform_pos(rng)) {
            np = kernel_uniform_scan(alpha, gamma, &dt, nnode, nedge, sigma, edges, propose, rng);
            kernel_update(np, propose, sigma);
        }

        int infected_number=0;
        for(int i=0; i<nnode; i++) {
            infected_number += sigma[i];
        }

        if(infected_number==0) {
            int ihist = (int)(gsl_rng_uniform_pos(rng)*M);
            for(int i=0;i<nnode;i++) {
                sigma[i] = sigma_history[nnode*ihist+i];
            }
            nreject++;
        } else {
            int ihist = i_thermal%M;
            for(int i=0;i<nnode;i++) {
                sigma_history[nnode*ihist+i] = sigma[i];
            }

            if((i_thermal+1)%100==0) {
                double percentage = (double)(i_thermal+1)*100/nthermal;

                end = clock();
                cpu_time_used = (double)(end-start)/CLOCKS_PER_SEC;

                double est_time = 0;
                if(percentage<0.1) {
                    est_time = 0;
                } else {
                    est_time = cpu_time_used/percentage*(100-percentage);
                }
                
                printf("%lf ( %d %d) est: %.2f s\n",percentage, nreject, infected_number, est_time);
            }
            i_thermal++;
            nreject=0;
        }
    }

    end = clock();
    cpu_time_used = (double)(end-start)/CLOCKS_PER_SEC;
    printf("end of thermalization...\n");
    printf("cpu time used : %lf s\n", cpu_time_used);
    start = clock();

    pr = 0.01;
    printf("start sampling...\n");
    nreject=0;
    for(int i_sample=0; i_sample<nsample;) {
        dt = 0;
        while(dt*pr<gsl_rng_uniform_pos(rng)) {
            np = kernel_uniform_scan(alpha, gamma, &dt, nnode, nedge, sigma, edges, propose, rng);
            kernel_update(np, propose, sigma);
        }

        int infected_number=0;
        for(int i=0; i<nnode; i++) {
            infected_number += sigma[i];
        }

        if(infected_number==0) {
            int ihist = (int)(gsl_rng_uniform_pos(rng)*M);
            for(int i=0;i<nnode;i++) {
                sigma[i] = sigma_history[nnode*ihist+i];
            }
            nreject++;
        } else {
            int ihist = (int)(gsl_rng_uniform_pos(rng)*M);
            for(int i=0;i<nnode;i++) {
                sigma_history[nnode*ihist+i] = sigma[i];
                distrib[i] += sigma[i];
            }

            end = clock();
            cpu_time_used = (double)(end-start)/CLOCKS_PER_SEC;

            FILE* sampling_file = fopen("sampling_data.txt","a");
            printf("%d %d %d %lf \n", i_sample, nreject, infected_number, cpu_time_used);
            fprintf(sampling_file, "%d %d %lf \n", nreject, infected_number, cpu_time_used);

            fclose(sampling_file);

            start = clock();

            i_sample++;
            nreject=0;
        }
    }

    FILE* distrib_file = fopen("distrib_data.txt","w");
    for(int i=0;i<nnode;i++) {
        fprintf(distrib_file,"%d ",distrib[i]);
    }
    fprintf(distrib_file,"\n");
    fclose(distrib_file);
    
}

void qs_sampling_largeTlimit(const int* edges, int nnode, int nedge, double alpha, double gamma, double ir, double t, int nsample, gsl_rng* rng) {
    int* sigma   = (int*)malloc(sizeof(int)*nnode);
    int* propose = (int*)malloc(sizeof(int)*nnode);

    for(int i=0; i<nnode; i++) {
        if(gsl_rng_uniform_pos(rng)<ir) {
            sigma[i] = 1;
        } else {
            sigma[i] = 0;
        }
    }


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
                if(gsl_rng_uniform_pos(rng)<ir) {
                    sigma[i] = 1;
                } else {
                    sigma[i] = 0;
                }
            }
            i_sample--;
        } else {
            end = clock();
            cpu_time_used = (double)(end-start)/CLOCKS_PER_SEC;

            FILE* sampling_file = fopen("sampling_data.txt","a");
            printf("%d %d %d %lf \n", i_sample, ntrial, infected_number, cpu_time_used);
            fprintf(sampling_file, "%d %d %lf \n", ntrial, infected_number, cpu_time_used);

            fclose(sampling_file);

            start = clock();
            ntrial=0;
        }
    }

    free(sigma);
    kernel_free();
}

int main(int argc, char* argv[]) {
    int N = 100;
    int mode = 1;
    double rho = 1.0;
    double alpha = 1.0;
    double gamma = 1.0;
    double initial_ratio = 0.5; 
    double T = 100;
    int nsample = 1000;
    unsigned long int seed = 0;


    int opt;
    while ((opt = getopt(argc, argv, "m:N:r:a:g:i:T:n:s:")) != -1) { // Fixed option string
        switch (opt) {
            case 'm':
                mode = atoi(optarg);
                break;
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

    if(mode) {
        qs_sampling_largeTlimit(edges, nnode, nedge, alpha, gamma, initial_ratio, T, nsample, rng) ;
    } else {
        qs_sampling_memory(edges, nnode, nedge, alpha, gamma, nsample, rng);
    }

    FILE* beta_file = (FILE*)fopen("ave_number_edge_per_node.dat","a");
    fprintf(beta_file, "%.12f \n", beta);
    fclose(beta_file);

    FILE* edge_file = (FILE*)fopen("edges.dat","w");
    for(int i=0; i<nedge; i++) {
        fprintf(edge_file, "%d %d \n", edges[i*2],edges[i*2+1]);
    }
    fclose(edge_file);

    free(edges);
    gsl_rng_free(rng);
    return 0;
}
