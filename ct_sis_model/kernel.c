#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>


int random_integer_with_distribution(const gsl_rng *rng, double *weight_dist, double total_weight, int len) {
    double u = gsl_rng_uniform(rng)*total_weight;
    double cdf;
    cdf = weight_dist[0];

    for (int i = 1; i < len; i++) {
        if (u <= cdf) {
            return i-1;
        }
        cdf = cdf + weight_dist[i];
    }

    return len-1;
}

static int* kernel_uniform_recover_list = NULL;
static int* kernel_uniform_infected_list = NULL;
int kernel_uniform_scan(double alpha, double gamma, double* dt, int nnode, int nedge, const int* sigma, const int* edges, int* propose, gsl_rng* rng) {
    int nr=0;
    int ni=0;
    int np=0;

    if(kernel_uniform_recover_list==NULL) {
        kernel_uniform_recover_list = (int*)malloc(sizeof(int)*nnode);
        if (kernel_uniform_recover_list == NULL) {
            fprintf(stderr, "Error: failed to allocate memory for kernel_uniform_recover_list.\n");
            exit(EXIT_FAILURE);
        }
    }

    if(kernel_uniform_infected_list==NULL) {
        kernel_uniform_infected_list = (int*)malloc(sizeof(int)*nedge);
        if (kernel_uniform_infected_list == NULL) {
            fprintf(stderr, "Error: failed to allocate memory for kernel_uniform_infected_list.\n");
            exit(EXIT_FAILURE);
        }
    }

    for(int i=0; i<nnode; i++) {
        if(sigma[i]) {
            kernel_uniform_recover_list[nr] = i;
            nr++;
        }
    }

    for(int i=0; i<nedge; i++) {
        int i_node = edges[i*2+0];
        int j_node = edges[i*2+1];

        if(sigma[i_node]^sigma[j_node]) {
            kernel_uniform_infected_list[ni] = i;
            ni++;
        }
    }

    if(nr==0) {
        *dt = DBL_MAX;
        return np;
    }

    if(nr*gamma > (nr*gamma+ni*alpha)*gsl_rng_uniform_pos(rng)) {
        int dis = nr*gsl_rng_uniform_pos(rng);
        int i_node = kernel_uniform_recover_list[dis];
        propose[np] = i_node;
        np++;
    } else {
        int dis = ni*gsl_rng_uniform_pos(rng);
        int i_edge = kernel_uniform_infected_list[dis];
        int i_node = edges[i_edge*2+0];
        int j_node = edges[i_edge*2+1];

        if(sigma[i_node]) {
            propose[np] = j_node;
        } else {
            propose[np] = i_node;
        }
        np++;
    }

    *dt = -log(gsl_rng_uniform_pos(rng))/(nr*gamma+ni*alpha);

    return np;
}

void kernel_update(int np, const int* propose, int* sigma) {
    for(int i=0;i<np;i++) {
        sigma[propose[i]] = sigma[propose[i]]^1;
    }
}

void kernel_free() {
    if(kernel_uniform_recover_list==NULL) {
        free(kernel_uniform_recover_list);
        kernel_uniform_recover_list = NULL;
    }

    if(kernel_uniform_infected_list==NULL) {
        free(kernel_uniform_infected_list);
        kernel_uniform_infected_list = NULL;
    }
}
