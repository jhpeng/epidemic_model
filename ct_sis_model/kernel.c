#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
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

int* kernel_uniform_recover_list = NULL;
int* kernel_uniform_infected_list = NULL;
void kernel_uniform(double alpha, double gamma, double* dt, int nnode, int nedge, int* sigma, int* edges, gsl_rng* rng) {
    int nr=0;
    int ni=0;

    if(kernel_uniform_recover_list==NULL) {
        kernel_uniform_recover_list = (int*)malloc(sizeof(int)*nnode);
    }

    if(kernel_uniform_infected_list==NULL) {
        kernel_uniform_infected_list = (int*)malloc(sizeof(int)*nedge);
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
        *dt = INT_MAX;
        return;
    }

    if(nr*gamma > (nr*gamma+ni*alpha)*gsl_rng_uniform_pos(rng)) {
        int dis = nr*gsl_rng_uniform_pos(rng);
        int i_node = kernel_uniform_recover_list[dis];
        sigma[i_node] = 0;
    } else {
        int dis = ni*gsl_rng_uniform_pos(rng);
        int i_edge = kernel_uniform_infected_list[dis];
        int i_node = edges[i_edge*2+0];
        int j_node = edges[i_edge*2+1];

        sigma[i_node] = 1;
        sigma[j_node] = 1;
    }

    *dt = -log(gsl_rng_uniform_pos(rng))/(nr*gamma+ni*alpha);

    return;
}
