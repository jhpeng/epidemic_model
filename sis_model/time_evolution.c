#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>


void kernel_linear(double alpha, double gamma, double dt, int nnode, int nedge, int* sigma, int* edges, gsl_rng* rng) {
    double pt = alpha*dt;
    double pr = gamma*dt;
    int i_node,j_node,i_sigma,j_sigma;

    for(int i_edge=0;i_edge<nedge;i_edge++) {
        i_node = edges[i_edge*2+0];
        j_node = edges[i_edge*2+1];

        i_sigma = sigma[i_node];
        j_sigma = sigma[j_node];
        if(i_sigma^j_sigma) {
            if(gsl_rng_uniform_pos(rng)<pt) {
                sigma[i_node] = 1;
                sigma[j_node] = 1;
            }
        }
    }

    for(i_node=0;i_node<nnode;i_node++) {
        i_sigma = sigma[i_node];
        if(i_sigma) {
            if(gsl_rng_uniform_pos(rng)<pr) {
                sigma[i_node] = 0;
            }
        }
    }
}

