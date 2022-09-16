#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "graph.h"
#include "time_evolution.h"

static int* initial_state(int nnode, double p, int n, int type, gsl_rng* rng) {
    int* sigma = (int*)malloc(sizeof(int)*nnode);

    if(type==0) {
        for(int i=0;i<nnode;i++) {
            sigma[i]=0;
            if(gsl_rng_uniform_pos(rng)<p) sigma[i]=1;
        }
    } else if(type==1) {
        for(int i=0;i<nnode;i++) {
            sigma[i]=0;
            if(i<n) sigma[i]=1;
        }
    }

    return sigma;
}

static double order_parameter(int* sigma, int nnode) {
    double I=0;
    for(int i=0;i<nnode;i++) {
        I+=sigma[i];
    }
    I = I/nnode;

    return I;
}

int main() {
    double gamma=1.0;
    double alpha=1.0;
    double dt=0.01;
    unsigned long int seed = 384932;

    int nthermal=100000;
    int nsample=1000;
    int nstep=(int)(1.0/dt);

    gsl_rng* rng=gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    
    int nnode,nedge;
    char filename[128] = "/home/alan/Works/path_sampling/networks/jupyters/N_100.edgelist";
    int* edges = read_edgelist(filename,&nnode,&nedge);

    int* sigma = initial_state(nnode,0.5,0,0,rng);

    double infected_ratio=0;
    int block_size=1000;
    for(int k=0;k<nsample;k++) {
        infected_ratio=0;
        for(int j=0;j<block_size;j++) {
            for(int i=0;i<nstep;i++) {
                kernel_linear(alpha,gamma,dt,nnode,nedge,sigma,edges,rng);
            }

            infected_ratio += order_parameter(sigma,nnode);
        }
        infected_ratio = infected_ratio/block_size;

        printf("%.12f\n",infected_ratio);
    }
}
