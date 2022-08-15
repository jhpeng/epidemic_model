#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
    }

    return sigma;
}

static void print_state(int* sigma, int nnode) {
    for(int i=0;i<nnode;i++) {
        printf("%d ",sigma[i]);
    }
    printf("\n");
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
    double alpha=1.0;
    double gamma=1.0;
    double dt=0.001;
    double T=100.0;
    double p=0.05;
    unsigned long int seed=8239127933;
    int nstep=(int)(T/dt);
    int nshow=(int)(0.1/dt);

    gsl_rng* rng=gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int nnode,nedge;
    char filename[128] = "../../networks/jupyters/test.edgelist";
    int* edges = read_edgelist(filename,&nnode,&nedge);

    int* sigma = initial_state(nnode,p,0,0,rng);

    print_state(sigma,nnode);
    double sigma_ave = order_parameter(sigma,nnode);
    for(int i=0;i<nstep;i++) {
        kernel_linear(alpha,gamma,dt,nnode,nedge,sigma,edges,rng);
        if((i+1)%nshow==0) {
            //sigma_ave = order_parameter(sigma,nnode);

            //printf("T=%.f sigma_ave=%.8lf \n",dt*(i+1),sigma_ave);
            print_state(sigma,nnode);
        }
    }
    


    free(edges);
}
