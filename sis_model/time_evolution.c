#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>


static int kernel_count=0;
static unsigned long int total_infected_time=0;
static unsigned long int total_ninfection=0;
static unsigned long int total_nrecover=0;
static int* kernel_tiks=NULL;
void kernel_init(int nnode) {
    kernel_count=0;
    total_infected_time=0;
    total_ninfection=0;
    total_nrecover=0;
    if(kernel_tiks==NULL) 
        kernel_tiks = (int*)malloc(sizeof(int)*nnode);

    for(int i=0;i<nnode;i++)
        kernel_tiks[i]=0;
}

void kernel_linear(double alpha, double gamma, double dt, int nnode, int nedge, int* sigma, int* edges, gsl_rng* rng) {
    double pt = alpha*dt;
    double pr = gamma*dt;
    int i_node,j_node,i_sigma,j_sigma;

    int ninfection = 0;
    for(int i_edge=0;i_edge<nedge;i_edge++) {
        i_node = edges[i_edge*2+0];
        j_node = edges[i_edge*2+1];

        i_sigma = sigma[i_node];
        j_sigma = sigma[j_node];
        if(i_sigma^j_sigma) {
            if(gsl_rng_uniform_pos(rng)<pt) {
                if(sigma[i_node]) {
                    kernel_tiks[j_node] = kernel_count;
                } else {
                    kernel_tiks[i_node] = kernel_count;
                }
                sigma[i_node] = 1;
                sigma[j_node] = 1;
                ninfection++;
            }
        }
    }
    total_ninfection += ninfection;

    int nrecover=0;
    for(i_node=0;i_node<nnode;i_node++) {
        i_sigma = sigma[i_node];
        if(i_sigma) {
            if(gsl_rng_uniform_pos(rng)<pr) {
                sigma[i_node] = 0;
                nrecover++;
                total_infected_time += kernel_count-kernel_tiks[i_node];
            }
        }
    }
    total_nrecover += nrecover;
    kernel_count++;
}

void kernel_close(int nnode, int* sigma) {
    for(int i=0;i<nnode;i++) {
        if(sigma[i]) 
            total_infected_time += kernel_count-kernel_tiks[i];
    }
}

double total_infected_time_value() {
    return (double)total_infected_time;
}

double total_ninfection_value() {
    return (double)total_ninfection;
}

double total_nrecover_value() {
    return (double)total_nrecover;
}
