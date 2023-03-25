#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "graphs.h"

int main() {
    unsigned long int seed = 398479375;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int nnode,nedge;
    char filename[128] = "/home/alan/Works/epidemic_model/ct_sis_model/graphs/test/test.edgelist";
    int* edges = read_edgelist(filename,&nnode,&nedge);

    for(int i=0;i<nedge;i++) {
        printf("(%d,%d)\n",edges[i*2+0],edges[i*2+1]);
    }
}
