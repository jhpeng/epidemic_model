#include <gsl/gsl_rng.h>

#ifndef graph_h
#define graph_h

int* read_edgelist(char* filename, int* nnode, int* nedge);

/** filename : input file name
***    nnode : output number of node
***    nedge : output number of edge
***   return : a list of edge
***/

int* kinship_graphs_generator(int* nnode, int* nedge, int N, double rho, gsl_rng* rng);

#endif
