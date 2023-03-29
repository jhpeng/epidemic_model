#include <gsl/gsl_rng.h>

#ifndef graph_h
#define graph_h

/**
 * @brief Reads the edge list from a file and returns an array of edges.
 *
 * This function reads an edge list from a file with the specified filename. It also
 * calculates the number of nodes and edges in the graph.
 *
 * @param filename The name of the file containing the edge list.
 * @param nnode Pointer to the integer representing the number of nodes.
 * @param nedge Pointer to the integer representing the number of edges.
 * @return Pointer to the integer array representing the edge list.
 */
int* read_edgelist(char* filename, int* nnode, int* nedge);

#endif
