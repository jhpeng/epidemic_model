#include <gsl/gsl_rng.h>

#ifndef graph_h
#define graph_h

/**
 * @brief Reads an edge list from a file and returns it as an integer array.
 *
 * This function reads a file containing an edge list, determines the number of nodes and edges,
 * and returns the edge list as a dynamically allocated integer array. The format of the file is
 * assumed to be an adjacency list, where each line contains a pair of connected node indices and
 * a string (which is ignored). The node indices are assumed to be zero-based.
 *
 * @param filename Pointer to a string containing the name of the file to read.
 * @param[out] nnode Pointer to an integer where the number of nodes in the graph will be stored.
 * @param[out] nedge Pointer to an integer where the number of edges in the graph will be stored.
 * @return Pointer to a dynamically allocated integer array containing the edge list, or NULL if the file could not be opened.
 */
int* read_edgelist(char* filename, int* nnode, int* nedge);

/**
 * @brief Generates a kinship graph with a given number of nodes and edges using a random number generator.
 *
 * This function generates a kinship graph based on the provided parameters. The graph is generated
 * by randomly connecting nodes according to a specified probability (rho). The function also
 * handles the allocation of memory for the generated edges.
 *
 * @param[out] nnode Pointer to an integer where the number of nodes in the graph will be stored.
 * @param[out] nedge Pointer to an integer where the number of edges in the graph will be stored.
 * @param N Number of groups in the graph.
 * @param rho Probability parameter used for determining node connections.
 * @param rng Pointer to a gsl_rng random number generator.
 * @return Pointer to a dynamically allocated integer array containing the edge list.
 */
int* kinship_graphs_generator(int* nnode, int* nedge, int N, double rho, gsl_rng* rng);

#endif
