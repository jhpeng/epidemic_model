#include <gsl/gsl_rng.h>

#ifndef time_evolution_h
#define time_evolution_h

/**
 * @brief Initializes the kernel for simulating an epidemic process.
 *
 * This function initializes the kernel and allocates memory for kernel_tiks. It also resets all
 * relevant counters and statistics.
 *
 * @param nnode Number of nodes in the graph.
 */
void kernel_init(int nnode);

/**
 * @brief Simulates one step of an epidemic process using a linear kernel.
 *
 * This function simulates one step of an epidemic process on a given graph with a linear kernel.
 * It updates the infection and recovery status of the nodes according to the provided parameters.
 *
 * @param alpha Infection rate parameter.
 * @param gamma Recovery rate parameter.
 * @param dt Time step for the simulation.
 * @param nnode Number of nodes in the graph.
 * @param nedge Number of edges in the graph.
 * @param sigma Integer array representing the infection status of each node.
 * @param edges Integer array representing the edge list of the graph.
 * @param rng Pointer to a gsl_rng random number generator.
 */
void kernel_linear(double alpha, double gamma, double dt, int nnode, int nedge, int* sigma, int* edges, gsl_rng* rng);

/**
 * @brief Closes the kernel and updates the total infected time.
 *
 * This function updates the total infected time for each infected node and should be called
 * after the simulation is completed.
 *
 * @param nnode Number of nodes in the graph.
 * @param sigma Integer array representing the infection status of each node.
 */
void kernel_close(int nnode, int* sigma);

/**
 * @brief Returns the total infected time.
 *
 * @return The total infected time as a double.
 */
double total_infected_time_value();

/**
 * @brief Returns the total number of infections.
 *
 * @return The total number of infections as a double.
 */
double total_ninfection_value();

/**
 * @brief Returns the total number of recoveries.
 *
 * @return The total number of recoveries as a double.
 */
double total_nrecover_value();

#endif
