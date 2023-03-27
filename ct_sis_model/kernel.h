#include <gsl/gsl_rng.h>

#ifndef kernel_h
#define kernel_h

/**
 * @brief Generates a random integer from a discrete distribution.
 *
 * @param[in] rng A pointer to a GSL random number generator.
 * @param[in] weight_dist An array containing the weights of the discrete distribution.
 * @param[in] total_weight The sum of the weights in the weight_dist array.
 * @param[in] len The length of the weight_dist array.
 *
 * @return A random integer selected according to the input discrete distribution.
 *
 * This function generates a random integer from a discrete distribution specified by the input
 * weight_dist array. The function uses the GSL random number generator to generate a uniform
 * random number, and then selects the corresponding integer based on the cumulative
 * distribution function (CDF) of the input discrete distribution.
 */
int random_integer_with_distribution(const gsl_rng *rng, double *weight_dist, double total_weight, int len);


/**
 * @brief Performs a uniform kernel scan on the input data.
 * 
 * @param[in] alpha The infection rate parameter.
 * @param[in] gamma The recovery rate parameter.
 * @param[out] dt The time step.
 * @param[in] nnode The number of nodes in the network.
 * @param[in] nedge The number of edges in the network.
 * @param[in] sigma An array representing the infection states of the nodes.
 * @param[in] edges An array representing the edges in the network.
 * @param[out] propose An array to store the nodes selected for infection or recovery.
 * @param[in] rng A GSL random number generator.
 * 
 * @return The number of proposed infection or recovery events.
 * 
 * This function performs a uniform kernel scan on the input data to identify
 * nodes to be infected or recovered based on the input parameters. The function
 * modifies the dt and propose arrays with the time step and the selected nodes,
 * respectively. The function returns the number of proposed infection or recovery
 * events.
 */
int kernel_uniform_scan(double alpha, double gamma, double* dt, int nnode, int nedge, const int* sigma, const int* edges, int* propose, gsl_rng* rng);

/**
 * @brief Updates the infection states of the nodes in the network.
 *
 * @param[in] np The number of proposed infection or recovery events.
 * @param[in] propose An array containing the indices of the nodes selected for infection or recovery.
 * @param[in,out] sigma An array representing the infection states of the nodes; this array is modified in-place to update the infection states.
 *
 * This function updates the infection states of the nodes in the network based on the proposed infection or recovery events.
 * The infection state of each node in the `propose` array is toggled (infected nodes become recovered and recovered nodes become infected).
 */
void kernel_update(int np, const int* propose, int* sigma);


/**
 * @brief Frees the memory allocated for kernel_uniform_recover_list and kernel_uniform_infected_list.
 *
 * This function checks if the global pointers kernel_uniform_recover_list and kernel_uniform_infected_list
 * are not NULL. If they are not NULL, it frees the memory allocated for these lists and sets their pointers
 * to NULL.
 */
void kernel_free();


#endif
