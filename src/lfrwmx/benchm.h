#ifndef _BENCHM_H_
#define _BENCHM_H_

#include "binary_benchm.h"

/**
 * @brief print_network
 * @param E
 * @param member_list
 * @param member_matrix
 * @param num_seq
 * @param neigh_weigh
 * @param beta
 * @param mu
 * @param mu0
 * @return
 */
int print_network(deque<set<int> > & E, const deque<deque<int> > & member_list, const deque<deque<int> > & member_matrix,
                  deque<int> & num_seq, deque<map <int, double > > & neigh_weigh, double beta, double mu, double mu0);

/**
 * @brief check_weights
 * @param neigh_weigh
 * @param member_list
 * @param wished
 * @param factual
 * @param tot_var
 * @param strs
 * @return
 */
int check_weights(deque<map <int, double > > & neigh_weigh, const deque<deque<int> > & member_list,
                  deque<deque<double> > & wished, deque<deque<double> > & factual, const double tot_var, double *strs);

/**
 * @brief propagate
 * @param neigh_weigh
 * @param member_list
 * @param wished
 * @param factual
 * @param i
 * @param tot_var
 * @param strs
 * @param internal_kin_top
 * @return
 */
int propagate(deque<map <int, double > > & neigh_weigh, const deque<deque<int> > & member_list,
              deque<deque<double> > & wished, deque<deque<double> > & factual, int i, double & tot_var, double *strs, const deque<int> & internal_kin_top);

/**
 * @brief weights
 * @param en
 * @param member_list
 * @param beta
 * @param mu
 * @param neigh_weigh
 * @return
 */
int weights(deque<set<int> > & en, const deque<deque<int> > & member_list, const double beta, const double mu, deque<map <int, double > > & neigh_weigh);

/**
 * @brief benchmark
 * @param excess
 * @param defect
 * @param num_nodes
 * @param average_k
 * @param max_degree
 * @param tau
 * @param tau2
 * @param mixing_parameter
 * @param mixing_parameter2
 * @param beta
 * @param overlapping_nodes
 * @param overlap_membership
 * @param nmin
 * @param nmax
 * @param fixed_range
 * @param ca
 * @param W
 * @return
 */
int benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2, double  mixing_parameter, double  mixing_parameter2, double  beta, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, bool  fixed_range, double ca, Eigen::MatrixXd &W, vector<int> &membership);

int benchmark_py(int excess, int defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2, double  mixing_parameter, double  mixing_parameter2, double  beta, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, int fixed_range, double ca, double *W, int *membership);

/**
 * @brief erase_file_if_exists
 * @param s
 */
void erase_file_if_exists(string s);

#endif
