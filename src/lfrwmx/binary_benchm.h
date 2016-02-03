#ifndef _BINARY_BENCHM_H_
#define _BINARY_BENCHM_H_

#include "standard_include.h"
#include "cc.h"
#include "random.h"
#include "combinatorics.h"
#include "set_parameters.h"

bool they_are_mate(int a, int b, const deque<deque<int> > & member_list);
// it computes the sum of a deque<int>
int deque_int_sum(const deque<int> & a);
// it computes the integral of a power law
double integral (double a, double b);

// it returns the average degree of a power law
double average_degree(const double &dmax, const double &dmin, const double &gamma);

//bisection method to find the inferior limit, in order to have the expected average degree
double solve_dmin(const double& dmax, const double &dmed, const double &gamma);

// it computes the correct (i.e. discrete) average of a power law
double integer_average (int n, int min, double tau);
// this function changes the community sizes merging the smallest communities
int change_community_size(deque<int> &seq);
int build_bipartite_network(deque<deque<int> >  & member_matrix, const deque<int> & member_numbers, const deque<int> &num_seq);
int internal_degree_and_membership (double mixing_parameter, int overlapping_nodes, int max_mem_num, int num_nodes, deque<deque<int> >  & member_matrix,
                                    bool excess, bool defect,  deque<int> & degree_seq, deque<int> &num_seq, deque<int> &internal_degree_seq, bool fixed_range, int nmin, int nmax, double tau2);

int build_subgraph(deque<set<int> > & E, const deque<int> & nodes, const deque<int> & degrees);
int build_subgraphs(deque<set<int> > & E, const deque<deque<int> > & member_matrix, deque<deque<int> > & member_list,
                    deque<deque<int> > & link_list, const deque<int> & internal_degree_seq, const deque<int> & degree_seq, const bool excess, const bool defect);

int connect_all_the_parts(deque<set<int> > & E, const deque<deque<int> > & member_list, const deque<deque<int> > & link_list);
int internal_kin(deque<set<int> > & E, const deque<deque<int> > & member_list, int i);
int internal_kin_only_one(set<int> & E, const deque<int> & member_matrix_j);
int erase_links(deque<set<int> > & E, const deque<deque<int> > & member_list, const bool excess, const bool defect, const double mixing_parameter);

#endif
