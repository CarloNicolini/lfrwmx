#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "standard_include.h"
#include "cast.h"

class Parameters
{
public:
    Parameters();
    ~Parameters() {}

    int num_nodes;
    double average_k;
    int max_degree;
    double tau_degree;
    double tau_commsize;
    double mixing_parameter_topological;
    double mixing_parameter_weights;
    double beta; //exponent
    int overlapping_nodes;
    int overlap_membership;
    int nmin;
    int nmax;
    bool fixed_range;
    bool excess;
    bool defect;
    bool randomf;
    double clustering_coeff;

    bool set(string &, string &);
    void set_random();
    int arrange();
    deque<string> command_flags;
    void print_parameters();
    int verbosity;
};

void print_usage();
bool set_from_file(string & file_name, Parameters & par1);
bool set_parameters(int argc, char * argv[], Parameters & par1);

#endif
