/* This file is part of PACO-PArtitioning Clustering Optimization a program
* to find network partitions using modular solvers and quality functions.
*
*  Copyright (C) 2015 Carlo Nicolini <carlo.nicolini@iit.it>
*
*  PACO is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 3 of the License, or (at your option) any later version.
*
*  Alternatively, you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  PACO is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License and a copy of the GNU General Public License along with
*  PACO. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <string.h>
#include <sstream>
#include <igraph.h>
#include <sys/time.h>

#include "benchm.h"
#include "set_parameters.h" // the LFR parameters
#include "../paco/igraph_utils.h" // to handle conversion from EigenMatrix to igraph object and then to mxArray

#ifdef __linux__
#include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;

void printUsage()
{
    mexPrintf("LFRW: Lancichinetti-Fortunato-Radicchi network generator\n");
    mexPrintf("This program produces weighted adjacency matrices for graph community detection benchmark\n");
    mexPrintf("[WeightedGraph, membership] = lfrw('argname',argvalue);\n");
    mexPrintf("Input:\n");
    mexPrintf("'N'\n");
    mexPrintf("\tDesidered number of nodes in the network\n");
    mexPrintf("'k'\n");
    mexPrintf("\tDesidered average degree of the network\n");
    mexPrintf("'maxk'\n");
    mexPrintf("\tDesidered max degree of the network\n");
    mexPrintf("'mut'\n");
    mexPrintf("\tDesidered community topological mixing coefficient (range [0,1])\n");
    mexPrintf("'muw'\n");
    mexPrintf("\tDesidered weights mixing coefficient (range [0,1])\n");
    mexPrintf("'t1'\n");
    mexPrintf("\tDesidered exponent for the degree distribution\n");
    mexPrintf("'t2'\n");
    mexPrintf("\tDesidered exponent for the community size distribution\n");
    mexPrintf("'beta'\n");
    mexPrintf("\tDesidered beta exponent\n");
    mexPrintf("'C'\n");
    mexPrintf("\tDesidered clustering coefficient\n");
}

enum error_type
{
    NO_ERROR = 0,
    ERROR_TOO_MANY_OUTPUT_ARGS = 1,
    ERROR_NOT_ENOUGH_ARGS = 2,
    ERROR_ARG_VALUE = 3,
    ERROR_ARG_TYPE = 4,
    ERROR_MATRIX = 5,
    ERROR_ARG_EMPTY=6,
    ERROR_ARG_UNKNOWN=7
};

static const char *error_strings[] =
{
    "",
    "Too many output arguments.",
    "Not enough input arguments.",
    "Non valid argument value.",
    "Non valid argument type.",
    "Non valid input adjacency matrix. Must be symmetric real dense-type square matrix.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};

/**
 * @brief parse_args
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 * @param pars
 * @param argposerr
 * @return
 */
error_type parse_args(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[], Parameters *pars, int *argposerr )
{
    if (nInputArgs < 1)
    {
        *argposerr = 0;
        return ERROR_NOT_ENOUGH_ARGS;
    }

    if (nOutputArgs>2)
    {
        *argposerr = 0;
        return ERROR_TOO_MANY_OUTPUT_ARGS;
    }

    // Iterate on function arguments
    int argcount=0;
    while (argcount<nInputArgs)
    {
        // Be sure that something exists after c-th argument
        if (argcount+1 >= nInputArgs)
        {
            *argposerr = argcount;
            return ERROR_ARG_EMPTY;
        }
        // Couple argument type - argument value
        const mxArray *partype = inputArgs[argcount];
        const mxArray *parval = inputArgs[argcount+1];
        char* cpartype;
        // To be a valid parameter specification it must be a couple ['char',real]
        if (mxIsChar(partype) && !mxIsChar(parval))
        {
            cpartype = mxArrayToString(partype);
            //mexPrintf("ARGUMENT: %s VALUE=%g\n", cpartype,*mxGetPr(parval));
            // Parse string value inputArgs[c]
            if ( strcasecmp(cpartype,"N")==0 )
            {
                pars->num_nodes = static_cast<int>((*mxGetPr(parval)));
                if (pars->num_nodes<0 )
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"k")==0 )
            {
                pars->average_k = static_cast<double>((*mxGetPr(parval)));
                if (pars->average_k <0 )
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"maxk")==0 )
            {
                pars->max_degree = static_cast<int>((*mxGetPr(parval)));
                if (pars->max_degree <0 )
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"mut")==0 )
            {
                pars->mixing_parameter_topological = static_cast<double>((*mxGetPr(parval)));
                if (pars->mixing_parameter_topological <0 || pars->mixing_parameter_topological > 1)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"muw")==0 )
            {
                pars->mixing_parameter_weights = static_cast<double>((*mxGetPr(parval)));
                if (pars->mixing_parameter_weights <0 || pars->mixing_parameter_weights > 1)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"t1")==0 )
            {
                pars->tau = static_cast<double>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"t2")==0 )
            {
                pars->tau2 = static_cast<double>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"minc")==0 )
            {
                pars->nmin = static_cast<int>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"maxc")==0 )
            {
                pars->nmax = static_cast<int>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"beta")==0 )
            {
                pars->beta = static_cast<double>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"on")==0 )
            {
                pars->overlapping_nodes = static_cast<int>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"om")==0 )
            {
                pars->overlap_membership = static_cast<int>((*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"C")==0 )
            {
                pars->clustering_coeff = static_cast<double>((*mxGetPr(parval)));
                argcount+=2;
            }
            else
            {
                *argposerr = argcount;
                return ERROR_ARG_UNKNOWN;
            }
        }
        else //else return the position of the argument and type of error
        {
            *argposerr = argcount;
            return ERROR_ARG_TYPE;
        }
        mxFree(cpartype); // free the converted argument
    }

    if (pars->arrange()==false)
        return ERROR_ARG_VALUE;

    return NO_ERROR;
}

/**
 * @brief mexFunction
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 */
void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    FILELog::ReportingLevel() = logDEBUG;
    // Set standard parameters argument
    Parameters p;

    // Check the arguments of the function
    int error_arg_pos=-1;
    error_type err = parse_args(nOutputArgs, outputArgs, nInputArgs, inputArgs, &p, &error_arg_pos);

    if (err!=NO_ERROR)
    {
        std::stringstream ss;
        ss << "Error at argument: " << error_arg_pos  << ": " << error_strings[err] ;
        if (err == ERROR_NOT_ENOUGH_ARGS)
            printUsage();
        mexErrMsgTxt(ss.str().c_str());
    }

    try
    {
        // Prepare output
        outputArgs[0] = mxCreateDoubleMatrix(p.num_nodes,p.num_nodes, mxREAL); // final weighted adjacency matrix
        outputArgs[1] = mxCreateDoubleMatrix(1,p.num_nodes, mxREAL); // membership of every vertex

        Eigen::MatrixXd W;
        vector<int> membership;
        // Call the main body of LFR weighted
        benchmark(p.excess, p.defect, p.num_nodes, p.average_k, p.max_degree, p.tau, p.tau2, p.mixing_parameter_topological,  p.mixing_parameter_weights,  p.beta, p.overlapping_nodes, p.overlap_membership, p.nmin, p.nmax, p.fixed_range, p.clustering_coeff, W, membership);

        // Copy the resulting matrix to output argument #0
        memcpy(mxGetPr(outputArgs[0]),W.data(),p.num_nodes*p.num_nodes*sizeof(double));

        // Copy the membership to output argument #1
        for (int i=0; i<p.num_nodes;++i)
            mxGetPr(outputArgs[1])[i]=membership.at(i);
        //memcpy(mxGetPr(outputArgs[1]),membership.data(),p.num_nodes*sizeof(int));
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
        mexErrMsgTxt(e.what());
    }

    // Finish the function
    return;
}
