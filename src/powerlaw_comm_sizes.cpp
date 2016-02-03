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
#include "combinatorics.h"
#include "print.h"

#ifdef __linux__
#include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;

void printUsage()
{
    mexPrintf("LFR Power law community size generator\n");
    mexPrintf("This program produces the community size assignment for given LFR parameters\n");
    mexPrintf("Input:\n");
    mexPrintf("powerlaw_comm_sizes(nmin,nmax,tau2,num_nodes)\n");
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
    if (nInputArgs<4)
    {
        printUsage();
        mexErrMsgTxt("Must specify 4 arguments!");
        return;
    }

    int nmin = *mxGetPr(inputArgs[0]);
    int nmax = *mxGetPr(inputArgs[1]);
    double tau2 = *mxGetPr(inputArgs[2]);
    int num_nodes = *mxGetPr(inputArgs[3]);

    srand4();
    deque<double> cumulative;
    deque<int> num_seq;

    int max_mem_num = 0;
    bool fixed_range=true;
    int max_degree_actual=0;
    int overlapping_nodes = 0;
    powerlaw(nmax,nmin,tau2,cumulative);

    if (num_seq.empty())
    {
        int _num_=0;
        if (!fixed_range && (max_degree_actual+1)>nmin)
        {
            _num_=max_degree_actual+1;			// this helps the assignment of the memberships (it assures that at least one module is big enough to host each node)
            num_seq.push_back(max_degree_actual+1);

        }

        while (true)
        {

            int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+nmin;

            if (nn+_num_<=num_nodes + overlapping_nodes * (max_mem_num-1) )
            {

                num_seq.push_back(nn);
                _num_+=nn;

            }
            else
                break;

        }

        num_seq[min_element(num_seq.begin(), num_seq.end()) - num_seq.begin()]+=num_nodes + overlapping_nodes * (max_mem_num-1) - _num_;
    }

    std::sort(num_seq.begin(),num_seq.end());

    vector<double>num_seq_dbl(num_seq.begin(),num_seq.end());
    outputArgs[0] = mxCreateDoubleMatrix(1,num_seq.size(),mxREAL);

    memcpy(mxGetPr(outputArgs[0]),num_seq_dbl.data(),num_seq_dbl.size()*sizeof(double));

    // Finish the function
    return;
}
