/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *    This program is free software; you can redistribute it and/or modify       *
 *  it under the terms of the GNU General Public License as published by         *
 *  the Free Software Foundation; either version 2 of the License, or            *
 *  (at your option) any later version.                                          *
 *                                                                               *
 *  This program is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 *  GNU General Public License for more details.                                 *
 *                                                                               *
 *  You should have received a copy of the GNU General Public License            *
 *  along with this program; if not, write to the Free Software                  *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *  Created by Andrea Lancichinetti on 7/01/09 (email: arg.lanci@gmail.com)      *
 *    Modified on 28/05/09                                                       *
 *    Collaborators: Santo Fortunato                                             *
 *  Location: ISI foundation, Turin, Italy                                       *
 *    Project: Benchmarking community detection programs                         *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

#include <Eigen/Core>
#include "binary_benchm.h"
#include "histograms.h"
#include "print.h"

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
                  deque<int> & num_seq, deque<map <int, double > > & neigh_weigh, double beta, double mu, double mu0, Eigen::MatrixXd &W)
{
    int edges=0;
    int num_nodes=member_list.size();
    deque<double> double_mixing;
    for (unsigned int i=0; i<E.size(); i++)
    {
        double one_minus_mu = double(internal_kin(E, member_list, i))/E[i].size();
        double_mixing.push_back(1.- one_minus_mu);
        edges+=E[i].size();
    }

    double density=0;
    double sparsity=0;
    for (unsigned int i=0; i<member_matrix.size(); i++)
    {
        double media_int=0;
        double media_est=0;
        for (unsigned int j=0; j<member_matrix[i].size(); j++)
        {
            double kinj = double(internal_kin_only_one(E[member_matrix[i][j]], member_matrix[i]));
            media_int+= kinj;
            media_est+=E[member_matrix[i][j]].size() - double(internal_kin_only_one(E[member_matrix[i][j]], member_matrix[i]));
        }
        double pair_num=(member_matrix[i].size()*(member_matrix[i].size()-1));
        double pair_num_e=((num_nodes-member_matrix[i].size())*(member_matrix[i].size()));

        if(pair_num!=0)
            density+=media_int/pair_num;
        if(pair_num_e!=0)
            sparsity+=media_est/pair_num_e;
        // CARLO Inserito ma non corretto
        // FILE_LOG(logINFO) << "Cluster average intra "<<  media_int;
        // FILE_LOG(logINFO) << "Cluster average extra "<<  media_est;
    }
    density=density/member_matrix.size();
    sparsity=sparsity/member_matrix.size();

    // CARLO - Inserted output matrix in adjacency matrix format
    W.setZero(num_nodes,num_nodes);

    // Fill the adjacency matrix
    unsigned int u = 0;
    for (deque< set<int> >::iterator it = E.begin(); it!=E.end();++it)
    {
        for (set<int>::iterator it2 = it->begin(); it2!=it->end();++it2)
        {
            int v = *it2;
            W.coeffRef(u,v) = neigh_weigh[u][v];
        }
        ++u;
    }

    FILE_LOG(logINFO) << "Network of "<<num_nodes<<" vertices and "<<edges/2<<" edges"<<";\t average degree = "<< double(edges)/num_nodes;
    FILE_LOG(logINFO) << "Average mixing parameter (topology): "<< average_func(double_mixing)<<" +/- "<<sqrt(variance_func(double_mixing));
    FILE_LOG(logINFO) << "p_in: " << density << "\tp_out: " << sparsity;

    deque<int> degree_seq;
    for (unsigned int i=0; i<E.size(); i++)
        degree_seq.push_back(E[i].size());

    /*
    ofstream statout("statistics.dat");
    statout<<"degree distribution (probability density function of the degree in logarithmic bins) " << endl;
    log_histogram(degree_seq, statout, 10);
    statout<<"\ndegree distribution (degree-occurrences) " << endl;
    int_histogram(degree_seq, statout);
    statout<<endl<<"--------------------------------------" << endl;

    statout<<"community distribution (size-occurrences)" << endl;
    int_histogram(num_seq, statout);
    statout<<endl<<"--------------------------------------" << endl;

    statout<<"mixing parameter (topology)" << endl;
    not_norm_histogram(double_mixing, statout, 20, 0, 0);
    statout<<endl<<"--------------------------------------" << endl;
    */

    deque<double> inwij;
    deque<double> outwij;

    //double csi=(1. - mu) / (1. - mu0);
    //double csi2=mu /mu0;
    double tstrength=0;
    deque<double> one_minus_mu2;

    for(unsigned int i=0; i<neigh_weigh.size(); i++)
    {

        double internal_strength_i=0;
        double strength_i=0;
        for(map<int, double>::iterator itm = neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++)
        {
            if(they_are_mate(i, itm->first, member_list))
            {
                inwij.push_back(itm->second);
                //inkij.push_back(csi * pow(E[i].size(), beta-1));
                internal_strength_i+=itm->second;
            }
            else
            {
                outwij.push_back(itm->second);
                //outkij.push_back(csi2 * pow(E[i].size(), beta-1));
            }
            tstrength+=itm->second;
            strength_i+=itm->second;
        }
        one_minus_mu2.push_back(1 - internal_strength_i/strength_i);
    }

    FILE_LOG(logINFO)<<"average mixing parameter (weights): "<<average_func(one_minus_mu2)<<" +/- "<<sqrt(variance_func(one_minus_mu2));
    //statout<<"mixing parameter (weights)" << endl;
    //not_norm_histogram(one_minus_mu2, statout, 20, 0, 0);
    //statout<<endl<<"--------------------------------------" << endl;

    FILE_LOG(logINFO)<<"average weight of an internal link "<<average_func(inwij)<<" +/- "<<sqrt(variance_func(inwij))<<endl;
    FILE_LOG(logINFO)<<"average weight of an external link "<<average_func(outwij)<<" +/- "<<sqrt(variance_func(outwij))<<endl;

    //statout<<"internal weights (weight-occurrences)" << endl;
    //not_norm_histogram(inwij, statout, 20, 0, 0);
    //statout<<endl<<"--------------------------------------" << endl;

    //statout<<"external weights (weight-occurrences)" << endl;
    //not_norm_histogram(outwij, statout, 20, 0, 0);

    return 0;
}

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
                  deque<deque<double> > & wished, deque<deque<double> > & factual, const double tot_var, double *strs)
{
    double d1t=0;
    double d2t=0;
    double d3t=0;
    double var_check=0;

    for (unsigned int k=0; k<member_list.size(); k++)
    {
        double in_s=0;
        double out_s=0;

        for(map<int, double>::iterator itm=neigh_weigh[k].begin(); itm!=neigh_weigh[k].end(); itm++)
        {
            if(itm->second<0)
                cherr(333);
            if(fabs(itm->second - neigh_weigh[itm->first][k]) > 1e-7)
                cherr(111);

            if (they_are_mate(k, itm->first, member_list))
                in_s+=itm->second;
            else
                out_s+=itm->second;
        }

        if (fabs(in_s - factual[k][0]) > 1e-7)
            cherr(in_s - factual[k][0]);
        //cout<<"ok1" << endl;

        if (fabs(out_s - factual[k][1]) > 1e-7)
            cherr(out_s - factual[k][1]);
        //cout<<"ok2" << endl;

        if (fabs(in_s + out_s + factual[k][2] - strs[k]) > 1e-7)
            cherr(in_s + out_s + factual[k][2] - strs[k]);
        //cout<<"ok3" << endl;

        double d1=(in_s - wished[k][0]);
        double d2=(out_s - wished[k][1]);
        double d3=(strs[k] - in_s - out_s);
        var_check+= d1*d1 + d2*d2 + d3*d3;

        d1t+=d1*d1;
        d2t+=d2*d2;
        d3t+=d3*d3;
    }

    FILE_LOG(logINFO)<<"tot_var "<<tot_var<<"\td1t "<<d1t<<"\td2t "<<d2t<<"\td3t "<<d3t;
    if (fabs(var_check - tot_var) > 1e-5)
        FILE_LOG(logWARNING)<<"found this difference in check "<<fabs(var_check - tot_var);
    else
        FILE_LOG(logINFO)<<"Ok: check passed";

    return 0;

}

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
              deque<deque<double> > & wished, deque<deque<double> > & factual, int i, double & tot_var, double *strs, const deque<int> & internal_kin_top)
{
    {
        // in this case I rewire the idle strength

        double change=factual[i][2]/neigh_weigh[i].size();

        double oldpartvar=0;
        for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) if(itm->second + change > 0)
            for (int bw=0; bw<3; bw++)
                oldpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);

        for (int bw=0; bw<3; bw++)
            oldpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);

        double newpartvar=0;

        for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++) if(itm->second + change > 0)
        {

            if (they_are_mate(i, itm->first, member_list))
            {
                factual[itm->first][0]+=change;
                factual[itm->first][2]-=change;

                factual[i][0]+=change;
                factual[i][2]-=change;
            }
            else
            {
                factual[itm->first][1]+=change;
                factual[itm->first][2]-=change;

                factual[i][1]+=change;
                factual[i][2]-=change;
            }
            for (int bw=0; bw<3; bw++)
                newpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
            itm->second+= change;
            neigh_weigh[itm->first][i]+=change;
        }

        for (int bw=0; bw<3; bw++)
            newpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);
        tot_var+= newpartvar - oldpartvar;
    }

    int internal_neigh=internal_kin_top[i];

    if(internal_neigh!=0)          // in this case I rewire the difference strength
    {
        double changenn=(factual[i][0] - wished[i][0]);
        double oldpartvar=0;
        for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++)
        {
            if(they_are_mate(i, itm->first, member_list))
            {
                double change = changenn/internal_neigh;
                if(itm->second - change > 0)
                    for (int bw=0; bw<3; bw++)
                        oldpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
            }
            else
            {
                double change = changenn/(neigh_weigh[i].size() - internal_neigh);

                if(itm->second + change > 0)
                    for (int bw=0; bw<3; bw++)
                        oldpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);
            }
        }

        for (int bw=0; bw<3; bw++)
            oldpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);

        double newpartvar=0;

        for(map<int, double>::iterator itm=neigh_weigh[i].begin(); itm!=neigh_weigh[i].end(); itm++)
        {
            if (they_are_mate(i, itm->first, member_list))
            {
                double change = changenn/internal_neigh;

                if(itm->second - change > 0)
                {
                    factual[itm->first][0]-=change;
                    factual[itm->first][2]+=change;

                    factual[i][0]-=change;
                    factual[i][2]+=change;

                    for (int bw=0; bw<3; bw++)
                        newpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);

                    itm->second-= change;
                    neigh_weigh[itm->first][i]-=change;
                }
            }

            else
            {
                double change = changenn/(neigh_weigh[i].size() - internal_neigh);
                if(itm->second + change > 0)
                {
                    factual[itm->first][1]+=change;
                    factual[itm->first][2]-=change;

                    factual[i][1]+=change;
                    factual[i][2]-=change;

                    for (int bw=0; bw<3; bw++)
                        newpartvar+= (factual[itm->first][bw] - wished[itm->first][bw]) * (factual[itm->first][bw] - wished[itm->first][bw]);

                    itm->second+= change;
                    neigh_weigh[itm->first][i]+=change;
                }
            }
        }

        for (int bw=0; bw<3; bw++)
            newpartvar+= (factual[i][bw] - wished[i][bw]) * (factual[i][bw] - wished[i][bw]);

        tot_var+=newpartvar - oldpartvar;
    }

    //check_weights(neigh_weigh, member_list, wished, factual, tot_var, strs);
    return 0;

}

/**
 * @brief weights
 * @param en
 * @param member_list
 * @param beta
 * @param mu
 * @param neigh_weigh
 * @return
 */
int weights(deque<set<int> > & en, const deque<deque<int> > & member_list, const double beta, const double mu, deque<map <int, double > > & neigh_weigh)
{
    double tstrength=0;
    for(unsigned int i=0; i<en.size(); i++)
        tstrength+=pow(en[i].size(), beta);

    double strs[en.size()]; // strength of the nodes
    // build a matrix like this: deque < map <int, double > > each row corresponds to link - weights
    for(unsigned int i=0; i<en.size(); i++)
    {

        map<int, double> new_map;
        neigh_weigh.push_back(new_map);

        for (set<int>::iterator its=en[i].begin(); its!=en[i].end(); its++)
            neigh_weigh[i].insert(make_pair(*its, 0.));

        strs[i]=pow(double(en[i].size()), beta);
    }

    deque<double> s_in_out_id_row(3);
    s_in_out_id_row[0]=0;
    s_in_out_id_row[1]=0;
    s_in_out_id_row[2]=0;

    // 3 numbers for each node: internal, idle and extra strength. the sum of the three is strs[i].
    // wished is the theoretical, factual the factual one.
    deque<deque<double> > wished;
    deque<deque<double> > factual;

    for (unsigned int i=0; i<en.size(); i++)
    {

        wished.push_back(s_in_out_id_row);
        factual.push_back(s_in_out_id_row);
    }
    double tot_var=0;
    for (unsigned int i=0; i<en.size(); i++)
    {

        wished[i][0]=(1. -mu)*strs[i];
        wished[i][1]=mu*strs[i];

        factual[i][2]=strs[i];

        tot_var+= wished[i][0] * wished[i][0] + wished[i][1] * wished[i][1] + strs[i] * strs[i];
    }
    deque<int> internal_kin_top;
    for(unsigned int i=0; i<en.size(); i++)
        internal_kin_top.push_back(internal_kin(en, member_list, i));

    double precision = 1e-9;
    double precision2 = 1e-2;
    double not_better_than = pow(tstrength, 2) * precision;

    int step=0;
    while (true)
    {

        //time_t t0=time(NULL);

        double pre_var=tot_var;
        for (unsigned int i=0; i<en.size(); i++)
            propagate(neigh_weigh, member_list, wished, factual, i, tot_var, strs, internal_kin_top);
        double relative_improvement=double(pre_var - tot_var)/pre_var;

        if (tot_var<not_better_than)
            break;

        if (relative_improvement < precision2)
            break;

        //time_t t1= time(NULL);
        //int deltat= t1 - t0;
        step++;
    }
    return 0;

}

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
 * @return
 */
int benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2,
              double  mixing_parameter, double  mixing_parameter2, double  beta, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, bool  fixed_range, double ca, Eigen::MatrixXd &W, vector<int> &membership)
{

    double dmin=solve_dmin(max_degree, average_k, -tau);
    if (dmin==-1)
        return -1;

    int min_degree=int(dmin);

    double media1=integer_average(max_degree, min_degree, tau);
    double media2=integer_average(max_degree, min_degree+1, tau);

    if (fabs(media1-average_k)>fabs(media2-average_k))
        min_degree++;

    // Range for the community sizes
    if (!fixed_range)
    {
        nmax=max_degree;
        nmin=max(int(min_degree), 3);
        FILE_LOG(logINFO) << "community size range automatically set equal to ["<<nmin<<" , "<<nmax<<"]" ;
    }

    // Degree sequence of the nodes
    deque <int> degree_seq ;
    deque <double> cumulative;
    powerlaw(max_degree, min_degree, tau, cumulative);

    for (int i=0; i<num_nodes; i++)
    {
        int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+min_degree;
        degree_seq.push_back(nn);
    }

    sort(degree_seq.begin(), degree_seq.end());

    if(deque_int_sum(degree_seq) % 2!=0)
        degree_seq[max_element(degree_seq.begin(), degree_seq.end()) - degree_seq.begin()]--;

    deque<deque<int> >  member_matrix;
    deque<int> num_seq;
    deque<int> internal_degree_seq;

    if(internal_degree_and_membership(mixing_parameter, overlapping_nodes, overlap_membership, num_nodes, member_matrix, excess, defect, degree_seq, num_seq, internal_degree_seq, fixed_range, nmin, nmax, tau2)==-1)
        return -1;

    deque<set<int> > E;                    // E is the adjacency matrix written in form of list of edges
    deque<deque<int> > member_list;        // row i cointains the memberships of node i
    deque<deque<int> > link_list;        // row i cointains degree of the node i respect to member_list[i][j]; there is one more number that is the external degree

    FILE_LOG(logINFO) << "Building communities... " ;
    if(build_subgraphs(E, member_matrix, member_list, link_list, internal_degree_seq, degree_seq, excess, defect)==-1)
        return -1;

    FILE_LOG(logINFO) << "Connecting communities... " ;
    connect_all_the_parts(E, member_list, link_list);

    if(erase_links(E, member_list, excess, defect, mixing_parameter)==-1)
        return -1;

    if(ca!=unlikely_value)
    {
        FILE_LOG(logINFO) << "Trying to approach an average clustering coefficient ... " << ca;
        cclu(E, member_list, member_matrix, ca);
    }

    deque<map <int, double > > neigh_weigh;
    weights(E, member_list, beta, mixing_parameter2, neigh_weigh);
    print_network(E, member_list, member_matrix, num_seq, neigh_weigh, beta, mixing_parameter2, mixing_parameter, W);

    // Copy the membership
    membership.resize(member_list.size());
    for (unsigned int i=0; i<member_list.size();++i)
    {
        //cout << "->" << member_list[i][0] << endl;
        membership[i]=member_list[i][0];
    }

    return 0;
}

int benchmark_py(int excess,
                 int defect,
                 int num_nodes,
                 double average_k,
                 int max_degree,
                 double tau,
                 double tau2,
                 double mixing_parameter,
                 double mixing_parameter2,
                 double beta,
                 int  overlapping_nodes,
                 int overlap_membership,
                 int nmin,
                 int nmax,
                 int fixed_range,
                 double ca,
                 double *W,
                 int *membership)
{
    Eigen::MatrixXd WM;
    vector<int> vmembership;
    int val = benchmark(excess, defect, num_nodes, average_k, max_degree, tau, tau2,
                        mixing_parameter, mixing_parameter2, beta, overlapping_nodes,
                        overlap_membership, nmin, nmax, fixed_range, ca, WM, vmembership);
    membership = vmembership.data();
    W = WM.data();
    //cout << WM << endl;
//    for (int i=0; i<vmembership.size();++i)
//        cout << membership[i] << endl;
    return val;
}

/**
 * @brief erase_file_if_exists
 * @param s
 */
void erase_file_if_exists(string s)
{
    char b[100];
    cast_string_to_char(s, b);
    ifstream in1(b);
    if(in1.is_open())
    {
        char rmb[120];
        sprintf(rmb, "rm %s", b);
        //int erase= system(rmb);
    }
}
