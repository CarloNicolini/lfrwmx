#ifndef _COMBINATORICS_H_
#define _COMBINATORICS_H_

#include <set>
#include <deque>
#include <string>
#include <vector>
#include "random.h"

using std::set;
using std::deque;
using std::string;
using std::vector;


double log_factorial (int num);
double log_combination (int n, int k);
double binomial(int n, int x, double p);  //	returns the binomial distribution, n trials, x successes, p probability
//to draw a number:
int binomial_cumulative(int n, double p, deque<double> &cum);
// this function sets "cumulative" as the cumulative function of (1/x)^tau, with range= [min, n]
//to draw a number:
int powerlaw (int n, int min, double tau, deque<double> &cumulative);
int distribution_from_cumulative(const deque<double> &cum, deque<double> &distr); // cum is the cumulative, distr is set equal to the distribution
int cumulative_from_distribution (deque<double> &cum, const deque<double> &distr);  // cum is set equal to the cumulative, distr is the distribution
double poisson (int x, double mu);
int shuffle_and_set(int *due, int dim);  		// it sets due as a random sequence of integers from 0 to dim-1
int shuffle_s(deque<int> & sq);
double compute_r(int x, int k, int kout, int m);
int add_factors (deque<double> & num, deque<double> &den, int  n, int k);
double compute_hypergeometric(int i, int k, int kout, int m);
int random_from_set(set<int> & s);

template <typename Seq>
double average_func(Seq &sq)
{
    if (sq.empty())
        return 0;

    double av=0;
    typename Seq::iterator it = sq.begin();
    while(it != sq.end())
        av+=*(it++);

    av=av/sq.size();
    return av;
}

template <typename Seq>
double variance_func(Seq &sq)
{

    if (sq.empty())
        return 0;

    double av=0;
    double var=0;

    typename Seq::iterator it = sq.begin();
    while(it != sq.end())
    {

        av+=*(it);
        var+=(*(it))*(*(it));
        it++;

    }

    av=av/sq.size();
    var=var/sq.size();
    var-=av*av;

    if(var<1e-7)
        return 0;

    return var;

}

// this returns the average of the discrete probability function stored in Seq
template <typename Seq>
double average_pf(Seq &sq)
{

    double av=0;
    int h=0;

    typename Seq::iterator it = sq.begin();
    while(it != sq.end())
    {

        av+=*(it)*h;
        it++;
        h++;

    }

    return av;

}

template <typename Seq>
double variance_pf(Seq &sq)
{

    double av=0;
    double var=0;
    int h=0;

    typename Seq::iterator it = sq.begin();
    while(it != sq.end())
    {

        av+=*(it) * h;
        var+=(*(it)) * h * h ;
        it++;
        h++;
    }

    var-=av*av;

    if(var<1e-7)
        return 0;

    return var;

}


template <typename type_>
int shuffle_s(type_ *a, int b)
{
    int siz=b;
    if(siz==0)
        return -1;
    for (int i=0; i<b; i++)
    {
        int random_pos=irand(siz-1);
        type_ random_card_=a[random_pos];
        a[random_pos]=a[siz-1];
        a[siz-1]=random_card_;
        siz--;
    }
    return 0;
}

#endif
