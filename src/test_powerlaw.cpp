#include <iostream>
#include <deque>
#include <algorithm>
#include "combinatorics.h"
#include "print.h"

using namespace std;

int main(int argc, char *argv[])
{
    srand4();
    int nmin = atoi(argv[1]);
    int nmax = atoi(argv[2]);
    double tau2 = atof(argv[3]);
    deque<double> cumulative;
    deque<int> num_seq;
    int num_nodes = atoi(argv[4]);
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
    cout << "|C|=" << num_seq.size() << endl;
    cout<<"community sizes: ";
    int total_nodes = 0;
    for (int i=0; i<num_seq.size(); i++)
    {
        cout<<num_seq[i]<<" ";
        total_nodes += num_seq[i];
    }
    cout<<endl;
    cout << "Total= " << total_nodes << " Asked= " << num_nodes << endl;
    ofstream outhist; outhist.open("bars.txt");
    
    for (int i=0; i<num_seq.size(); i++)
    {
        outhist << num_seq[i] << "\n";
    }

    return 0;
}
