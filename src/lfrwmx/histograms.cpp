#include "histograms.h"
#include "print.h"

void int_histogram (vector <int> &c, ostream & out)
{

    map<int, double> hist;

    double freq=1/double(c.size());

    for (int i=0; i<c.size(); i++)
    {

        map<int, double>::iterator itf=hist.find(c[i]);
        if (itf==hist.end())
            hist.insert(make_pair(c[i], 1.));
        else
            itf->second++;

    }

    for (map<int, double>::iterator it=hist.begin(); it!=hist.end(); it++)
        it->second=it->second*freq;

    prints(hist, out);
}

void int_histogram (deque <int> &c, ostream & out)
{

    map<int, double> hist;

    double freq=1/double(c.size());

    for (int i=0; i<c.size(); i++)
    {

        map<int, double>::iterator itf=hist.find(c[i]);
        if (itf==hist.end())
            hist.insert(make_pair(c[i], 1.));
        else
            itf->second++;

    }

    for (map<int, double>::iterator it=hist.begin(); it!=hist.end(); it++)
        it->second=it->second*freq;

    prints(hist, out);

}

