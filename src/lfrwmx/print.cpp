#include "print.h"

int cherr()
{

    cerr<<"the check failed"<<endl;
    int e;
    cin>>e;
    return e;
}


int cherr(double a)
{

    cerr<<"the check failed because of "<<a<<endl;
    int e;
    cin>>e;
    return e;
}
