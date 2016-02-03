#include "cast.h"

/*
bool cast_string_to_double (std::string &b, double &h)
{
    // set h= the number written in b[];
    // return false if there is an error
    h=0;
    if(b.size()==0)
        return false;

    int sign=1;
    if (b[0]=='-')
    {
        b[0]='0';
        sign=-1;
    }

    unsigned int digits_before=0;
    for(unsigned int i=0; i<b.size(); i++)
        if(b[i]!='.')
            digits_before++;
        else
            break;

    unsigned int j=0;
    while (j!=digits_before)
    {

        int number=(int(b[j])-48);
        h+=number*pow(10, digits_before-j-1);
        if (number<0 || number>9)
            return false;
        j++;
    }

    j=digits_before+1;
    while (j<b.size())
    {
        int number=(int(b[j])-48);
        h+=number*pow(10, digits_before-j);

        if (number<0 || number>9)
            return false;
        j++;
    }
    h=sign*h;
    return true;
}
*/

bool cast_string_to_double (std::string &s, double &r)
{
    const char *p = s.c_str();
    r = 0.0;
    bool neg = false;
    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (*p >= '0' && *p <= '9') {
        r = (r*10.0) + (*p - '0');
        ++p;
    }
    if (*p == '.') {
        double f = 0.0;
        int n = 0;
        ++p;
        while (*p >= '0' && *p <= '9') {
            f = (f*10.0) + (*p - '0');
            ++p;
            ++n;
        }
        r += f / std::pow(10.0, n);
    }
    if (neg) {
        r = -r;
    }
    return true;
}

int cast_int(double u)
{
    int a=int(u);
    if (u - a > 0.5)
        a++;
    return a;
}


int cast_string_to_char(std::string &file_name, char *b)
{
    for (unsigned int i=0; i<file_name.size(); i++)
        b[i]=file_name[i];
    b[file_name.size()]='\0';
    return 0;
}
