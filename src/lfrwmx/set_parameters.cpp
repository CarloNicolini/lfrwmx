#include "set_parameters.h"

Parameters::Parameters()
{

    num_nodes=unlikely_value;
    average_k=unlikely_value;
    max_degree=unlikely_value;

    tau_degree=2;
    tau_commsize=1;

    mixing_parameter_topological=unlikely_value;
    mixing_parameter_weights=unlikely_value;

    beta=1.5;

    overlapping_nodes=0;
    overlap_membership=0;

    nmin=unlikely_value;
    nmax=unlikely_value;

    randomf=false;
    fixed_range=false;
    excess=false;
    defect=false;

    clustering_coeff=unlikely_value;
    verbosity=0;

    command_flags.push_back("-N");			//0
    command_flags.push_back("-k");			//1
    command_flags.push_back("-maxk");		//2
    command_flags.push_back("-mut");		//3
    command_flags.push_back("-taudegree");			//4
    command_flags.push_back("-taucommsize");			//5
    command_flags.push_back("-minc");		//6
    command_flags.push_back("-maxc");		//7
    command_flags.push_back("-on");			//8
    command_flags.push_back("-om");			//9
    command_flags.push_back("-beta");		//10
    command_flags.push_back("-muw");		//11
    command_flags.push_back("-C");			//12
    command_flags.push_back("-verb");			//13
}

/**
 * @brief Parameters::print
 */
void Parameters::print_parameters()
{
    cerr << "Num Nodes=" <<
    "num_nodes "<< num_nodes << endl <<
    "average_k "<<average_k << endl <<
    "max_degree "<< max_degree << endl <<
    "tau_degree "<< tau_degree << endl <<
    "tau_commsize "<<tau_commsize << endl <<
    "mixing_parameter_topological "<< mixing_parameter_topological << endl <<
    "mixing_parameter_weights " << mixing_parameter_weights << endl <<
    "beta "<< beta << endl <<
    "overlapping_nodes " << overlapping_nodes << endl <<
    "overlap_membership "<< overlap_membership << endl <<
    "nmin" << nmin << endl <<
    "nmax" << nmax << endl <<
    "fixed_range" << fixed_range << endl <<
    "excess"<< excess << endl <<
    "defect" << defect << endl <<
    "randomf" << randomf << endl <<
    "clustering_coeff=" << clustering_coeff << endl;
}

/**
 * @brief Parameters::set_random
 */
void Parameters::set_random()
{
    FILE_LOG(logINFO)<<"This is a random network";
    mixing_parameter_topological=0;
    mixing_parameter_weights=0;
    overlapping_nodes=0;
    overlap_membership=0;
    nmax=num_nodes;
    nmin=num_nodes;

    fixed_range=true;
    excess=false;
    defect=false;
}

/**
 * @brief Parameters::arrange
 * @return
 */
int Parameters::arrange()
{
    if(randomf)
        set_random();

    if (num_nodes==unlikely_value)
    {
        throw std::logic_error("ERROR: number of nodes unspecified");
        FILE_LOG(logERROR) << "ERROR: number of nodes unspecified";
        return false;
    }

    if (average_k==unlikely_value)
    {
        throw std::logic_error("ERROR: average degree unspecified");
        FILE_LOG(logERROR) << "ERROR: average degree unspecified";
        return false;
    }

    if (max_degree==unlikely_value)
    {
        throw std::logic_error("ERROR: maximum degree unspecified");
        FILE_LOG(logERROR) << "ERROR: maximum degree unspecified";
        return false;
    }

    if(mixing_parameter_topological==unlikely_value)
    {
        mixing_parameter_topological=mixing_parameter_weights;
    }

    if (mixing_parameter_weights==unlikely_value)
    {
        throw std::logic_error("ERROR: weight mixing parameter (option -muw) unspecified");
        FILE_LOG(logERROR) << "ERROR: weight mixing parameter (option -muw) unspecified";
        return false;
    }


    if(mixing_parameter_topological > 1 || mixing_parameter_topological <0 )
    {
        throw std::logic_error("ERROR:mixing parameter > 1 (must be between 0 and 1)");
        FILE_LOG(logERROR)<<"ERROR:mixing parameter > 1 (must be between 0 and 1)";
        return false;
    }

    if(overlapping_nodes<0 || overlap_membership<0)
    {
        throw std::logic_error("ERROR:some positive parameters are negative");
        FILE_LOG(logERROR)<<"ERROR:some positive parameters are negative";
        return false;
    }

    if (num_nodes<=0 || average_k<=0 || max_degree<=0 || mixing_parameter_topological<0 || mixing_parameter_weights<0 || (nmax<=0 && nmax!=unlikely_value) || (nmin<=0 && nmin!=unlikely_value) )
    {
        throw std::logic_error("ERROR:some positive parameters are negative");
        FILE_LOG(logERROR)<<"ERROR:some positive parameters are negative";
        return false;
    }


    if(mixing_parameter_weights > 1 || mixing_parameter_weights <0 )
    {
        throw std::logic_error("ERROR:mixing2 parameter > 1 (must be between 0 and 1)");
        FILE_LOG(logERROR)<<"ERROR:mixing2 parameter > 1 (must be between 0 and 1)";
        return false;
    }

    if(nmax!= unlikely_value && nmin!=unlikely_value)
        fixed_range=true;
    else
        fixed_range=false;

    if(excess && defect)
    {
        throw std::logic_error("ERROR:both options -inf and -sup cannot be used at the same time");
        FILE_LOG(logERROR)<<"ERROR:both options -inf and -sup cannot be used at the same time";
        return false;
    }

    FILE_LOG(logINFO)<<"number of nodes:"<<num_nodes;
    FILE_LOG(logINFO)<<"average degree:"<<average_k;
    FILE_LOG(logINFO)<<"maximum degree:"<<max_degree;
    FILE_LOG(logINFO)<<"exponent for the degree distribution:"<<tau_degree;
    FILE_LOG(logINFO)<<"exponent for the community size distribution:"<<tau_commsize;
    FILE_LOG(logINFO)<<"mixing parameter (topology):"<<mixing_parameter_topological;
    FILE_LOG(logINFO)<<"mixing parameter (weights):"<<mixing_parameter_weights;
    FILE_LOG(logINFO)<<"beta exponent:"<<beta;
    FILE_LOG(logINFO)<<"number of overlapping nodes:"<<overlapping_nodes;
    FILE_LOG(logINFO)<<"number of memberships of the overlapping nodes:"<<overlap_membership;

    if(clustering_coeff!=unlikely_value)
    {
        FILE_LOG(logINFO) << "Average clustering coefficient: "<<clustering_coeff;
    }

    if (fixed_range)
    {
        FILE_LOG(logINFO)<<"Community size range set equal to ["<<nmin<<" , "<<nmax<<"]";
        if (nmin>nmax)
        {
            throw std::logic_error("ERROR: Inverted community size bounds");
            FILE_LOG(logERROR)<<"ERROR: Inverted community size bounds";
            return false;
        }

        if(nmax>num_nodes)
        {
            throw std::logic_error("ERROR: maxc bigger than the number of nodes");
            FILE_LOG(logERROR)<<"ERROR: maxc bigger than the number of nodes";
            return false;
        }
    }
    return true;
}

/**
 * @brief Parameters::set
 * @param flag
 * @param num
 * @return
 */
bool Parameters::set(string & flag, string & num)
{
    // False if something goes wrong
    FILE_LOG(logINFO)<<"Setting... " << flag << " " << num;

    double value;
    if(!cast_string_to_double(num, value))
    {
        throw std::logic_error("ERROR while reading parameters");
        FILE_LOG(logERROR) << "ERROR while reading parameters";
        return false;
    }

    if (flag==command_flags[0])
    {
        if (fabs(value-int(value))>1e-8)
        {
            throw std::logic_error("ERROR: number of nodes must be an integer");
            FILE_LOG(logERROR) << "ERROR: number of nodes must be an integer";
            return false;
        }
        num_nodes=cast_int(value);
    }
    else if(flag==command_flags[1])
    {
        average_k=value;
    }
    else if(flag==command_flags[2])
    {
        max_degree=cast_int(value);
    }
    else if(flag==command_flags[3])
    {
        mixing_parameter_topological=value;
    }
    else if(flag==command_flags[11])
    {
        mixing_parameter_weights=value;
    }
    else if(flag==command_flags[10])
    {
        beta=value;
    }
    else if(flag==command_flags[4])
    {
        tau_degree=value;
    }
    else if(flag==command_flags[5])
    {
        tau_commsize=value;
    }

    else if(flag==command_flags[6])
    {
        if (fabs(value-int (value))>1e-8)
        {
            throw std::logic_error("ERROR: the minumum community size must be an integer");
            FILE_LOG(logERROR)<<"ERROR: the minumum community size must be an integer";
            return false;
        }
        nmin=cast_int(value);
    }
    else if(flag==command_flags[7])
    {
        if (fabs(value-int (value))>1e-8)
        {
            throw std::logic_error("ERROR: the maximum community size must be an integer");
            FILE_LOG(logERROR) << "ERROR: the maximum community size must be an integer";
            return false;
        }
        nmax=cast_int(value);
    }
    else if(flag==command_flags[8])
    {
        if (fabs(value-int (value))>1e-8)
        {
            throw std::logic_error("ERROR: the number of overlapping nodes must be an integer");
            FILE_LOG(logERROR) << "ERROR: the number of overlapping nodes must be an integer";
            return false;
        }
        overlapping_nodes=cast_int(value);
    }
    else if(flag==command_flags[9])
    {
        if (fabs(value-int (value))>1e-8)
        {
            throw std::logic_error("ERROR: the number of membership of the overlapping nodes must be an integer");
            FILE_LOG(logERROR) << "ERROR: the number of membership of the overlapping nodes must be an integer";
            return false;
        }
        overlap_membership=cast_int(value);
    }
    else if(flag==command_flags[12])
    {
        clustering_coeff=value;
    }
    else
    {
        stringstream ss; ss << "ERROR while reading parameters: "<<flag<<" is an unknown option";
        throw std::logic_error(ss.str());
        //FILE_LOG(logERROR) << ss;
        return false;
    }
    return true;
}

/**
 * @brief print_usage
 */
void print_usage()
{

    FILE_LOG(logINFO)<<"\nTo run the program type \n./benchmark [FLAG] [P]";
    FILE_LOG(logINFO)<<"\n----------------------\n";
    FILE_LOG(logINFO)<<"To set the parameters, type:"<<endl;
    FILE_LOG(logINFO)<<"-N[number of nodes]";
    FILE_LOG(logINFO)<<"-k[average degree]";
    FILE_LOG(logINFO)<<"-maxk[maximum degree]";
    FILE_LOG(logINFO)<<"-mut[mixing parameter for the topology]";
    FILE_LOG(logINFO)<<"-muw[mixing parameter for the weights]";
    FILE_LOG(logINFO)<<"-beta[exponent for the weight distribution]";
    FILE_LOG(logINFO)<<"-t1[minus exponent for the degree sequence]";
    FILE_LOG(logINFO)<<"-t2[minus exponent for the community size distribution]";
    FILE_LOG(logINFO)<<"-minc[minimum for the community sizes]";
    FILE_LOG(logINFO)<<"-maxc[maximum for the community sizes]";
    FILE_LOG(logINFO)<<"-on[number of overlapping nodes]";
    FILE_LOG(logINFO)<<"-om[number of memberships of the overlapping nodes]";
    FILE_LOG(logINFO)<<"-C[Average clustering coefficient]";
    FILE_LOG(logINFO)<<"----------------------\n";
    FILE_LOG(logINFO)<<"It is also possible to set the parameters writing flags and relative numbers in a file. To specify the file, use the option:";
    FILE_LOG(logINFO)<<"-f[filename]";
    FILE_LOG(logINFO)<<"You can set the parameters both writing some of them in the file, and using flags from the command line for others."<<endl;
    FILE_LOG(logINFO)<<"-N, -k, -maxk, -muw have to be specified. For the others, the program can use default values:";
    FILE_LOG(logINFO)<<"t1=2, t2=1, on=0, om=0, beta=1.5, mut=muw, minc and maxc will be chosen close to the degree sequence extremes.";
    FILE_LOG(logINFO)<<"If you don't specify -C the rewiring process for raising the average clustering coefficient will not be performed";
    FILE_LOG(logINFO)<<"If you set a parameter twice, the latter one will be taken.";
    FILE_LOG(logINFO)<<"\n-------------------- Other options ---------------------------\n";
    FILE_LOG(logINFO)<<"To have a random network use:";
    FILE_LOG(logINFO)<<"-rand";
    FILE_LOG(logINFO)<<"Using this option will set muw=0, mut=0, and minc=maxc=N, i.e. there will be one only community.";
    FILE_LOG(logINFO)<<"Use option -sup (-inf) if you want to produce a benchmark whose distribution of the ratio of external degree/total degree ";
    FILE_LOG(logINFO)<<"is superiorly (inferiorly) bounded by the mixing parameter.";
    FILE_LOG(logINFO)<<"\n-------------------- Examples ---------------------------\n";
    FILE_LOG(logINFO)<<"Example1:";
    FILE_LOG(logINFO)<<"./benchmark -N 1000 -k 15 -maxk 50 -muw 0.1 -minc 20 -maxc 50";
    FILE_LOG(logINFO)<<"Example2:";
    FILE_LOG(logINFO)<<"./benchmark -f flags.dat -t1 3";
    FILE_LOG(logINFO)<<"\n-------------------- Other info ---------------------------\n";
    FILE_LOG(logINFO)<<"Read file ReadMe.txt for more info."<<endl;
}

/**
 * @brief set_from_file
 * @param file_name
 * @param par1
 * @return
 */
bool set_from_file(string & file_name, Parameters & par1)
{
    int h= file_name.size();
    char b[h+1];
    cast_string_to_char(file_name, b);

    ifstream in(b);
    if (!in.is_open())
    {
        FILE_LOG(logERROR)<<"File "<<file_name<<" not found. Where is it?";
        return false;
    }

    string temp;
    while(in>>temp)  			// input file name
    {
        if(temp=="-rand")
            par1.randomf=true;
        else if(temp=="-sup")
            par1.excess=true;
        else if(temp=="-inf")
            par1.defect=true;
        else
        {
            string temp2;
            in>>temp2;
            if(temp2.size()>0)
            {
                if(temp=="-f" && temp2!=file_name)
                {
                    if(set_from_file(temp2, par1)==false)
                        return false;
                }
                if(temp!="-f")
                {
                    if(par1.set(temp, temp2)==false)
                        return false;
                }
            }
            else
            {
                FILE_LOG(logERROR)<<"Errow while reading parameters";
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief set_parameters
 * @param argc
 * @param argv
 * @param par1
 * @return
 */
bool set_parameters(int argc, char * argv[], Parameters & par1)
{
    if (argc <= 1)   // if no arguments, return statement about program usage.
    {
        print_usage();
        return false;
    }

    int argct = 0;
    string temp;
    while (++argct < argc)  			// input file name
    {
        temp = argv[argct];
        if(temp=="-rand")
            par1.randomf=true;
        else if(temp=="-sup")
            par1.excess=true;
        else if(temp=="-inf")
            par1.defect=true;
        else
        {
            argct++;
            string temp2;
            if(argct<argc)
            {
                temp2 = argv[argct];
                if(temp=="-f")
                {
                    if(set_from_file(temp2, par1)==false)
                        return false;
                }
                if(temp!="-f")
                {
                    if(par1.set(temp, temp2)==false)
                        return false;
                }
            }
            else
            {
                FILE_LOG(logERROR)<<"ERROR while reading parameters";
                return false;
            }
        }
    }

    if(par1.arrange()==false)
        return false;

    return true;
}

