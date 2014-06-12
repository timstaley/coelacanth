#include <iostream>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

#include "getoptpp/getopt_pp.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;


struct ProgramOptionSet
{
    double drizzle_scale_factor, drizzle_pixel_fraction;
    string input_list_path;
    string output_path;
};

ProgramOptionSet get_default_options()
{
    ProgramOptionSet defaultopts;
    defaultopts.drizzle_scale_factor=1.0;
    defaultopts.drizzle_pixel_fraction=0.45;
    defaultopts.output_path = "./drizcombine_";
    return defaultopts;
}

void print_usage_and_exit()
{
    ProgramOptionSet def_opts = get_default_options();
    cerr<<"Usage, with [-optional_arg default_value]:\n"

        <<"drizzle_combine "
        <<"[-s " << def_opts.drizzle_scale_factor << "] "
        <<"[-p " << def_opts.drizzle_pixel_fraction << "] "
        <<"[-o " << def_opts.output_path << "] "
        <<"input_frame_list.txt\n"

        <<" \nOptions:\n"
        <<"\t"<<"-s scale_factor\n"
        <<"\t"<<"-p pixel_frac\n"
        <<"\t"<<"-o output_basepath\n"
        <<endl;
    std::exit(0);
}

ProgramOptionSet load_options_from_command_line(int argc, char** argv)
{
    using namespace GetOpt;
    ProgramOptionSet opts = get_default_options();
    try{
        GetOpt_pp gopt(argc, argv);
        gopt.exceptions(std::ios::failbit | std::ios::eofbit);
        if (gopt >> OptionPresent('h', "help")){
            print_usage_and_exit();
        }
        gopt >> GetOpt::Option('s', "scale",
                opts.drizzle_scale_factor, opts.drizzle_scale_factor);
        gopt >> GetOpt::Option('p', "pixfrac",
                opts.drizzle_pixel_fraction, opts.drizzle_pixel_fraction);
        gopt >> GetOpt::Option('o', "out",
                opts.output_path,opts.output_path);
        gopt >> GetOpt::GlobalOption(opts.input_list_path);
    }
    catch(...){
        cerr << "Error in arguments" <<endl;
        print_usage_and_exit();
    }
    return opts;
}


int main(int argc, char** argv)
{
    ProgramOptionSet opts = load_options_from_command_line(argc,argv);
    cout<<"Drizzling with scale factor "<< opts.drizzle_scale_factor
        <<" and pixel fraction " << opts.drizzle_pixel_fraction <<endl;
    cout<<"Output path: " << opts.output_path<<endl;
    cout<<"Loading frames listed at: " << opts.input_list_path<<endl;
    return 0;
}
