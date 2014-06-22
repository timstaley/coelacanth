#include <iostream>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <set>
#include<map>

#include <boost/filesystem.hpp>
#include "getoptpp/getopt_pp.h"

#include "coela_luckypipe/src/file_info.h"
#include "coela_utility/src/string_utils.h"
#include "coela_core/src/drizzle.h"

//using std::string;
//using std::cout;
//using std::cerr;
//using std::endl;
using namespace std;
using namespace coela;

struct ProgramOptionSet {
  double drizzle_scale_factor, drizzle_pixel_fraction;
  string input_directory;
  string output_folder_path;
};

ProgramOptionSet get_default_options() {
  ProgramOptionSet defaultopts;
  defaultopts.drizzle_scale_factor = 1.0;
  defaultopts.drizzle_pixel_fraction = 0.45;
  defaultopts.output_folder_path = "./drizzled/";
  return defaultopts;
}

void print_usage_and_exit() {
  ProgramOptionSet def_opts = get_default_options();
  cerr << "Usage, with [-optional_arg default_value]:\n"

  << "drizzle_combine "
       << "[-s " << def_opts.drizzle_scale_factor << "] " << "[-p "
       << def_opts.drizzle_pixel_fraction << "] " << "[-o "
       << def_opts.output_folder_path << "] " << "input_frame_list.txt\n"

       << " \nOptions:\n"
       << "\t" << "-s scale_factor\n" << "\t" << "-p pixel_frac\n" << "\t"
       << "-o output_basepath\n" << endl;
  std::exit(0);
}

ProgramOptionSet load_options_from_command_line(int argc, char** argv) {
  using namespace GetOpt;
  ProgramOptionSet opts = get_default_options();
  try {
    GetOpt_pp gopt(argc, argv);
    gopt.exceptions(std::ios::failbit | std::ios::eofbit);
    if (gopt >> OptionPresent('h', "help")) {
      print_usage_and_exit();
    }
    gopt
        >> GetOpt::Option('s', "scale", opts.drizzle_scale_factor,
                          opts.drizzle_scale_factor);
    gopt
        >> GetOpt::Option('p', "pixfrac", opts.drizzle_pixel_fraction,
                          opts.drizzle_pixel_fraction);
    gopt >> GetOpt::Option('o', "out", opts.output_folder_path, opts.output_folder_path);
    gopt >> GetOpt::GlobalOption(opts.input_directory);
  } catch (...) {
    cerr << "Error in arguments" << endl;
    print_usage_and_exit();
  }
//  Always terminate the output path with /
  opts.output_folder_path+='/';
  return opts;
}

struct DrizzleInputGroup {
  string name;
  string output_pathstem;
  list<pair<string, PixelShift> > files_and_offsets;
};

list<DrizzleInputGroup> get_drizzle_groups(const ProgramOptionSet& opts) {

//  Quick and dirty:
//  Hard-code the pixel-shift offsets by frame number:
  map<int, PixelShift> offsets;
//  Offsets as provided to me
  offsets[1] = PixelShift(-2.5, -2.5);
  offsets[2] = PixelShift(-4.5, 2.0);
  offsets[3] = PixelShift(2, 2.5);
  offsets[4] = PixelShift(4, -2);

//  Negate them, to give the desired translation:
  for (auto& index_off : offsets) {
    index_off.second *= -1.0;
  }

  list<FileInfo> input_files = FileInfo::get_image_file_list(
      opts.input_directory, "", 0, ".fits");

//  cout << "Files:" << endl;
  map<string, list<string>> file_groups;
  for (auto f : input_files) {
    string prefix = string_utils::pull_filestem(f.file_path, true);
//    cout <<"File: "<<f.file_path<<", prefix: "<< prefix<<endl;
    file_groups[prefix].push_back(f.file_path);
  }

  list<DrizzleInputGroup> driz_groups;
  for (auto& iter : file_groups) {
    DrizzleInputGroup group;
    group.name = iter.first;
//    cout <<"Group:"<<group.name<<endl;
    for (auto filepath : iter.second) {
      int file_num = int(string_utils::pull_file_number(filepath));
      group.files_and_offsets.push_back(make_pair(filepath, offsets[file_num]));
    }
//    cout <<"files: "<< group.files_and_offsets.size()<<endl;
    group.output_pathstem = opts.output_folder_path +'/'+ group.name;
    driz_groups.push_back(group);
  }
  return driz_groups;
}

int main(int argc, char** argv) {
  ProgramOptionSet opts = load_options_from_command_line(argc, argv);
  cout << "Drizzling with scale factor " << opts.drizzle_scale_factor
       << " and pixel fraction " << opts.drizzle_pixel_fraction << endl;
  cout << "Output path: " << opts.output_folder_path << endl;
  cout << "Loading frames listed at: " << opts.input_directory << endl;

  boost::filesystem::create_directories(opts.output_folder_path);
  auto driz_groups = get_drizzle_groups(opts);

  double zoom = 1.0 / opts.drizzle_scale_factor;

  for (auto grp : driz_groups) {
    cout << endl;
    cout << "Group: " << grp.name << endl;
    cout << "Output: " << grp.output_pathstem << endl;
    auto first_filepath = grp.files_and_offsets.front().first;
    auto first_img = PixelArray2d<float>::load_from_file(first_filepath);
    auto img_size = first_img.range();
    auto drizzle_weighted_vals = PixelArray2d<float>(
        img_size.x_dim() * zoom, img_size.y_dim() * zoom, 0.0);
    auto drizzle_weights(drizzle_weighted_vals);
    auto input_weights = PixelArray2d<float>(
        img_size.x_dim(),img_size.y_dim(),1.0);

    FitsHeader additional_info;
    additional_info.set_keyword("PIXFRAC",
                            string_utils::ftoa(opts.drizzle_pixel_fraction,4));
    additional_info.set_keyword("PIXSCALE",
                              string_utils::ftoa(opts.drizzle_scale_factor,4));
//    additional_info.add_comment("INPUT FILES:");

    for (auto file_offset : grp.files_and_offsets) {
      auto fpath = file_offset.first;
//      additional_info.add_comment(string_utils::pull_filename(fpath));
      auto offset = file_offset.second;
      cout << fpath << " : " << offset << endl;
      auto img =  PixelArray2d<float>::load_from_file(fpath);
      drizzle::translate_and_drizzle_frame(
        img, input_weights,
        drizzle_weighted_vals, drizzle_weights,
        offset,
        opts.drizzle_scale_factor, opts.drizzle_pixel_fraction);
    }
    auto result = drizzle::unweight_drizzle_results(drizzle_weighted_vals,
                                                    drizzle_weights);

    result.write_to_file(grp.output_pathstem+"_driz.fits",
                         additional_info);
//    drizzle_weights.write_to_file(grp.output_pathstem+"_weights.fits",
//                                  additional_info);

  }

return 0;
}
