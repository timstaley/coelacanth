

#include <string>
#include "../file_info.h"
//#include "fits_header.h"
#include "coela_utility/src/string_utils.h"
#include "coela_utility/src/simple_serialization.h"
#include <boost/filesystem.hpp>
#include <stdexcept>
#include <iostream>
#include <fstream>
//#include <dirent.h>
#include <algorithm>
using namespace std;
namespace bfs=boost::filesystem;

namespace coela {
using namespace string_utils;

FileInfo::FileInfo():
    file_path("p"),derived_lcc_filename("p"),
    header_byte_offset(0), byte_size(0),
    file_is_cube_FITS(false),
    ccd_id(-42), header_frame_id(-1), corrected_frame_id(-1),
    folder_number(-1), subfile_number(-1), header_timestamp(-1) {}


FileInfo FileInfo::basic_fits_file(const std::string& fname)
{
    FileInfo fi; //use default constructor
    fi.file_path=fname;
    return fi;
}

std::ostream& operator<<(std::ostream& os, const FileInfo& file)
{
    os<<file.ccd_id<<","<<file.folder_number<<","<<pull_file_number_string(
          file.file_path)<<","<<file.subfile_number
      <<","<<file.header_byte_offset<<","<<
      file.byte_size<<","<<file.header_frame_id<<","<<file.header_timestamp;
    return os;
}

void FileInfo::load_key_vals_from_header_table(const FitsHeader& kvc_table)
{
    ccd_id=atoi(kvc_table.get_key_value("CAMERAID"));
    header_timestamp = atol(kvc_table.get_key_value("FRAMETIM"));
    header_frame_id = atoi(kvc_table.get_key_value("FRAMENUM"));
}


void FileInfo::get_subfiles_info(std::list<string>& lcm_list,
                                  std::list<FileInfo>& subfiles, const std::string& extension)
{
    if (extension==".lcm"||extension=="lcm") { get_lcm_subfiles_info(lcm_list, subfiles); }
    else if (extension==".lcz") { subfiles = lcz_utils::build_lcz_subfile_index(lcm_list); }
    else { throw runtime_error("get_subfiles_info - extension not recognised"); }
}
void FileInfo::get_lcm_subfiles_info(const std::list<string>& lcm_list,
                                      std::list<FileInfo>& subfiles)
{
    //iterate through the lcm filenames
    for (list<string>::const_iterator iter = lcm_list.begin(); iter!=lcm_list.end(); ++iter) {

        vector<size_t> offsets(parse_lcm_header(
                                   *iter));    //this is where the work is really done - A vector of byte offsets is returned
        FileInfo subfile;
        size_t total_files_in_lcm=
            offsets.size();           //this line mostly just improves readability

        for (vector<size_t>::size_type j=0; j<total_files_in_lcm; ++j) {
            subfile.file_path=*iter;
            subfile.derived_lcc_filename=generate_unpacked_lcc_filename(*iter, j, total_files_in_lcm);
//                subfile.file_num_in_lcm= j;
            subfile.total_files_in_lcm=total_files_in_lcm;
            subfile.header_byte_offset= offsets[j];
            subfiles.push_back(subfile);
        }
    }
}


bool FileInfo::filename_number_predicate(const string& first_filename,
        const string& second_filename)
{
    return (string_utils::pull_file_number(first_filename) < string_utils::pull_file_number(
                second_filename));
}
bool FileInfo::subfile_number_predicate(const FileInfo& first, const FileInfo& second)
{
    return (first.subfile_number < second.subfile_number);
}

bool FileInfo::frame_id_predicate(const FileInfo& first, const FileInfo& second)
{
    return (first.corrected_frame_id < second.corrected_frame_id);
}

bool FileInfo::subfile_numbers_equal(const FileInfo& first, const FileInfo& second)
{
    return (first.subfile_number == second.subfile_number);
}








vector<size_t> FileInfo::parse_lcm_header(const string& filename)
{
    vector<size_t> byte_offset_values;
    string line;
    int num_files;
    ifstream lcm(filename.c_str());
    try {
        if (lcm.is_open()) {
            getline(lcm, line);  //Standard "LuckyCam mult filebuf" description - (ignore)
            getline(lcm, line);  // "X Frames"

            num_files=string_utils::atoi(line.substr(0, line.find(' ')));        //Returns X

            for (int i=0; i<num_files; ++i) {
                getline(lcm, line);   //e.g. "3 12345"
                byte_offset_values.push_back(string_utils::atoi(line.substr(
                                                 line.find(' '))));       //Returns "12345" part
            }
        }

        else { throw std::runtime_error(filename); }
        lcm.close();
    } catch (std::runtime_error fname)
    {cerr << "Fatal Error accessing file: " <<fname.what() <<endl; exit(0); }

    return byte_offset_values;
}

std::list<FileInfo> FileInfo::get_image_file_list(const std::string& dir_in,
        const std::string& filestem, int CCD_id, const std::string& extension)
{
    if (extension.find("lcz")!=string::npos) { cout <<"Warning - this file lister does not check for lcz indices."<<endl; }
    string dir(dir_in);
    if (dir[dir.size()-1]!='/') { dir+='/'; }
    //First list actual files
    cout <<"Getting file list..."<<flush;
    list<string> matching_files;
    matching_files = file_utils::get_matching_files_recursively(dir,
                     extension); //returns a list of filenames matching the name criteria
//        matching_files.sort(file_utils::compare_filenames); //This may(?) actually improve speed because they tend to be stored on the disk in order
    if (matching_files.empty()) { throw runtime_error("No matching files for: "+ dir_in +"; "+filestem +"; "+extension); }

    //Now get subfiles info if applicable
    list<FileInfo> files_info;

    if ((extension==".lcm" || extension=="lcm" ||extension ==".lcz"))  {
        FileInfo::get_subfiles_info(matching_files, files_info,
                                     extension); //if the files are LCM, interrogate the headers to create subfile info
    } else {
        //dealing with regular fits files, assume they have no index file and no header data on frame id, etc.
        //Assign header ordering stamps, (assuming that alphabetically sorting the files == sort into chronological order)
        size_t file_counter=1;
        matching_files.sort(FileInfo::filename_number_predicate);
        for (list<string>::iterator iter=matching_files.begin(); iter!=matching_files.end();
                ++iter) {
            //otherwise just initialise each FileInfo as normal.
            if (iter->find(filestem)==0 || iter->find('/'+filestem)!=string::npos) {
                //ASSUMPTION - directory markers are fwd slash ("/", unix!)
                FileInfo f=FileInfo::basic_fits_file(*iter);
                f.subfile_number=file_counter;
                f.header_frame_id=file_counter;
                f.corrected_frame_id=file_counter;
                f.derived_lcc_filename="file_"+itoa(file_counter)+".fits";
                f.ccd_id=CCD_id;
                file_counter++;
                files_info.push_back(f);  //initialises each lcm_subfile with the filename string.
            }
        }
    }

    //Check if the FITs files are in fact cube files
    //currently assumes we have a homogenous dataset.
    FitsHeader temp_table(files_info.front().file_path,
                          files_info.front().header_byte_offset);
    if (temp_table.key_exists("NAXIS3")
            && string_utils::atoi(temp_table.get_key_value("NAXIS3"))!=0) {
        for (list<FileInfo>::iterator iter=files_info.begin(); iter!=files_info.end(); ++iter) {
            iter->file_is_cube_FITS=true;
        }
    }


    cout <<"done"<<endl;
    FileInfo::convert_paths_to_system_complete(files_info);
//        files_info.sort();
    return files_info;
}


string FileInfo::generate_unpacked_lcc_filename(const string& filename, int file_num,
        int files_in_lcm)   //filename is e.g. dir/blah_X.lcm
{
    string end_string(filename.substr(filename.find_last_of('_') +
                                      1));         //returns "X.lcm"
    int end_num(string_utils::atoi(end_string.substr(0,
                                   end_string.find('.'))));             //returns X
    string base_name(filename.substr(0,
                                     filename.find_last_of('_')));           //returns dir/blah
    return base_name + "_" + itoa(end_num - (files_in_lcm-1) + file_num) +
           ".fits";                //returns dir/blah_(X+y).lcc
}

int FileInfo::read_num_files_in_lcm(const string& lcm_filename)
{
    string line;
    int num_files;
    ifstream lcm(lcm_filename.c_str());
    try {
        if (lcm.is_open()) {
            getline(lcm, line);  //Standard "LuckyCam mult filebuf" description - (ignore)
            getline(lcm, line);  // "X Frames"
            num_files=string_utils::atoi(string(line.substr(0,
                                                line.find(' '))).c_str());               //Returns X
        }

        else { throw runtime_error(lcm_filename); }
        lcm.close();
    } catch (runtime_error fname)
    {cerr << "Fatal Error accessing file: " <<fname.what() <<endl; exit(0); }
    return num_files;
}

FileInfo FileInfo::load_from_lcz_index_line(const std::string& index_line,
        const std::string& base_folder,
        const std::string& ccd_stem, const std::string& subfolder_stem,
        const std::string& filestem,
        const std::string& file_extension)
{
    vector<string> line_segments = tokenize(index_line,",");
    return load_from_lcz_index_line_segments(line_segments, base_folder, ccd_stem,
            subfolder_stem, filestem, file_extension);
}

FileInfo FileInfo::load_from_lcz_index_line_segments(const vector<string>&
        line_segments,
        const std::string& base_folder,
        const std::string& ccd_stem, const std::string& subfolder_stem,
        const std::string& filestem,
        const std::string& file_extension)
{
    assert(line_segments.size()==8);
    FileInfo subfile;
    subfile.ccd_id = atoi(line_segments[0]);
    subfile.folder_number = atoi(line_segments[1]);
    string lcz_number=line_segments[2];
    subfile.subfile_number = atoi(line_segments[3]);
    subfile.header_byte_offset= atoi(line_segments[4]);
    subfile.byte_size =atol(line_segments[5]);
    subfile.header_frame_id = atoi(line_segments[6]);
    subfile.header_timestamp = atol(line_segments[7]);

    subfile.file_path = base_folder + ccd_stem+"_"+itoa(subfile.ccd_id) + "/"
                        +subfolder_stem+"_"+itoa(subfile.folder_number)+"/"
                        +filestem+"_"+itoa(subfile.ccd_id)+"_"+lcz_number+file_extension;

    subfile.derived_lcc_filename = filestem+"_"+itoa(subfile.ccd_id)+"_"
                                   +itoa(subfile.subfile_number)+".fits";

    return subfile;
}

void FileInfo::convert_paths_to_system_complete(std::list<FileInfo>& file_list)
{
    for (list<FileInfo>::iterator it=file_list.begin(); it!=file_list.end(); ++it) {
        it->file_path = boost::filesystem::system_complete(it->file_path).string();
    }
}

///NB: This lists files in the directory passed as an argument, whose filenames include 'pattern' string and end with 'extension' string. It does not list contents of subdirectories.
/*list<string> file_utils::list_directory_contents(const string& dir, const string& pattern,const string& extension)
{
   list<string> files_in_dir;

    DIR *dp;
    struct dirent *ep;

    dp = opendir (dir.c_str());
    if (dp != NULL)
    {
        string fname;
        while ( (ep = readdir (dp))  ){
         fname=ep->d_name;
          // if (extension pattern we're looking for is actually AT THE END of the filename...) && (extension pattern exists in filename) && (stem exists in filename)
         if ( ( fname.size() - fname.rfind(extension) == extension.size())  &&  (    fname.rfind(extension)!=string::npos   )  && (fname.find(pattern) != string::npos) )
         {  files_in_dir.push_back(dir + fname); }
        }

        (void) closedir (dp);
    }
    else
        throw runtime_error("Couldn't open the directory: "+ dir);

    return files_in_dir;
}
*/



namespace lcz_utils {

list<FileInfo> build_lcz_subfile_index(const std::list<string>& lcz_list)
{
    list<FileInfo> subfiles;
    int counter=0;
    for (list<string>::const_iterator iter = lcz_list.begin(); iter!=lcz_list.end(); ++iter) {
        cout <<"\rParsing file "<<++counter<<" of " <<lcz_list.size()<<": " << *iter <<"\t\t";
        list<FileInfo> lcz_files= parse_lcz_file(*iter);
        subfiles.splice(subfiles.end(), lcz_files);
    }
//            subfiles= parse_lcz_file(lcz_list.front()); //debug version
    return subfiles;
}

bool empty_files_present(const std::list<string>& lcz_list)
{
    list<FileInfo> subfiles;
    for (list<string>::const_iterator iter = lcz_list.begin(); iter!=lcz_list.end(); ++iter) {
        if (bfs::file_size(*iter)==0) { return true; }
    }
    return false;
}

list<FileInfo> get_lcz_subfile_list(const std::string& ccd_dir,
                                     const std::string& index_output_dir)
{
    string lcz_extension=".lcz";
    list<string> matching_filenames = file_utils::get_matching_files_recursively(ccd_dir,
                                      lcz_extension);

    list<FileInfo> lcc_subfiles;

    //-----------------------------------------------------------------------------------------------------------
    //Calculate the index filename and base dir:

    bfs::path absolute_ccd_dir = bfs::system_complete(ccd_dir);
    string ccd_id = absolute_ccd_dir.string();
    //knock any trailing '/' off the path so the next parent operation is correct
    if (ccd_id[ccd_id.size()-1]=='/') { absolute_ccd_dir=absolute_ccd_dir.parent_path(); }
    ccd_id = absolute_ccd_dir.string();
//            if (ccd_id.find("[")!=string::npos){
//                ccd_id = ccd_id.substr(ccd_id.find_last_of('['));
//                ccd_id = strip_lead_tail_chars(ccd_id, "[]");
//            }
//            else
    ccd_id = ccd_id.substr(ccd_id.find_last_of('_')+1);

    bfs::path target_base_dir = absolute_ccd_dir.parent_path();

    cout <<"Basedir: \""<<target_base_dir<<"\""<<endl;
    string target_name = target_base_dir.filename().string();
    cout <<"Run: \""<<target_name<<"\""<<endl;
    string target_num = target_name.substr(0, target_name.find('_'));
//            cout <<target_num<<endl;
    cout <<"CCD: \""<<ccd_id<<"\""<<endl;

    string index_filename = "/run"+target_num+"_Index_"+ ccd_id +".log";

    //-----------------------------------------------------------------------------------------------------------
    //Determine if the index is present

    string index_path, data_dir;
    list<FileInfo> subfile_list;
    cout <<"Looking for index at "<<index_output_dir+index_filename<<endl;
    if (bfs::exists(index_output_dir+index_filename)) {
        //Check for index in output dir first
        cout <<"Index found at " << index_output_dir+index_filename<<endl;
        index_path = index_output_dir
                     +index_filename; //assume an index in the output folder is valid!
    } else if (bfs::exists(target_base_dir.string()+index_filename)) {
        //Then check in base dir of data
        cout <<"Index found at " << target_base_dir.string()+index_filename<<endl;
        if (check_lcz_index_is_standard_format(target_base_dir.string()+index_filename)) {
            index_path =target_base_dir.string()+index_filename;
        } else { cout<<"Non-standard format, rebuild."<<endl; }

        cout <<"Checking for empty files..."<<flush;
        if (empty_files_present(matching_filenames)) {
            cout <<" empty files present, index probably corrupt, will rebuild."<<endl;
            index_path.clear();
        } else { cout <<" Done"<<endl; }
    }

    //-----------------------------------------------------------------------------------------------------------
    //Proceed to load the index, or generate one...

    if (index_path.empty()) {
        //if no index exists
        cout <<"No index file found, building index..."<<endl;
        subfile_list=build_lcz_subfile_index(matching_filenames);
        subfile_list.sort(FileInfo::subfile_number_predicate);
        if (subfile_list.empty()) { throw runtime_error("No valid files found for this folder: "+absolute_ccd_dir.string()); }
        cout <<"Writing new index to file \""<<index_output_dir+index_filename<<"\"";
        write_lcz_index(subfile_list, index_output_dir+index_filename);
        cout <<" ...done"<<endl;
    } else {
        subfile_list=parse_lcz_index(index_path, target_base_dir.string()+"/");
    }
    return subfile_list;
}


string mangle_lcz_filename(const string& filename)
{
    string name=filename;
    name = coela::string_utils::find_and_replace(name,"]","_");
    name = coela::string_utils::find_and_replace(name,"[","_");
    name = coela::string_utils::find_and_replace(name,"__","_");
    name = coela::string_utils::find_and_replace(name,"_.",".");
    return name;
}







std::list<int> get_lcz_nums(const vector<string>& text)
{
    list<int> lcz_list;
    int prev_file_num, this_file_num;
    //line[0] is a header
    vector<string> first_line_segments=tokenize(text[1],",");
    assert(first_line_segments.size()==8);
    prev_file_num = atoi(first_line_segments[2]);
    lcz_list.push_back(prev_file_num); //lcz number 0
    for (size_t i=2; i!=text.size(); ++i) {
        vector<string> line_segments=tokenize(text[i],",");
        this_file_num =atoi(line_segments[2]);
        if (this_file_num!=prev_file_num) {
            lcz_list.push_back(this_file_num);
            prev_file_num = this_file_num;
        }
    }
    return lcz_list;
}




std::list<int> get_lcz_nums(const std::string& index_filename)
{
    cout <<"Getting indexed lcz nums..."<<flush;
    ifstream index_file(index_filename.c_str());
    if (!index_file.is_open()) { throw runtime_error("Can't open index at "+index_filename); }
    vector<string> text = simple_serialization::convert_istream_to_string_vec(index_file);
    if (text.size()<2) { throw runtime_error("Not a valid index file:" + index_filename); }
    index_file.close();
    cout <<"Done"<<endl;
    return get_lcz_nums(text);
}

void write_lcz_index(const std::list<FileInfo>& subfiles, const string& index_filename)
{
    assert(!subfiles.empty());

    bfs::path file_path(bfs::system_complete(subfiles.front().file_path));
    string file_path_str = file_path.string();
    cout <<"\nFirst file path:"<<file_path_str<<endl;
    string filestem = pull_filestem(file_path_str);

    if (file_path_str.find('[')!=string::npos) {
        //temp modification to deal with problematic early files.
        bfs::path base_folder = file_path.parent_path().parent_path().parent_path();
        write_lcz_index(subfiles, index_filename, base_folder.string(), "[]", "[]", filestem);
    } else {
        bfs::path sub_folder = file_path.parent_path();
        cout <<"subfolder:" <<sub_folder<<endl;
        string sub_folder_stem = sub_folder.filename().string();
        sub_folder_stem = sub_folder_stem.substr(0,sub_folder_stem.find('_'));
        cout <<"subfolderstem: " <<sub_folder_stem<<endl;
        bfs::path cam_folder = sub_folder.parent_path();
        string cam_folder_stem = cam_folder.filename().string();
        cam_folder_stem = cam_folder_stem.substr(0,cam_folder_stem.find('_'));
        cout <<"Camfolder stem" <<cam_folder_stem<<endl;
        bfs::path base_folder = cam_folder.parent_path();
        cout <<"Basefolder:"<<base_folder<<endl;
        write_lcz_index(subfiles, index_filename, base_folder.string(),
                        cam_folder_stem, sub_folder_stem, filestem);
    }

}
void write_lcz_index(const std::list<FileInfo>& subfiles, const string& index_filename,
                     const string& root_folder, const string& cam_folder_stem, const string& sub_folder_stem,
                     const string& filestem)
{
    ofstream index_file(index_filename.c_str());
    if (!index_file) { throw runtime_error("Could not write to file: "+index_filename); }
    index_file
            <<root_folder<<","<<cam_folder_stem<<","<<sub_folder_stem<<","<<filestem<<",lcz\n"; //nb root folder actually ignored, just a placeholder
    for (list<FileInfo>::const_iterator it=subfiles.begin(); it!=subfiles.end(); ++it) {
        index_file <<*it<<"\n";
    }
}

std::list<FileInfo> parse_lcz_file(const std::string& lcz_filename)
{
    ifstream input_file(lcz_filename.c_str(), ifstream::binary |ifstream::ate);
    if (!input_file.is_open()) {throw runtime_error("Could not open file: \""+lcz_filename+"\"");}
    if (!input_file) {throw runtime_error("Could not read file: "+lcz_filename);}
    size_t byte_size = input_file.tellg();
    input_file.seekg(0);

    char* buffer_char_ptr = new char[byte_size];
    input_file.read((char*)buffer_char_ptr, byte_size);
    input_file.close();

    uint32_t * buffer_int_ptr = (unsigned int *)(buffer_char_ptr);
    char start_text[5]="<<<<";
    char end_text[5]=">>>>";
    uint32_t start_int = *(uint32_t*)&start_text;
    uint32_t end_int = *(uint32_t*)&end_text;

    list<FileInfo> index;

    //NB Here, "header" and "tail" refer to <<<<blah.lcc   >>>>> (NOT the fits header.)

    //Fixme: Once we've found the header, we could advance 'compressedbits'*n_pixels mod whatever. (after checking the data is there)
    //But this would require reading the FITS Header, so there is a tradeoff.

    uint32_t *header_text_begin, *header_text_end, *tail_text_begin, *tail_text_end,
             *buffer_end;
    buffer_end = buffer_int_ptr + byte_size/sizeof(uint32_t);
    tail_text_end= buffer_int_ptr;
    while (tail_text_end != buffer_end) {

        header_text_begin = find(tail_text_end, buffer_end, start_int);

        header_text_end = 1 + find(header_text_begin, buffer_end,
                                   end_int)  ; //NB Marks the next 4 bytes AFTER ">>>>" , so +1
        size_t subfile_offset = (header_text_end - buffer_int_ptr)*sizeof(
                                    uint32_t); //NB Offset is recorded in bytes

        tail_text_begin =find(header_text_end, buffer_end, start_int);

        //check that we really have found the next marker, and not just some confusing data, by checking that the first 4 bytes after <<<< match with the header.
        //NB Implemented this after coming across some data with exactly that random byte sequence...
        //NB This is still not bulletproof, but chances of failure should be astronomically small. (Could still get a random occurrence matching 8 header bytes)
        while (*(tail_text_begin+1) != *(header_text_begin+1) && tail_text_begin!=buffer_end) {
            //if it doesn't match, try searching again
            tail_text_begin = find(tail_text_begin +1, buffer_end,
                                   start_int) ; //+1 to start searching from the next 4 bytes
        }
        tail_text_end = 1 + find(tail_text_begin, buffer_end, end_int);

        if ((header_text_end < header_text_begin) || (tail_text_end < tail_text_begin)) {
            throw logic_error("Symbol mismatch in file "+lcz_filename);
        }
        string header_text((char*)header_text_begin, (char*)header_text_end),
               tail_text((char*)tail_text_begin, (char*)tail_text_end);

        if (header_text!=tail_text) {
            throw runtime_error("LCZ File error: "+ lcz_filename+ " -header does not match tail: "
                                + header_text + " / " + tail_text +"\n");
        }
        string::size_type name_start = header_text.find_first_not_of('<');
        string subfile_name = header_text.substr(name_start,
                              header_text.find_last_not_of(" >") - name_start + 1);

        FileInfo subfile=FileInfo::basic_fits_file(lcz_filename);
        //get rid of any square bracketed filenames
        if (subfile_name.find('[')!=string::npos) { subfile_name = lcz_utils::mangle_lcz_filename(subfile_name); }
        subfile.derived_lcc_filename= subfile_name;
        subfile.header_byte_offset = subfile_offset;
        subfile.byte_size = (char*)tail_text_begin - (char*)header_text_end;
        string folder_number_str =lcz_filename.substr(0, lcz_filename.find_last_of('/'));
        folder_number_str= folder_number_str.substr(folder_number_str.find_last_of('_')+1);
        subfile.folder_number = string_utils::atoi(folder_number_str);
        subfile.subfile_number=string_utils::pull_file_number(subfile.derived_lcc_filename);
        FitsHeader temp_hdr(subfile.file_path, subfile.header_byte_offset);
        subfile.load_key_vals_from_header_table(temp_hdr);
        index.push_back(subfile);
    }
    delete[] buffer_char_ptr;
    return index;
}

//deprecated:
/*
void get_lcz_subfiles_info(std::list<string>& lcm_list, std::list<FileInfo>& subfiles){
    string ccd_folder_stem, subfolder_stem;
    string target_root_folder = lcm_list.front(); //probably contains /blah/target/Camera_0/Folder_0/Image_0_0.lcz (but could be Folder0/image_0_1.lcz if already in Cam dir.)
    cout <<"First file path: " <<target_root_folder<<endl;
    target_root_folder = target_root_folder.substr(0, target_root_folder.find_last_of('/')); // /blah/target/Camera_0/Folder_0 (Strip filename)
//        target_root_folder = target_root_folder.substr(0, target_root_folder.find_last_of('/')); // /blah/target/Camera_0
//        target_root_folder = target_root_folder.substr(0, target_root_folder.find_last_of('/')); // /blah/target

    if (target_root_folder.find("Folder")!=string::npos){
        target_root_folder = target_root_folder.substr(0, target_root_folder.find("Folder"));  ///blah/target/Camera_0 or blank
        subfolder_stem ="Folder";
    }
    else throw runtime_error("Non-standard subfolder naming");
    cout <<"Cam folder resolved to: "<<target_root_folder<<endl;
//
    if (target_root_folder.find("Camera")!=string::npos){
        target_root_folder = target_root_folder.substr(0, target_root_folder.find("Camera"));  ///blah/target/
        ccd_folder_stem="Camera";
    }
    else if (target_root_folder=="./"){}//That's ok. but will need to check dir above for index
    else throw runtime_error("Non-standard CCD folder naming");
    if (target_root_folder.empty()) target_root_folder="./";
    cout <<"Root folder resolves to " <<target_root_folder<<endl;


    string lcz_path = lcm_list.front(); // contains /blah/target/Camera_0/Folder_0/Image_0_0.lcz
    string filestem = string_utils::pull_filestem(lcz_path);
    string camera_id_str=string_utils::pull_camera_id_string(lcz_path);

    string index_filename(target_root_folder+filestem+"_Index_" + camera_id_str +".log");
    string alt_index_filename(target_root_folder+"../"+filestem+"_Index_" + camera_id_str +".log");

    cout <<"Attempting to find index file at:\n"<<index_filename <<" ,  "<<alt_index_filename<<"../"<<endl;

    ifstream index_file(index_filename.c_str(), ios::ate | ios::binary);
    list<int> indexed_file_nums;
    bool index_present(false), valid_index(false);
    if (index_file.is_open() && index_file.tellg()!=0) {
        cout <<"Found index file at " <<index_filename<<endl;
        index_file.close();
        indexed_file_nums = get_lcz_nums(index_filename);
        index_present=true;
    }
    else {
        index_file.open( alt_index_filename.c_str() ,ios::ate | ios::binary);
        if (index_file.is_open() && index_file.tellg()!=0) {
            cout <<"Found index file at " <<alt_index_filename<<endl;
            index_filename=alt_index_filename;
            index_file.close();
            indexed_file_nums = get_lcz_nums(index_filename);
            index_present=true;
        }

    }

    //FIXME:: should really do a more robust check, but this will do for now...
    if (indexed_file_nums.size()==lcm_list.size()) {valid_index=true; cout<<"Index appears to be ok"<<endl;}
    else if (indexed_file_nums.size() < lcm_list.size() && index_present ) {
        //Often there is an extra copy of a penultimate lcz
        cout <<"Scanning " <<lcm_list.back() <<" to determine if it is a mis-copy to be ignored"<<endl;
        list<string>::const_iterator it=lcm_list.end(); --it;
        list<FileInfo> suspect_lcz_subfiles= parse_lcz_file( *it);
        --it;
        list<FileInfo> one_from_last= parse_lcz_file(*it);
        --it;
        list<FileInfo> two_from_last= parse_lcz_file(*it);

        if (suspect_lcz_subfiles.front().derived_lcc_filename==two_from_last.front().derived_lcc_filename){
            cout <<lcm_list.back() <<" Is a duplicate of "<< two_from_last.front().file_path<<", ignoring..."<<endl;
            lcm_list.pop_back();
            valid_index=true;
        }
        else if (suspect_lcz_subfiles.front().derived_lcc_filename==one_from_last.front().derived_lcc_filename ){
            cout <<lcm_list.back() <<" Is a duplicate of "<< one_from_last.front().file_path<<", ignoring..."<<endl;
            lcm_list.pop_back();
            valid_index=true;
        }
        else if ( find(indexed_file_nums.begin(), indexed_file_nums.end(), (int)string_utils::pull_file_number(lcm_list.back()) ) == indexed_file_nums.end()   ){
            //Sometimes there is simply an extra, un-indexed file
            cout <<lcm_list.back() <<" Appears to be an un-indexed file. Will add to list."<<endl;
            subfiles.splice(subfiles.end(), suspect_lcz_subfiles);
            valid_index =true;
        }
    }
    else cout <<"Index filesize does not match file-list ( "<< indexed_file_nums.size() <<" != "<< lcm_list.size() <<" ), rebuilding..."<<endl;

    if (valid_index) {
        cout <<"Loading subfile info from index file...";
        list<FileInfo> indexed_files = parse_lcz_index(index_filename);
        subfiles.splice(subfiles.end(), indexed_files);

        assert(!subfiles.empty());
        cout <<"Done"<<endl;
    }
    else{
        cout <<"Building index file" <<endl;
        for (list<string>::const_iterator iter = lcm_list.begin(); iter!=lcm_list.end(); ++iter){
            cout <<"\rParsing " << *iter <<"...            ";
            list<FileInfo> lcz_files= parse_lcz_file(*iter);
            subfiles.splice(subfiles.end(), lcz_files);
        }
        write_lcz_index(subfiles, index_filename,
            target_root_folder, ccd_folder_stem, subfolder_stem, filestem );
    }
    size_t before_check_unique=subfiles.size();
    subfiles.sort(FileInfo::subfile_number_predicate);
    subfiles.unique(FileInfo::subfile_numbers_equal);
    if (subfiles.size()!= before_check_unique) throw runtime_error("Most troubling - duplicate subfile numbers detected");
}
*/
bool check_lcz_index_is_standard_format(const string& index_filename)
{
    ifstream index_file(index_filename.c_str());
    vector<string> text = simple_serialization::convert_istream_to_string_vec(index_file);
    index_file.close();
    //Load header line
    if (text.size()<2) { throw runtime_error("Unknown index format: "+index_filename); }
    vector<string> line_segments;
    line_segments = tokenize_and_strip_spaces(text.front(),",");
    if (line_segments.size()!=5) { throw runtime_error("Unknown index format: "+index_filename); }
    line_segments = tokenize_and_strip_spaces(text[1],",");
    if (line_segments.size()==8) { return true; }
    else { return false; }
}


std::list<FileInfo> parse_lcz_index(const string& index_filename,
                                     const string& data_root_dir)
{
    using namespace string_utils;
    list<FileInfo> lcc_files;
    ifstream index_file(index_filename.c_str());
    if (!index_file) { throw runtime_error("Can't open index for parsing at" + index_filename); }
    vector<string> text = simple_serialization::convert_istream_to_string_vec(index_file);
    index_file.close();

    //Load header line
    vector<string> line_segments;
    line_segments = tokenize_and_strip_spaces(text.front(),",");
    assert(line_segments.size()==
           5); //first line format: <<root_folder(originally...)<<",Camera,Folder,"<<filestem<<",.lcz\n";
    string cam_stem(line_segments[1]), subfolder_stem(line_segments[2]),
           filestem(line_segments[3]),
           file_extension(line_segments[4]);

    bfs::path absolute_path = bfs::system_complete(index_filename);
    string root_folder = absolute_path.parent_path().string()+"/";
    //if the index file is stored elsewhere than the data:
    if (!data_root_dir.empty()) { root_folder=data_root_dir; }
    if (!file_extension.empty() && file_extension[0]!='.') { file_extension= "."+file_extension; }

    list<int> lcz_filenums_present = get_lcz_nums(text);
    bool reached_last_lcz=false;
    cout <<"Loading index file with root folder " <<root_folder<<endl;
    for (size_t i=1; i<text.size(); ++i) {
        //(start at 1 to avoid header line)
        FileInfo indexed_subfile = FileInfo::load_from_lcz_index_line(text[i],
                                    root_folder,
                                    cam_stem, subfolder_stem, filestem, file_extension
                                                                       );
        string lcz_file_num=string_utils::pull_file_number_string(indexed_subfile.file_path);
        //Deal with a bug in Frank's Pixcel output - sometimes the last lcz file is in fact a copy of the last but one
        //(so the last few frames get dropped, and the file indexes are wrong)
        if (!reached_last_lcz && atoi(lcz_file_num)==lcz_filenums_present.back()) {
            //last lcz file in the run.
            cout <<"Checking if last file is a miscopy..."<<flush;
            list<FileInfo> last_lcz_subfiles = parse_lcz_file(indexed_subfile.file_path);
            reached_last_lcz=true;
            if (last_lcz_subfiles.front().derived_lcc_filename!=
                    indexed_subfile.derived_lcc_filename) {
                //If the index information is wrong (buggy file)
                cout <<"LCZ file "
                     <<indexed_subfile.file_path<<" is a mis-copy, will not bother loading "<<endl;
                cout <<"( Found "<< last_lcz_subfiles.front().derived_lcc_filename <<" ; expected " <<
                     indexed_subfile.derived_lcc_filename <<" )" <<endl;
                text.erase(text.begin()+i, text.end()); //remove all references to this last file.
            } else {
                lcc_files.push_back(indexed_subfile);
                cout<<"clean"<<endl;
            }
        } else { lcc_files.push_back(indexed_subfile); }
    }
    return lcc_files;
}


bool check_headers_match_index(const std::list<FileInfo>& files)
{
    bool headers_match=true;
    for (list<FileInfo>::const_iterator it=files.begin(); it!=files.end(); ++it) {
        cout <<"\rChecking " <<it->derived_lcc_filename<<"\t\t"<<flush;
        FitsHeader hdr_tbl(it->file_path, it->header_byte_offset);
        FileInfo hdr_inf(*it);
        hdr_inf.load_key_vals_from_header_table(hdr_tbl);
        if (hdr_inf.header_frame_id!= it->header_frame_id) {
            cerr<<"Frame ID mismatch for subfile " <<it->subfile_number <<endl; headers_match=false;
        }
        if (hdr_inf.ccd_id != it->ccd_id) {
            cerr<<"CCD ID mismatch for subfile " <<it->subfile_number <<endl; headers_match=false;
        }
        if (hdr_inf.header_timestamp != it->header_timestamp) {
            cerr<<"header_timestamp mismatch for subfile " <<it->subfile_number
                <<" index/header " <<it->header_timestamp<<" / " << hdr_inf.header_timestamp
                <<", diff: "<< it->header_timestamp - hdr_inf.header_timestamp <<endl;
            headers_match=false;
        }
    }
    return headers_match;
}


int check_frames_stored_sequentially(const std::list<FileInfo>& files)
{
    int frames_dropped=0;
    int prev_num = files.front().corrected_frame_id-1;
    for (list<FileInfo>::const_iterator it=files.begin(); it!=files.end(); ++it) {
        int current_num=it->corrected_frame_id;
        if (current_num!= prev_num+1) {
            cout <<"Non sequitur between frame ids " <<prev_num<<" and " <<current_num<<endl;
            cout <<"(Current file number: " <<it->subfile_number<<" )"<<endl;
            frames_dropped+= current_num - prev_num -1;
        }
        prev_num=current_num;
    }
    return frames_dropped;
}

int check_for_dropped_frames(std::list<FileInfo>& files)
{
    int frames_dropped=0;
    files.sort(FileInfo::frame_id_predicate);
    int prev_num = files.front().corrected_frame_id-1;
    for (list<FileInfo>::const_iterator it=files.begin(); it!=files.end(); ++it) {
        int current_num=it->corrected_frame_id;
        if (current_num!= prev_num+1) {
            cout <<"Frame ids missing between " <<prev_num<<" and "
                 <<current_num<<" ("<<it->header_frame_id<<")" <<endl;
            cout <<"( Current file number: " <<it->subfile_number<<", prev file num";
            --it;
            cout<< it->subfile_number<<" )"<<endl;
            ++it;
            frames_dropped+= current_num - prev_num -1;
        }
        prev_num=current_num;
    }
    return frames_dropped;
}

int check_for_duplicate_frame_ids(std::list<FileInfo>& files)
{
    int frames_duped=0;
    files.sort(FileInfo::frame_id_predicate);
    int prev_num = files.front().corrected_frame_id-1;
    for (list<FileInfo>::const_iterator it=files.begin(); it!=files.end(); ++it) {
        int current_num=it->corrected_frame_id;
        if (current_num== prev_num) {
            //                cout <<"Duplicate Frame id: " <<current_num<<endl;
            //                cout <<"(Current file number: " <<it->subfile_number<<", prev file num";
            --it;
            //                cout<< it->subfile_number<<" )"<<endl;
            ++it;
            frames_duped++;
        }
        prev_num=current_num;
    }
    return frames_duped;
}

void set_correct_frame_ids(std::list<FileInfo>& files)
{
    int loop_number=0;
    int prev_hdr_id=0;
    for (list<FileInfo>::iterator it=files.begin(); it!=files.end(); ++it) {
        if (it->header_frame_id==-1) {
            FitsHeader tmp_hdr_tbl(it->file_path, it->header_byte_offset);
            if (tmp_hdr_tbl.key_exists("FRAMENUM")) {
                it->header_frame_id = string_utils::atoi(tmp_hdr_tbl.get_key_value("FRAMENUM"));
            }
        }

        if (it->header_frame_id==0  && it!=files.begin()) {
            if (prev_hdr_id==16383) {
                loop_number++;
            } else { throw runtime_error("Non-sequitur caused correct frame id confusion"); }
        } else if (it->header_frame_id< prev_hdr_id -14000) { throw runtime_error("Non-sequitur caused correct frame id confusion"); }

        it->corrected_frame_id = it->header_frame_id + loop_number*16384; //2^14
        prev_hdr_id = it->header_frame_id;
    }
}


void remove_duplicate_frame_ids(std::list<FileInfo>& files, const bool rigorous_check)
{
    files.sort(FileInfo::frame_id_predicate);
    int prev_num = files.front().corrected_frame_id-1;
    for (list<FileInfo>::iterator it=files.begin(); it!=files.end(); ++it) {
        int current_num=it->corrected_frame_id;

        if (current_num== prev_num) {
            cout <<"Removing duplicate frame id " <<current_num<<"("<<it->header_frame_id<<")"
                 <<" stored in file nums ";
            it--;
            cout<<it->subfile_number <<" ("<<pull_file_number_string(it->file_path)<<")";
            it++;
            cout<<" , "<<it->subfile_number<<" ("<<pull_file_number_string(it->file_path)<<")"<<endl;
            if (rigorous_check) {
                cout<<"Actual file checking not yet implemented"<<endl; //FIXME
            }

            it=files.erase(it); //return the next one after the erased file_inf
            it--;  //step back one before we increment, so we don't jump over the next one.
        }
        prev_num=current_num;
    }
}

void sanitise_list(std::list<FileInfo>& files)
{
    files.sort(FileInfo::subfile_number_predicate);
    set_correct_frame_ids(files);
    int dupes = check_for_duplicate_frame_ids(files);
    if (dupes) {
        cerr<<"Removing "<<dupes<<" duplicate frames"<<endl;
        remove_duplicate_frame_ids(files);
    }
    int dropped=check_for_dropped_frames(files);
    if (dropped) { cerr<<"Missing "  <<dropped<<" frames "<<endl; }

}


} //end lcz utils


}//end namespace coela
