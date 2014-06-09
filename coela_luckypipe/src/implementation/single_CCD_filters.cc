/*
 * File:   single_CCD_filters.cc
 * Author: ts337
 *
 * Created on 22 February 2011, 17:53
 */

#include "../single_CCD_filters.h"
#include "../image_cleanup.h"
#include "../gain_utils.h"
#include "coela_core/src/image_utils.h"
#include "coela_core/src/pixel_array_routines.h"
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <stdexcept>
#include <utility>

using std::cerr; using std::cout; using std::endl;
using std::runtime_error;

namespace coela {
namespace single_CCD_filters {
int filter_debug=0;
//=========================================================================================================


void frame_cleanup_token::clear()
{
    frame_inf = FrameInfo();
    //bmp = image<double>(); //Actually, don't blank the bmp. By re-assigning it we can avoid de-allocating and re-allocating memory chunks
    frame_is_good=true;

//    original_bias_pedestal=-1e10; //Best I can think of for an "error" value
    cosmic_rays_confirmed_due_to_edge_case=false;
    cosmic_rays_rejected_spatially=false;
    cosmic_rays_rejected_temporally=false;
    cosmic_ray_candidates.clear();


}
//=========================================================================================================
//Begin template definitions:


//=========================================================================================================


File_Timestamp_Imprinting_Filter::File_Timestamp_Imprinting_Filter(
    const FrameInfo& first_frame, const double milliseconds_per_frame_):
    filter(serial_in_order), //NB made serial so std::gmtime does not produce memory conflicts from separate threads (it assigns to a single static memory)
    milliseconds_per_frame(milliseconds_per_frame_)
{
    init_time = boost::filesystem::last_write_time(first_frame.file_path);
    init_frame_number = first_frame.corrected_frame_id;
}

void* File_Timestamp_Imprinting_Filter::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);

    double milliseconds_since_init = milliseconds_per_frame* (double)(
                                         b.frame_inf.corrected_frame_id - init_frame_number);
    using namespace boost::posix_time;
    using namespace boost::gregorian;

    tm * ptm;
    ptm =std::gmtime(&init_time);
    ptime init_ptime =ptime_from_tm(*ptm);

    date frame_date = init_ptime.date();

//             string date_string = string_utils::find_and_replace( to_iso_extended_string(frame_date), "-" ," " );
    string date_string = to_iso_extended_string(frame_date);

    b.fht.add_keyword("UTCFDATE", date_string,"Estimated from file timestamps");

    time_duration frame_time_of_day = init_ptime.time_of_day();

//             cerr<<"Time double "<< milliseconds_since_init<<endl;
//             cerr<<"Add time: "<< to_simple_string( nanosec(milliseconds_since_init*1000))<<endl;
    frame_time_of_day += microsec(milliseconds_since_init*1000);

//             cerr<<"Frame time: "<< to_simple_string(frame_time_of_day)<<endl;

    string time_string = to_simple_string(frame_time_of_day);
//             string time_string = string_utils::find_and_replace( to_simple_string(frame_time_of_day),":" ," " );

    b.fht.add_keyword("UTCFTIME",time_string,"Estimated from file timestamps");
    return &b;
}


void* Timestamp_Load_Filter::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    b.frame_inf.header_timestamp=string_utils::atoi(b.fht.get_key_value("FRAMETIM"));
    return &b;
}


void* Frame_Count_Display_Filter::operator()(void* item)
{
    if (filter_debug) { cerr<<"Frame count filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);

    if (frame_predicate_ref(b)) { counter++; }
    if (display_frame_count) {
        cout <<"\r"<<counter<<": Loaded frame named "<<b.frame_inf.derived_lcc_filename<<" from file "<<
             string_utils::pull_filename(b.frame_inf.file_path)<<"\t";
//                            cout <<counter<<": Loaded frame "<<b.frame_inf.derived_lcc_filename<<endl;
    }

    if (throw_limit) {
        if (counter > throw_limit) {
            cerr<<"\nFlag count limit reached at"<<counter <<" throwing\n"<<endl;
            throw runtime_error("Flag count limit reached, throwing");
        }
    }
    return &b;
}

template < class T>
Sequential_File_Buffer_Filter<T>::Sequential_File_Buffer_Filter(
    const vector<FrameInfo>& frame_info_vec, size_t n_frames_to_load):
    filter(serial_in_order),
    frm_inf_vec(frame_info_vec),
    next_buffer(0), bytes_off_disk(0), counter(0)
{
    token_buffers=frame_cleanup_token::initialize_token_vec<T>(lucky_n_tokens);
    token_buffer_ptrs = vector<frame_cleanup_token*>(token_buffers.size());
    for (size_t i=0; i!= token_buffers.size(); i++) {
        token_buffer_ptrs[i]=&(token_buffers[i]);
    }
    current_posn_index = 0;
    if (n_frames_to_load!=0 && n_frames_to_load < frame_info_vec.size()) {
        frm_inf_vec.resize(n_frames_to_load) ;
    }
}

template < class T>
void Sequential_File_Buffer_Filter<T>::reset()
{
    current_posn_index=0;
    next_buffer=0;
    bytes_off_disk=0;
    counter=0;
}

template < class T>
void* Sequential_File_Buffer_Filter<T>::operator()(void*)
{
//            if (filter_debug) cerr<<"Buffer filter...\n";
    frame_cleanup_token& b = *token_buffer_ptrs[next_buffer];
    next_buffer = (next_buffer+1) % lucky_n_tokens;

    //b.clear();    //NB this gets called by the initialization and decommission filters, so is unnecessary (although harmless) here.

    if (current_posn_index==frm_inf_vec.size()) {
        return NULL;
    } else {
        b.frame_inf = frm_inf_vec[current_posn_index];
        //load the header table so we know what to do next:
//        b.fht = FitsHeader(b.frame_inf.file_path, b.frame_inf.header_byte_offset);

        string file_extension=string_utils::pull_extension(b.frame_inf.file_path);

//        if (b.frame_inf.file_is_cube_FITS){
//            PixelCubeHeader cube_header(b.fht);
//            b.array_hdr = cube_header.get_xy_array2d_header();
//
//            size_t data_start_offset=b.frame_inf.header_byte_offset+b.fht.header_file_length_in_bytes()
//                    + b.array_hdr.byte_size_of_uncompressed_data()*(b.frame_inf.cube_FITS_z_index-1);
//
//            b.buffer.buffer_file_portion(b.frame_inf.file_path,
//                                         data_start_offset,
//                                         b.array_hdr.byte_size_of_uncompressed_data() );
//
//        }
//        else
        if (file_extension==".lcz") {

            b.buffer.buffer_file_portion(b.frame_inf.file_path,
                                         b.frame_inf.header_byte_offset,
                                         b.frame_inf.byte_size);
        } else if (b.frame_inf.byte_size==0 && (file_extension==".fits"
                                                || file_extension==".lcc")) {

            b.buffer.buffer_whole_file(b.frame_inf.file_path);
        } else { throw runtime_error("Sequential_File_Buffer_Filter<T>::operator() - unrecognised data type"); }


        bytes_off_disk+=b.buffer.data_vec().size();
        counter++;
        current_posn_index++;
        return &b;
    }
}

//Declare templated types for compilation:
template class Sequential_File_Buffer_Filter<frame_cleanup_token>;


void* Decompress_Filter::operator()(void* item)
{
    if (filter_debug) { cerr<<"Decompress filter..."; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
//     if (b.frame_inf.file_is_cube_FITS) {
//         b.raw_image=image<int>::load_image_from_buffered_cube_slice(b.buffer, b.fht, b.array_hdr, 0);
//     }
//     else {

    b.fht = FitsHeader(b.frame_inf.file_path, b.frame_inf.header_byte_offset);
    b.img=CCDImage<float>::load_image_from_buffered_data(b.buffer, b.fht,
            b.fht.header_file_length_in_bytes());
    b.img.initialize_CCD_grid_for_raw_data();
    if (filter_debug) { cerr<<"Done\n"; }
    return &b;
}

void* Frame_Output_Filter::operator()(void* item)
{
    if (filter_debug) { cerr<<"Frame output filter..."; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    if (frame_predicate_ref(b)) {
        try {
            b.img.write_to_file(output_folder+"/"+b.frame_inf.derived_lcc_filename);
            n_frames_output++;
        } catch (runtime_error& e) {
            cerr<<"Output filter caught runtime error:\n"<<e.what()
                <<"\nOutput folder was set to: "<<output_folder
                <<"\nexiting..."<<endl;
            exit(0);
        }
    }
    return &b;
}

void* Frame_Sequential_Write_Filter::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    b.img.write_to_file(output_dir+"/"+b.frame_inf.derived_lcc_filename);
    return &b;
}

void* Frame_Info_Collection_Filter::operator()(void* item)
{
    if (filter_debug) { cerr<<"List collation filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    if (b.frame_is_good) { frames.push_back(b.frame_inf); }
    return &b;
}

void* Serial_Decommission_Filter::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    b.clear(); //Set the token blank so it cannot be accidentally used until properly re-initialised (NB virtual func, works for derived class "drizzle_token" also)
    return NULL;
}





void* FrameCropFilter::operator()(void* item)
{
    if (filter_debug) { cerr<<"FrameCropFilter filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    b.img=CCDImage<float>::sub_image(b.img, crop_rgn);
    return &b;
}

void* BiasDriftTracker::operator()(void* item)
{
    if (filter_debug) { cerr<<"Temporal debias filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    if (b.raw_data_zeroing_histogram_space.pixel_count()) {
        b.raw_data_zeroing_histogram_space.clear();
    }
    if (b.raw_data_zeroing_histogram_space.minimum_has_been_set()) {
        b.raw_data_zeroing_histogram_space.set_min_value(0);
    }



    b.frame_inf.bias_pedestal =
        image_cleanup::determine_bias_pedestal_from_box_in_raw_image(
            b.img.pix,
            hist_region_box_,
            b.raw_data_zeroing_histogram_space
        );
//            //Fixme: temporary hard coding of limits:
//            if (b.original_bias_pedestal<500.0 || b.original_bias_pedestal> 1200.0){
//                b.frame_is_good=false;
//                b.bias_pedestal_anomaly=true;
//            }
    return &b;
}

//void* Float_To_Double_Filter::operator()( void* item ){
//    if (filter_debug) cerr<<"Float to Double filter...\n";
//    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
//    b.debiased_input= image<double>(b.raw_input);
//    return &b;
//}


//---------------------------------------------------------------------------------------------------------------


//=========================================================================================================
//---------------------------------------------------------------------------------------------------------------

Pixel_Time_Series_Record::Pixel_Time_Series_Record(
    const vector<PixelIndex>& pixels_to_watch,
    const size_t n_frames_estimate): tbb::filter(serial_in_order)
{
    pixel_indices=pixels_to_watch;
    for (size_t i=0; i!=pixels_to_watch.size(); ++i) {
        pixval_series.push_back(vector<pixel_event>()) ;
        pixval_series.back().reserve(n_frames_estimate);
    }
}

void* Pixel_Time_Series_Record::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);

    for (size_t i=0; i!=pixel_indices.size(); ++i) {
        pixval_series[i].push_back(
            pixel_event(b.frame_inf.corrected_frame_id, b.img.pix(pixel_indices[i]))

        );
    }

    return &b;
}

using std::pair;
vector< std::pair< PixelIndex, vector<Pixel_Time_Series_Record::pixel_event> > >
Pixel_Time_Series_Record::get_time_series() const
{
    vector< pair< PixelIndex, vector<pixel_event> > > v;
    for (size_t i=0; i!=pixel_indices.size(); ++i) {
        v.push_back(
            pair<PixelIndex, vector<pixel_event> >(pixel_indices[i], pixval_series[i])
        );
    }
    return v;
}

//==================================================================================================
Col_Histogram_Gather::Col_Histogram_Gather(const PixelRange& image_layout,
        int y_range_low, int y_range_high):
    filter(serial_in_order), y_lo(y_range_low), y_hi(y_range_high)
{
    col_histograms =  vector<HistogramContainer14bit>(image_layout.x_dim());
}

void Col_Histogram_Gather::set_hist_min_value(const int minval)
{
    for (size_t i=0; i!= col_histograms.size(); ++i) {
        col_histograms[i].set_min_value(minval);
    }
}

void* Col_Histogram_Gather::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    CCDImage<float>& img=b.img;


    PixelRange img_col_hist_rgn(1,y_lo, img.pix.range().x_dim(), y_hi);
    for (PixelIterator pix(img_col_hist_rgn); pix!=pix.end; ++pix) {
        col_histograms[pix.x-1].count_rounded_value(img.pix(pix));
    }

    //possibly better optimized - do the cols 8 at a time for better output caching - probably doesn't work now I think about it.
    //input in groups of 8 cols for better row caching.
//    for (int x0=1; x0<=img.pix.range().x_dim(); x0+=8){
//        int xhi = min(x0+7ul, img.pix.range().x_dim());
//        pixel_box img_col_hist_rgn(x0,y_lo, xhi, y_hi);
//        for (PixelIterator pix(img_col_hist_rgn); pix!=pix.end; ++pix){
//            col_histograms[pix.x-1].count_rounded_value(image(pix));
//        }
//    }
    return &b;
}

//==================================================================================================
RowHistogram_Gather::RowHistogram_Gather(const PixelRange& image_layout,
        int x_range_low):
    filter(serial_in_order), x_lo_(x_range_low)
{
    row_histograms =  vector<HistogramContainer14bit>(image_layout.y_dim());
}

void RowHistogram_Gather::set_hist_min_value(const int minval)
{
    for (size_t i=0; i!= row_histograms.size(); ++i) {
        row_histograms[i].set_min_value(minval);
    }
}

void* RowHistogram_Gather::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    CCDImage<float>& img=b.img;


    PixelRange img_row_hist_rgn(img.pix.range());
    img_row_hist_rgn.low.x = x_lo_;

    for (PixelIterator pix(img_row_hist_rgn); pix!=pix.end; ++pix) {
        row_histograms[pix.y-1].count_rounded_value(img.pix(pix));
    }
    return &b;
}

//==================================================================================================
void* UniformDebias::operator()(void* item)
{
    if (filter_debug) { cerr<<"Uniform DB filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    b.img.pix-=b.frame_inf.bias_pedestal;
//    b.fht.add_comment("Uniform bias pedestal removed. ");
    return &b;
}


void* BiasFrameSubtractor::operator()(void* item)
{
    if (filter_debug) { cerr<<"Static DB filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    b.img.pix-=bias_frame.pix;
    b.fht.add_comment("Col debiased using bias frame ");
    return &b;
}


void* Histogram14BitBuildFilter::operator()(void* item)
{
    if (filter_debug) { cerr<<"Hist record filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    if (b.frame_is_good) { gain_utils::append_histogram_data_from_float_bitmap_region(b.img, hist_box_, histogram_); }
    return &b;
}

//void* Histogram10BitBuildFilter::operator()( void* item ){
//    if (filter_debug) cerr<<"Hist record filter...\n";
//    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
//    if (b.frame_is_good)    gain_utils::append_histogram_data_from_float_bitmap_region(b.img, hist_box_, histogram_);
//    return &b;
//}

using psf_models::reference_psf;
CrossCorrelator::
CrossCorrelator(
    const vector<CCD_BoxRegion>& gs_regions_,
    const reference_psf& kernel,
    const double input_resample_factor,
    const double convolution_threshold_factor,
    const double ray_proximity_limit_in_pixels,
    const bool use_parabola_fit_to_xcorr):
    filter(parallel),
    use_parabola_fitting_(use_parabola_fit_to_xcorr),
    ray_proximity_limit_(ray_proximity_limit_in_pixels)
{
    for (size_t i=0; i!=lucky_n_tokens; ++i) {
        reg_data_copies_[i].analysis_resample_factor=input_resample_factor;
        reg_data_copies_[i].convolution_threshold=convolution_threshold_factor;
        reg_data_copies_[i].gs_regions=gs_regions_;
        reg_data_copies_[i].ref_psf=kernel;
    }
}

void* CrossCorrelator::operator()(void* item)
{
    if (filter_debug) { cerr<<"Gs reg filter...\n"; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);


    if (!b.frame_inf.guide_star_estimates.empty()) {
        throw runtime_error("Frame " + b.frame_inf.derived_lcc_filename
                            +" appears to already be registered, buffer id is " +string_utils::itoa(
                                b.token_buffer_id));
    }
//            assert(b.frame_inf.guide_star_estimates.empty());
    registration_data_struct& rd=reg_data_copies_[b.token_buffer_id];

    assert(!rd.gs_regions.empty());


    for (size_t star_num=0; star_num!=rd.gs_regions.size(); star_num++) {
//        cerr<<"Guiding on star "<<star_num<<endl;
        frame_registration::gs_lock<CCD_Position> star_posn_est =
            frame_registration::find_best_psf_match(
                b.img,
                rd.gs_regions[star_num],
                rd.ref_psf,
                rd.analysis_resample_factor,
                rd.convolution_threshold,
                use_parabola_fitting_
            );

        b.frame_inf.guide_star_estimates.push_back(star_posn_est);

        for (size_t i=0; i!=b.frame_inf.confirmed_cosmic_rays.size(); ++i) {
            const CCD_Position& ray = b.frame_inf.confirmed_cosmic_rays[i];
            double ray_to_star_lock_dist = coord_distance(ray,
                                           b.frame_inf.guide_star_estimates.back().Position);
            if (ray_to_star_lock_dist < ray_proximity_limit_) {
                b.frame_is_good=false;
                return &b;
            }
        }
    }


    return &b;
}


void* Frame_Summation::operator()(void* item)
{
    if (filter_debug) { cerr<<"Sum filter..."; }
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    if (b.frame_is_good) {
        sum.pix+=b.img.pix;
        n_frames_summed++;
    }
    if (filter_debug) { cerr<<"Done\n"; }
    return &b;
}

CCDImage<double> Frame_Summation::Avg()
{
    CCDImage<double> avg(sum);
    avg.pix/=(double)n_frames_summed;
    return avg;
}

void* Dynamic_Row_Debias::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    image_cleanup::row_10p_debias_preserving_bg_level(&(b.img.pix)) ;
    return &b;
}

void* Region_Centroid::operator()(void* item)
{
    frame_cleanup_token& b = *static_cast<frame_cleanup_token*>(item);
    frames.push_back(b.frame_inf);
    PixelPosition frm_centroid = pixel_array_routines::centroid(b.img.pix , GS_box);
    CCD_Position ccd_centroid = b.img.CCD_grid.corresponding_grid_Position(frm_centroid);
    frames.back().guide_star_estimates.clear();
    frames.back().guide_star_estimates.push_back(
        frame_registration::gs_lock<CCD_Position>(ccd_centroid, 1.0)
    );
    return &b;
}




void* Cosmic_Ray_Candidate_Detection::operator()(void* item)
{
    if (filter_debug) { cerr<<"Cosmic ray candidate detection filter...\n"; }
    frame_cleanup_token& ft = *static_cast<frame_cleanup_token*>(item);
    for (PixelIterator i(ft.img.pix.range()); i!=i.end; ++i) {
        if (ft.img.pix(i) > threshold) {
            ft.cosmic_ray_candidates.push_back(i);
        }
    }
    if (!ft.cosmic_ray_candidates.empty()) { n_frames_with_candidates++; }
    return &ft;
}


void* Cosmic_Ray_Confirmation::operator()(void* item)
{
    frame_cleanup_token& ft = *static_cast<frame_cleanup_token*>(item);

    vector<PixelIndex> rays_to_check = ft.cosmic_ray_candidates;
    while (!rays_to_check.empty()) {

        Pixel<float> max_pix(rays_to_check.front(), ft.img.pix(rays_to_check.front()));
        bool candidate_rejected=false;

        for (size_t i=0; i!=rays_to_check.size(); ++i) {
            if (ft.img.pix(rays_to_check[i]) > max_pix.value) {
                max_pix = Pixel<float>(rays_to_check[i], ft.img.pix(rays_to_check[i]));
            }
        }

        if (ft.frame_inf.corrected_frame_id != (prev_frame_inf.corrected_frame_id+1)
                || first_frame) {
            ft.cosmic_rays_confirmed_due_to_edge_case=
                true; //edge case, assumed cosmic rays to be on safe side
        } else if (prev_bmp.pix(max_pix) *temporal_reduction_factor > ft.img.pix(max_pix)) {
            //for a cosmic ray we expect a much steeper increase from frame to frame, so then previous val * 10 < current val
            ft.cosmic_rays_rejected_temporally=true;
            candidate_rejected=true;
        }


        PixelRange bordered_rgn = PixelRange::pad(ft.img.pix.range(), -3);
        if (!bordered_rgn.contains_pixel(max_pix)) {
            ft.cosmic_rays_confirmed_due_to_edge_case=true;
        } else if ((ft.img.pix(max_pix.x,
                               max_pix.y +2) * spatial_reduction_factor) > ft.img.pix(max_pix)
                   && (ft.img.pix(max_pix.x, max_pix.y - 2)*spatial_reduction_factor)> ft.img.pix(max_pix)) {
            ft.cosmic_rays_rejected_spatially = true;
            candidate_rejected=true;
        }


        if (candidate_rejected==false) {
            ft.frame_inf.confirmed_cosmic_rays.push_back(
                ft.img.CCD_grid.corresponding_grid_Position(
                    PixelPosition::centre_of_pixel(max_pix))
            );
        }

        PixelRange local_window(max_pix, max_pix);
        local_window = PixelRange::pad(local_window, 20);
        vector<PixelIndex> rays_not_in_vicinity;
        for (size_t i=0; i!=rays_to_check.size(); ++i) {
            if (! local_window.contains_pixel(rays_to_check[i])) {
                rays_not_in_vicinity.push_back(rays_to_check[i]);
            }
        }
        rays_to_check=rays_not_in_vicinity;

    }

    if (!ft.frame_inf.confirmed_cosmic_rays.empty()) { n_frames_with_confirmed_rays++; }

    prev_frame_inf = ft.frame_inf;
    prev_bmp=ft.img;
    if (first_frame) { first_frame=false; }
    return &ft;
}



void* Cosmic_Ray_Diagnosis_Output::operator()(void* item)
{
    if (filter_debug) { cerr<<"Cosmic_Ray_Diagnosis_Output filter...\n"; }
    frame_cleanup_token& ft = *static_cast<frame_cleanup_token*>(item);
    if (!ft.cosmic_ray_candidates.empty()) {

        try {
            PixelRange candidates_window(ft.cosmic_ray_candidates.front(),
                                         ft.cosmic_ray_candidates.front());
            for (size_t i=1; i!=ft.cosmic_ray_candidates.size(); ++i) {
                PixelIndex p = ft.cosmic_ray_candidates[i];
                candidates_window = PixelRange::stretch_to_pixel(candidates_window, p);
            }
            candidates_window = PixelRange::pad(candidates_window, 20);


            PixelRange confirmed_rays_window;
            if (!ft.frame_inf.confirmed_cosmic_rays.empty()) {
                PixelIndex first_ray_pixel =
                    PixelPosition::pixel_centred_at(ft.img.CCD_grid.corresponding_pixel_Position(
                                                        ft.frame_inf.confirmed_cosmic_rays.front()));
                confirmed_rays_window = PixelRange(first_ray_pixel, first_ray_pixel);

                for (size_t i=1; i!=ft.frame_inf.confirmed_cosmic_rays.size(); ++i) {
                    PixelIndex p = PixelPosition::pixel_centred_at(
                                       ft.img.CCD_grid.corresponding_pixel_Position(ft.frame_inf.confirmed_cosmic_rays[i])
                                   );
                    confirmed_rays_window = PixelRange::stretch_to_pixel(confirmed_rays_window, p);
                }
                confirmed_rays_window = PixelRange::pad(confirmed_rays_window, 20);

            }

            if (!ft.frame_inf.confirmed_cosmic_rays.empty()
                    && !ft.cosmic_rays_confirmed_due_to_edge_case) {
                CCDImage<float> window_bmp=
                    CCDImage<float>::sub_image(
                        ft.img,
                        PixelRange::overlap(confirmed_rays_window, ft.img.pix.range())
                    );

                window_bmp.write_to_file(output_folder+"/confirmed_clearly/"
                                         +ft.frame_inf.derived_lcc_filename);
                confirmed_peak_vals.push_back(window_bmp.pix(window_bmp.pix.max_PixelIndex()));
            } else if (!ft.frame_inf.confirmed_cosmic_rays.empty()
                       && ft.cosmic_rays_confirmed_due_to_edge_case) {
                CCDImage<float> window_bmp=
                    CCDImage<float>::sub_image(ft.img, PixelRange::overlap(confirmed_rays_window,
                                               ft.img.pix.range()));

                window_bmp.write_to_file(output_folder+"/confirmed_due_to_edge_case/"
                                         +ft.frame_inf.derived_lcc_filename);
                confirmed_peak_vals.push_back(window_bmp.pix(window_bmp.pix.max_PixelIndex()));
            } else if (ft.cosmic_rays_rejected_spatially && !ft.cosmic_rays_rejected_temporally) {
                CCDImage<float> window_bmp=
                    CCDImage<float>::sub_image(ft.img, PixelRange::overlap(candidates_window,
                                               ft.img.pix.range()));

                window_bmp.write_to_file(output_folder+"/rejected_spatially/"
                                         +ft.frame_inf.derived_lcc_filename);
            } else if (!ft.cosmic_rays_rejected_spatially && ft.cosmic_rays_rejected_temporally) {
                CCDImage<float> window_bmp=
                    CCDImage<float>::sub_image(ft.img, PixelRange::overlap(candidates_window,
                                               ft.img.pix.range()));

                window_bmp.write_to_file(output_folder+"/rejected_temporally/"
                                         +ft.frame_inf.derived_lcc_filename);
            } else if (ft.cosmic_rays_rejected_spatially && ft.cosmic_rays_rejected_temporally) {
                CCDImage<float> window_bmp=
                    CCDImage<float>::sub_image(ft.img, PixelRange::overlap(candidates_window,
                                               ft.img.pix.range()));

                window_bmp.write_to_file(output_folder+"/rejected_both/"
                                         +ft.frame_inf.derived_lcc_filename);
            }

            n_frames_output++;
        } catch (runtime_error& e) {
            cerr<<"Cosmic_Ray_Diagnosis_Output_Filter caught runtime error:\n"<<e.what()
                <<"\nOutput folder was set to: "<<output_folder
                <<"\nexiting..."<<endl;
            exit(0);
        }
    }
    return &ft;

}


//---------------------------------------------------------------------------------------------------------------


} //end namespace coela::single_CCD_filters
}//end namespace coela



