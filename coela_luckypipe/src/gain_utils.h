/*
 * File:   gain_utils.h
 * Author: ts337
 *
 * Created on 11 April 2009, 21:32
 */

#ifndef COELA_GAIN_UTILS_H
#define COELA_GAIN_UTILS_H
#include "coela_core/src/ccd_image.h"

#include "coela_utility/src/histogram_container.h"
#include <fstream>
#include <utility>
#include <map>
#include <Minuit2/FCNBase.h>
using std::map;
namespace coela {
namespace gain_utils {
//    typedef image<double> double_bitmap;
//======================================================================================================================

struct gain_info {
    double bias_pedestal, readout_sigma;
    double gain;
    double N_dark_pix, N_light_pix, N_serial_CIC_pix;
    double photon_event_freq, serial_CIC_rate;
//    double serial_CIC_gain;
    long actual_number_pix_events_recorded;
    gain_info(); ///<initializes all members to negative values so it's obvious what is undetermined
    gain_info(const string& filename);

    void write_to_file(const std::string& filename);
    static std::string get_column_headers();
};
std::ostream& operator<<(std::ostream& os, const gain_info& g_inf);

//======================================================================================================================

//======================================================================================================================


//to do: templatize (worth the effort?)
void append_histogram_data_from_float_bitmap_region(const CCDImage<float>& bmp,
        const PixelRange& rgn, HistogramContainer14bit& hist);
void append_histogram_data_from_float_bitmap_region(const CCDImage<float>& bmp,
        const PixelRange& rgn, HistogramContainer10bit& hist);


template <class T>
map<int,T> sum_histograms(const map<int,T>& hist1, const map<int,T>& hist2);

template <class T>
void write_histogram_data_to_file(const map<int,T>& histogram, const string& filename);

template <class T>
double sum_counts_in_histogram(const map<int,T>& histogram);

template <class T>
double sum_counts_in_histogram_above_threshold(const map<int,T>& histogram,
        const int threshold);

map<int,long> load_histogram_data_from_file(const string& filename);


double histogram_differences_squared(const map<int,long>& data,
                                     const map<int,double>& model);
//        double linear_histogram_differences(const map<int,long>& data,const map<int,double>& model);
///Better for fitting power laws:
double poisson_noise_weighted_histogram_difference(const map<int,long>& data,
        const map<int,double>& model);

map<int, double> convolve_histogram_with_gaussian(const map<int, double>& model_hist,
        const double readout_sigma,
        const int kernel_half_width);
//==================================================================================
class gaussian_histogram_fit_FCN:public ROOT::Minuit2::FCNBase {

public:
    gaussian_histogram_fit_FCN(const map<int,long>& histogram):data(histogram),
        error_def(1.0) {}
    ~gaussian_histogram_fit_FCN() {}
    virtual double  Up() const {return error_def;}
    virtual double operator()(const std::vector<double>& gauss_pars)
    const;  //Params: zero value, peak_count, sigma.
    map<int,double> get_model_histogram(const std::vector<double>& gauss_pars) const;

private:
    map<int,long> data;
    double error_def;
};
//==================================================================================

class exponential_histogram_fit_FCN:public ROOT::Minuit2::FCNBase {
public:
    exponential_histogram_fit_FCN(const map<int,long>& histogram):
        error_def(1.0), data(histogram) {}
    ~exponential_histogram_fit_FCN() {}
    virtual double  Up() const {return error_def;}
    virtual double operator()(const std::vector<double>& exp_pars)
    const;  //Params:  amplitude, gain
    map<int,double> get_model_histogram(const std::vector<double>& exp_pars) const;
private:
    double error_def;
    map<int,long> data;
};

//==================================================================================
//struct full_EMCCD_low_light_level_histogram_model{
//    static map<int,double> get_full_model_histogram(
//            const double bias_pedestal,
//            const double readout_sigma,
//            const double EM_gain,
//            const double N_dark_pix,
//            const double N_CICIR_pix);
//};

class full_EMCCD_histogram_fit_FCN: public ROOT::Minuit2::FCNBase {
public:
    full_EMCCD_histogram_fit_FCN(const map<int,long>& histogram);
    full_EMCCD_histogram_fit_FCN(const map<int,long>& histogram,
                                 const int fit_region_min, const int fit_region_max);


    virtual double  Up() const {return error_def;}

    ///Params: bias_pedestal, N_dark_pix, RO_sigma, (from gaussian fit)
    ///         N_photon_events, photon_gain, (from tail fit)
    ///         N_serial_CIC_pix.
    virtual double operator()(const std::vector<double>& pars) const;
    map<int,double> get_model_histogram(const std::vector<double>& pars) const;

    static map<int,double> get_CICIR_model_histogram(const double gain,
            const double N_CICIR_pix,
            const int N_EM_stages,
            const int max_val_to_generate);

    static vector<double> generate_pars_vec(const gain_info&);


private:
    const static int N_EM_serial_register_stages=604;
    double error_def;
    map<int,long> data_;
    int fit_rgn_min_, fit_rgn_max_;
    static const size_t npars_=6;
};

//double calculate_thresholded_SNR(const double readout_noise_in_ADU,
//                                const double mean_photon_gain_in_ADU,
//                                const double CIC_event_rate_per_pixel_readout,
//                                const double light_level,
//                                const double threshold_in_photo_electrons,
//                                const bool output_to_screen=false);
//==================================================================================
class Thresholded_SNR_Calculator {
public:
    Thresholded_SNR_Calculator(const double readout_noise_in_ADU,
                               const double mean_photon_gain_in_ADU,
                               const double light_level,
                               const double CICIR_event_rate_per_pixel_readout
                              );

    double calc_SNR_at_threshold(const double threshold_in_photo_electrons) const;
    double find_optimum_threshold(
        const double search_min = 0.01,
        const double search_max = 1.0,
        const double stepsize = 0.01) const;
private:
    double gain_, light_level_;
    map<int,double> photon_distribution_, full_distribution_;
};

double calculate_ideal_photon_counter_SNR(const double light_level);

double calculate_linear_mode_SNR(const double readout_noise_in_ADU,
                                 const double mean_photon_gain_in_ADU,
                                 const double light_level,
                                 const double CICIR_event_rate_per_pixel_readout
                                );


//class ThresholdOptimizingFCN: public ROOT::Minuit2::FCNBase{
//public:
//    ThresholdOptimizingFCN( const Thresholded_SNR_Calculator& x): calc_(x){}
//
//    virtual double operator()(const std::vector<double>& pars) const{
//        assert(pars.size()==1);
//        return -1.0*calc_.calc_SNR_at_threshold(pars[0]); //NB minimizer function, so invert it to get maxima
//    };
//
//    virtual double Up() const {return 1.0;}
//
//private:
//    Thresholded_SNR_Calculator calc_;
//};

//==================================================================================

std::pair<int,long> get_mode(const map<int,long>& histogram_data);



float get_gaussian_HWHM(const map<int,long>& histogram_data);

map<int, long> get_gaussian_fitting_section(const map<int,long>& histogram_data);
map<int, long> get_tail_fitting_section(const map<int,long>& histogram_data,
                                        const double gaussian_sigma_estimate, const bool output_to_screen=false);

map<int,long> get_serial_CIC_model_fitting_section(const map<int,long>& histogram_data,
        const double bias_pedestal,
        const double readout_sigma);

map<int,long> get_full_model_fitting_section(const map<int,long>& histogram_data,
        const double bias_pedestal,
        const double readout_sigma);

map<int, long> get_points_at_values_greater_than(const map<int,long>& histogram_data,
        int zero_point);
map<int, long> threshold_above_count(const map<int,long>& histogram_data,
                                     int count_threshold);

///Returns gauss_par vector: peak value (zero value), peak count, sigma (readout sigma);
vector<double> fit_readout_gaussian(const map<int,long>& histogram_data,
                                    const bool output_to_screen=false);

///Returns amplitude, gain
vector<double> fit_exponential_tail(const map<int,long>& histogram_data,
                                    const double gaussian_sigma_estimate, const bool output_to_screen=false);
vector<double> fit_exponential_tail_fixed_gain(const map<int,long>& histogram_data,
        const double gaussian_sigma_estimate, const double gain_estimate,
        const bool output_to_screen=false);

gain_info fit_full_CCD_model(const map<int,long>& histogram_data,
                             const gain_info& approx_fit,
                             const bool output_to_screen
//,
//        const bool sqrt_serial_CIC_gain_mode=false
                            );

gain_info fit_CCD_histogram(const map<int,long>& histogram_bins,
                            const bool perform_advanced_fit=false,
                            const bool fixed_gain_mode=false, const double gain_value=0.0,
                            const bool output_to_screen=false
//        , const bool sqrt_serial_CIC_gain_mode=false
                           );

//==================================================================================


template<typename input_datatype>
CCDImage<input_datatype>& normalise_CCD_with_uniform_gain(CCDImage<input_datatype>& input,
        const gain_info& det_info);
//==================================================================================

//image<float>& normalise_CCD_with_gain_map(double_bitmap& input, const gain_info& det_info, const double_bitmap& gain_map);

CCDImage<float> threshold_bitmap(const CCDImage<float>& normalized_input,
                                 const double normalized_threshold);



///Threshold mask can be considered boolean: 1 --> Apply thresholding here, 0--> Bright region, No thresholding
CCDImage<float> create_threshold_mask(const CCDImage<float>& raw_input,
                                      const float threshold_count);


//gain_info analyse_histogram(const map<int, long>& full_data_histogram , const string& output_dir);


//define templates:
template <class T>
void write_histogram_data_to_file(const std::map<int,T>& histogram,
                                  const std::string& filename)
{
    std::ofstream datafile(filename.c_str());
    typedef typename std::map<int,T>::const_iterator ConstIterator;
    for (ConstIterator it=histogram.begin(); it!=histogram.end(); ++it) {
        datafile.precision(12);
        datafile << it->first <<"  " << it->second <<"\n";
    }
    datafile.close();
};

template <class T>
double sum_counts_in_histogram(const map<int,T>& hist)
{
    double sum=0.0;
    typedef typename std::map<int,T>::const_iterator ConstIterator;
    for (ConstIterator it=hist.begin(); it!=hist.end(); ++it) {
        sum+= it->second;
    }
    return sum;
}

template <class T>
double sum_counts_in_histogram_above_threshold(const map<int,T>& hist,
        const int threshold)
{
    double sum=0.0;
    typedef typename std::map<int,T>::const_iterator ConstIterator;
    for (ConstIterator it=hist.begin(); it!=hist.end(); ++it) {
        if (it->first > threshold) { sum+= it->second; }
    }
    return sum;

}

template <class T>
map<int,T> sum_histograms(const map<int,T>& hist1, const map<int,T>& hist2)
{
    map<int,T>merged(hist1);
    typedef typename std::map<int,T>::const_iterator ConstIterator;
    for (ConstIterator it=hist2.begin(); it!=hist2.end(); ++it) {
        if (merged.count(it->first)==0) { merged[it->first]=it->second; }
        else { merged[it->first]+=it->second; }
    }
    return merged;
}

template<typename T>
std::pair<int,T> get_min_value_pair(const map<int,T>& histogram_data)
{
    using std::pair;
    int min_value(histogram_data.begin()->first);
    T corresponding_count(histogram_data.begin()->second);
    typedef typename std::map<int,T>::const_iterator ConstIterator;
    for (ConstIterator it=histogram_data.begin(); it!=histogram_data.end(); ++it) {
        if (it->first < min_value) {
            min_value = it->first; corresponding_count = it->second;
        }
    }
    pair<int,T> min_pair(min_value,corresponding_count);
    return min_pair;
}

template<typename T>
std::pair<int,T> get_max_value_pair(const map<int,T>& histogram_data)
{
    using std::pair;
    int max_value(histogram_data.begin()->first);
    T corresponding_count(histogram_data.begin()->second);
    typedef typename std::map<int,T>::const_iterator ConstIterator;
    for (ConstIterator it=histogram_data.begin(); it!=histogram_data.end(); ++it) {
        if (it->first > max_value) {
            max_value =it->first; corresponding_count = it->second;
        }
    }
    pair<int,T> max_pair(max_value,corresponding_count);
    return max_pair;
}


//======================================================================================================================
}//end namespace gain utils
}


#endif  /* _GAIN_UTILS_H */

