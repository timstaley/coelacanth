

#include "../gain_utils.h"
#include "coela_utility/src/string_utils.h"

#include "coela_utility/src/misc_math.h"

#include "coela_utility/src/microstats.h"

#include <boost/math/distributions/poisson.hpp>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>
#include <iostream>
#include <stdexcept>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnScan.h>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

using namespace std;
namespace coela {
namespace gain_utils {
//=====================================================================================================================

void append_histogram_data_from_float_bitmap_region(const CcdImage<float>& bmp,
        const PixelRange& rgn, HistogramContainer14bit& hist)
{
    assert(hist.minimum_has_been_set());
    assert(rgn.is_valid());
    for (int y=rgn.low.y; y<=rgn.high.y; ++y) {
        for (int x=rgn.low.x; x<=rgn.high.x; ++x) {
            hist.count_rounded_value(bmp.pix(x,y));
        }
    }//end pixel scan
    return;
}

void append_histogram_data_from_float_bitmap_region(const CcdImage<float>& bmp,
        const PixelRange& rgn, HistogramContainer10bit& hist)
{
    assert(hist.minimum_has_been_set());
    assert(rgn.is_valid());
    for (int y=rgn.low.y; y<=rgn.high.y; ++y) {
        for (int x=rgn.low.x; x<=rgn.high.x; ++x) {
            hist.count_rounded_value(bmp.pix(x,y));
        }
    }//end pixel scan
    return;
}




map<int,long> load_histogram_data_from_file(const string& filename)
{
    ifstream datafile(filename.c_str());
    if (!datafile) { throw runtime_error("Could not load histogram data from \""+filename+"\""); }
    int value; long count;
    map <int,long> histogram_data;
    while ((datafile>>value)) {
        datafile>>count;
        histogram_data[value] = count;
    }
    return histogram_data;
}




double histogram_differences_squared(const map<int,long>& data,
                                     const map<int,double>& model)
{
    assert(data.size()==model.size());
    double diff_sq_sum=0.0;
    double diff;
    for (map<int,long>::const_iterator it=data.begin(); it!=data.end(); ++it) {
        map<int,double>::const_iterator model_it = model.find(it->first);
        if (model_it==model.end()) { throw runtime_error("Gain utils: model/data mismatch during fitting"); }
        diff = it->second - model_it->second;
        diff_sq_sum += diff*diff;
    }
    return diff_sq_sum;
}

double poisson_noise_weighted_histogram_difference(const map<int,long>& data,
        const map<int,double>& model)
{
    assert(data.size()>=model.size());
    double diff_sq_sum=0.0;
    double diff;

//    pair<int,long> data_min_pair = get_min_value_pair(data);
//    pair<int,double> model_min_pair = get_min_value_pair(model);
//    cerr<<"Data / model min "<< data_min_pair.first<<" / "<< model_min_pair.first<<endl;

    for (map<int,double>::const_iterator model_it=model.begin(); model_it!=model.end();
            ++model_it) {
//        map<int,double>::const_iterator model_it = model.find(it->first);
//        if (model_it==model.end()) throw runtime_error("Gain utils: model/data mismatch during fitting");
        map<int,long>::const_iterator data_it = data.find(model_it->first);
        if (data_it==data.end()) throw
            runtime_error("Gain utils: model/data mismatch during fitting;"
                          "could not find datapoint for value " + string_utils::itoa(model_it->first));
        diff = (data_it->second - model_it->second);
        diff_sq_sum += diff*diff / max(data_it->second,
                                       1l);  //Poisson distribution for each histogram entry, so variance == probability
//        diff_sq_sum += pow(fabs(diff),1.5) / sqrt(model_it->second); //This method less sensitive to noise spikes at the end of the histogram section
    }
    return diff_sq_sum;
}

double log_space_histogram_difference(const map<int,long>& data,
                                      const map<int,double>& model)
{
    assert(data.size()>=model.size());
    double diff_sq_sum=0.0;
    double diff;
    for (map<int,double>::const_iterator model_it=model.begin(); model_it!=model.end();
            ++model_it) {
        map<int,long>::const_iterator data_it = data.find(model_it->first);
        if (data_it==data.end()) throw
            runtime_error("Gain utils: model/data mismatch during fitting;"
                          "could not find datapoint for value " + string_utils::itoa(model_it->first));
        diff = (log(data_it->second) - log(model_it->second));
        diff_sq_sum += diff*diff ;
    }
    return diff_sq_sum;
}

map<int, double> convolve_histogram_with_gaussian(const map<int, double>& input_hist,
        const double readout_sigma,
        const int kernel_half_width)
{
    map<int, double> readout_kernel;

    const double unity_area_gaussian_height = 1.0 / (sqrt(2.0*M_PI)*readout_sigma);

    const int kernel_min(-kernel_half_width), kernel_max(kernel_half_width);
    for (int value = kernel_min ; value <= kernel_max; ++value) {
        readout_kernel[value] = misc_math::gaussian_1d_function(value,
                                unity_area_gaussian_height,
                                readout_sigma);
    }

    const int input_min= get_min_value_pair(input_hist).first ;
    const int input_max= get_max_value_pair(input_hist).first ;

//    cerr<<"Kernel min / max" << kernel_min <<" / "<< kernel_max<<endl;
//    cerr<<"Input min / max" << input_min <<" / "<< input_max<<endl;


    map<int, double> convolved_hist;
    for (map<int,double>::const_iterator input_it=input_hist.begin();
            input_it!=input_hist.end(); ++input_it) {
        const int input_value= input_it->first;
        if (input_value + kernel_min >= input_min &&
                input_value + kernel_max <= input_max) {
            if (convolved_hist.find(input_value) == input_hist.end()) {
                convolved_hist[input_value ] = 0.0;
            }
            for (map<int,double>::const_iterator kernel_it=readout_kernel.begin();
                    kernel_it!=readout_kernel.end(); ++kernel_it) {

                int kernel_offset_input = kernel_it->first + input_value;
                convolved_hist[ input_value ] +=
                    input_hist.find(kernel_offset_input)->second * kernel_it->second;

            }
        }
//        else{
//            cerr<<"rejecting input val " << input_value<<endl;
//        }

    }
    return convolved_hist;
}

//==================================================================================

double GaussianHistogramFitFCN::operator()(const std::vector<double>& gauss_pars) const
{
    assert(gauss_pars.size()==3);
    map<int,double> model=get_model_histogram(gauss_pars);
//            return histogram_differences_squared(data, model);
//            cerr<<"Pars: "<< gauss_pars[0]<<" " <<gauss_pars[1]<<" "<<gauss_pars[2]<<endl;
    return poisson_noise_weighted_histogram_difference(data, model);

}

map<int,double> GaussianHistogramFitFCN::get_model_histogram(
    const std::vector<double>& gauss_pars) const
{
    map<int,double>model;
    for (map<int,long>::const_iterator it=data.begin(); it!=data.end(); ++it) {
        model[it->first] = misc_math::gaussian_1d_function((double)it->first - gauss_pars[0],
                           gauss_pars[1],
                           gauss_pars[2]);  //FIXME -- should perhaps be measuring distance from centre of bin? (so +0.5 to it->first). but not sure.
    }
    return model;
}

//==================================================================================
map<int,double> ExponentialHistogramFitFCN::get_model_histogram(
    const std::vector<double>& exp_pars) const
{
    assert(exp_pars.size()==2);
    map<int,double>model;
    const double& amplitude = exp_pars[0];
    const double& gain = exp_pars[1];
    for (map<int,long>::const_iterator it=data.begin(); it!=data.end(); ++it) {
        model[it->first] = amplitude*exp(-1.0*it->first/gain);
    }
    return model;
}

double  ExponentialHistogramFitFCN::operator()(const std::vector<double>& exp_pars)
const
{
    assert(exp_pars.size()==2);
    map<int,double> model=get_model_histogram(exp_pars);
//    return log_space_histogram_difference(data, model);
    return poisson_noise_weighted_histogram_difference(data, model) / (model.size() - 1);

}
//==================================================================================
//map<int,double> full_EMCCD_low_light_level_histogram_model::get_model_histogram(
//        const GainData& inf){
//
//}
//

//==================================================================================
FullEmccdHistogramFitFCN::FullEmccdHistogramFitFCN(const map<int,long>& histogram)
    :
    error_def(1.0), data_(histogram)
{
    fit_rgn_min_ = get_min_value_pair(data_).first;
    fit_rgn_max_ = get_max_value_pair(data_).first;
}

FullEmccdHistogramFitFCN::FullEmccdHistogramFitFCN(const map<int,long>& histogram,
        const int fit_region_min, const int fit_region_max):
    error_def(1.0), data_(histogram),
    fit_rgn_min_(fit_region_min), fit_rgn_max_(fit_region_max)
{}

double FullEmccdHistogramFitFCN::operator()(const std::vector<double>& pars) const
{
    assert(pars.size()==npars_);

//    cout<<"Params: ";
//    for (size_t i=0; i!=pars.size();++i){
//        cout<<pars[i]<<" ";
//    }
//
    map<int,double> model = get_model_histogram(pars);
    double val = poisson_noise_weighted_histogram_difference(data_, model);
//
//    cout<<" : "<<val;
//    cout<<endl;
//
    double normalised = val / (model.size() - (npars_ - 1));

    return normalised;
}
vector<double> FullEmccdHistogramFitFCN::generate_pars_vec(const GainData& inf)
{
    vector<double> pars;
    pars.push_back(inf.bias_pedestal);
    pars.push_back(inf.N_dark_pix);
    pars.push_back(inf.readout_sigma);
    pars.push_back(inf.N_light_pix);
    pars.push_back(inf.gain);
    pars.push_back(inf.N_serial_CIC_pix);
//    pars.push_back(inf.serial_CIC_gain);
    return pars;
}


map<int,double> FullEmccdHistogramFitFCN::get_CICIR_model_histogram(
    const double gain,
    const double N_CICIR_pix,
    const int N_EM_stages,
    const int max_val_to_generate)
{

    int min_N_stages=1;
    double prob_per_stage = pow(gain, 1.0/N_EM_stages);


    map<int,double> gain_for_N_stages;
    for (int stages_passed = 1; stages_passed<=N_EM_stages; ++stages_passed) {
        gain_for_N_stages[stages_passed] = pow(prob_per_stage, stages_passed);
    }

//    cerr<<"gain at 604:"<< gain_for_N_stages[604]<<endl;
//    const double frequency_normalisation =  N_CICIR_pix ;
    const double frequency_normalisation = N_CICIR_pix / (N_EM_stages - min_N_stages);

    map<int,double> CICIR_histogram;
    for (int register_output_val=1; register_output_val<=max_val_to_generate;
            ++register_output_val) {
        double output_frequency=0.0;
        for (int CICIR_stages_passed = min_N_stages; CICIR_stages_passed<=N_EM_stages;
                ++CICIR_stages_passed) {
            const double & g = gain_for_N_stages[CICIR_stages_passed];

            output_frequency +=  exp((0.5-register_output_val)/g) / g;
        }
        CICIR_histogram[register_output_val]= frequency_normalisation * output_frequency;
    }
    return CICIR_histogram;
}

map<int,double> FullEmccdHistogramFitFCN::get_model_histogram(
    const std::vector<double>& pars) const
{
    assert(pars.size()==npars_);

    const double& bias_pedestal = pars[0];
    const double& dark_pix = pars[1];
    const double& RO_sigma = pars[2];
    const double& N_photon_events = pars[3];
    const double& phot_gain = pars[4];
    const double& N_serial_CIC_pix = pars[5];
//    const double& serial_CIC_gain = pars[6];
//    const double serial_CIC_gain = sqrt(phot_gain);

    int convolution_half_width=60;

//   int data_start = get_min_value_pair(data).first;
    int data_end = get_max_value_pair(data_).first;
    int tail_range = data_end - bias_pedestal;

    double total_pix_calc = N_photon_events + N_serial_CIC_pix + dark_pix;
    double photon_freq_est = N_photon_events / total_pix_calc;
//   cerr<<"Photon freq est: " << photon_freq_est<<", calc total: " << total_pix_calc<<endl;

    //See Basden et al, 2003 for equations


//   double serial_CIC_amplitude =
//    N_serial_CIC_pix * exp(0.5/serial_CIC_gain)/ serial_CIC_gain;

    double p1_amplitude, p2_amplitude;
    if (photon_freq_est>0) {
        boost::math::poisson_distribution<> pois_dist(photon_freq_est);
        double p1_prob = boost::math::pdf(pois_dist, 1);
        double p2_prob = boost::math::pdf(pois_dist, 2);

        //Note exp( 0.5 / phot_gain) term corrects ~1% normalisation error;
        // since integral from 1 to +inf of the PDF otherwise  =  exp(-1/g)
        p1_amplitude = (p1_prob / (p1_prob +p2_prob)) *
                       N_photon_events * exp(0.5/phot_gain) /  phot_gain;

        //Note normalisation correction not reqd here, more complex to calc. but also much smaller effect.
        p2_amplitude = (p2_prob / (p1_prob +p2_prob)) *
                       N_photon_events / (phot_gain*phot_gain);
    } else {
        p1_amplitude=0.0;
        p2_amplitude=0.0;
    }

    map<int,double> register_output_model;
    int fit_rgn_extension_below_bias_pedestal = bias_pedestal - fit_rgn_min_;
    if (fit_rgn_extension_below_bias_pedestal > 0) {
        for (int i= -(convolution_half_width + fit_rgn_extension_below_bias_pedestal); i<0; ++i) {
            register_output_model[i]=0;
        }
    } else {
        for (int i= -convolution_half_width; i<0; ++i) {
            register_output_model[i]=0;
        }
    }

    register_output_model[0]  = dark_pix;

    map<int,double> CICIR_register_component =
        FullEmccdHistogramFitFCN::get_CICIR_model_histogram(
            phot_gain, N_serial_CIC_pix,
            FullEmccdHistogramFitFCN::N_EM_serial_register_stages,
            tail_range + convolution_half_width
        );

    for (int register_output_val=1; register_output_val<=tail_range + convolution_half_width;
            ++register_output_val) {
        double photon_exponential =  exp(-1.0* register_output_val /phot_gain);
        //Contribution from 1 and 2 photo-electron events
        register_output_model[register_output_val] =
            p1_amplitude*photon_exponential;
        register_output_model[register_output_val] +=
            p2_amplitude*(register_output_val - 1)*photon_exponential;

        //Contribution from CICIR events
//        register_output_model[register_output_val] +=
//                serial_CIC_amplitude*exp(  -1.0* register_output_val /serial_CIC_gain );
        register_output_model[register_output_val] +=
            CICIR_register_component[register_output_val];
    }

    map<int,double> RO_model =
        convolve_histogram_with_gaussian(register_output_model, RO_sigma,
                                         convolution_half_width);

    map<int,double> biased_RO_model;
    for (map<int,double>::const_iterator it=RO_model.begin(); it!=RO_model.end(); ++it) {
        int output_value = it->first + bias_pedestal;
        if (output_value>=fit_rgn_min_ && output_value <= fit_rgn_max_) {
            biased_RO_model[output_value] = it->second;
        }
    }

    return biased_RO_model;
}

//==================================================================================

//double calculate_thresholded_SNR(const double readout_noise_in_ADU,
//                                const double mean_photon_gain_in_ADU,
//                                const double CIC_event_rate_per_pixel_readout,
//                                const double light_level,
//                                const double threshold_in_photo_electrons,
//                                const bool output_to_screen){
//    map<int,long> histogram_min_max;
//    histogram_min_max[-readout_noise_in_ADU*5]=0;
//    histogram_min_max[10000]=0;
//
//    FullEmccdHistogramFitFCN full_model(histogram_min_max);
//    GainData full_pars;
//    full_pars.bias_pedestal = 0.0;
//    full_pars.readout_sigma = readout_noise_in_ADU;
//    full_pars.gain = mean_photon_gain_in_ADU;
//    full_pars.N_light_pix = 0.0;
//    full_pars.N_serial_CIC_pix = 0.0;
//    full_pars.N_dark_pix = 0.0;
////                - mean_photon_gain_in_ADU*CIC_event_rate_per_pixel_readout);
//
//    GainData photon_hist_pars(full_pars);
//    photon_hist_pars.N_light_pix=1.0;
//    map<int,double> light_model = full_model.get_model_histogram(
//        full_model.generate_pars_vec(photon_hist_pars));
//
//    double all_counts = gain_utils::sum_counts_in_histogram(light_model);
//    double thresh_counts = gain_utils::sum_counts_in_histogram_above_threshold(
//            light_model, threshold_in_photo_electrons*mean_photon_gain_in_ADU);
//
//    double pass_rate_for_photon_events = thresh_counts / all_counts;
//
//    if (output_to_screen){
////        gain_utils::write_histogram_data_to_file(light_model, "photons_only.txt");
//        cerr<<"photons hist sum:" << all_counts << endl;
//        cerr<<"Filtered photon fraction: "<< pass_rate_for_photon_events<<endl;
//    }
//
////    GainData CICIR_hist_pars(full_pars);
////    CICIR_hist_pars.N_serial_CIC_pix=1.0;
////    map<int,double> CICIR_model = full_model.get_model_histogram(
////        full_model.generate_pars_vec(CICIR_hist_pars));
////
////    all_counts = gain_utils::sum_counts_in_histogram(CICIR_model);
////    thresh_counts = gain_utils::sum_counts_in_histogram_above_threshold(
////            CICIR_model, threshold_in_photo_electrons*mean_photon_gain_in_ADU);
////
////    double CICIR_fraction = thresh_counts / all_counts;
////
////    if (output_to_screen){
////        gain_utils::write_histogram_data_to_file(CICIR_model, "CICIR_only.txt");
////        cerr<<"CICIR hist sum:" << all_counts << endl;
////        cerr<<"Filtered CICIR fraction: "<< CICIR_fraction<<endl;
////    }
////
////    GainData RO_hist_pars(full_pars);
////    RO_hist_pars.N_dark_pix=1.0;
////    map<int,double> RO_model = full_model.get_model_histogram(
////        full_model.generate_pars_vec(RO_hist_pars));
////
////    all_counts = gain_utils::sum_counts_in_histogram(RO_model);
////    thresh_counts = gain_utils::sum_counts_in_histogram_above_threshold(
////            RO_model, threshold_in_photo_electrons*mean_photon_gain_in_ADU);
////
////    double RO_fraction = thresh_counts / all_counts;
////
////    if (output_to_screen){
////        gain_utils::write_histogram_data_to_file(RO_model, "RO_only.txt");
////        cerr<<"RO hist sum:" << all_counts << endl;
////        cerr<<"Filtered RO fraction: "<< RO_fraction<<endl;
////    }
//
//    full_pars.N_light_pix = light_level;
//    full_pars.N_serial_CIC_pix = CIC_event_rate_per_pixel_readout;
//    full_pars.N_dark_pix = 1.0 - (light_level + CIC_event_rate_per_pixel_readout);
//
//    map<int,double> full_model_hist = full_model.get_model_histogram(
//        full_model.generate_pars_vec(full_pars));
//
//    all_counts = gain_utils::sum_counts_in_histogram(full_model_hist);
//    thresh_counts = gain_utils::sum_counts_in_histogram_above_threshold(
//            full_model_hist, threshold_in_photo_electrons*mean_photon_gain_in_ADU);
//
//    double pass_rate_for_all_events = thresh_counts / all_counts;
//    if (output_to_screen){
//        cerr<<"Noise hist sum:" << all_counts << endl;
//        cerr<<"Filtered noise fraction: "<< pass_rate_for_all_events <<endl;
//    }
//
//    double SNR = pass_rate_for_photon_events*light_level /
//        sqrt( pass_rate_for_all_events );
//
//    if (output_to_screen){
//        cerr<<"SNR = " << SNR<<endl;
//    }
//    return SNR;
//}

//==================================================================================

ThresholdedSnrCalculator::ThresholdedSnrCalculator(
    const double RO,
    const double gain,
    const double light_level,
    const double CICIR_rate
):gain_(gain), light_level_(light_level)
{
    if ((light_level) > 0.5) {
        throw runtime_error("ThresholdedSnrCalculator() ---"
                            "Light level too high");
    }

    map<int,long> histogram_min_max;
    histogram_min_max[-RO*5]=0;
    histogram_min_max[10000]=0;

    FullEmccdHistogramFitFCN full_model(histogram_min_max);
    GainData full_pars;
    full_pars.bias_pedestal = 0.0;
    full_pars.readout_sigma = RO;
    full_pars.gain = gain;
    full_pars.N_light_pix = 0.0;
    full_pars.N_serial_CIC_pix = 0.0;
    full_pars.N_dark_pix = 0.0;
//                - mean_photon_gain_in_ADU*CIC_event_rate_per_pixel_readout);

    GainData photon_hist_pars(full_pars);
    photon_hist_pars.N_light_pix=1.0;

    photon_distribution_ = full_model.get_model_histogram(
                               full_model.generate_pars_vec(photon_hist_pars));

    full_pars.N_light_pix = light_level;
    full_pars.N_serial_CIC_pix = CICIR_rate;
    full_pars.N_dark_pix = 1.0 - (light_level + CICIR_rate);

    full_distribution_ = full_model.get_model_histogram(
                             full_model.generate_pars_vec(full_pars));

//    gain_utils::write_histogram_data_to_file(photon_distribution_, "photons.txt");
//    gain_utils::write_histogram_data_to_file(full_distribution_, "full_dist.txt");

}

double ThresholdedSnrCalculator::calc_SNR_at_threshold(
    const double threshold_in_photo_electrons) const
{

    double all_counts = gain_utils::sum_counts_in_histogram(photon_distribution_);
    double thresh_counts = gain_utils::sum_counts_in_histogram_above_threshold(
                               photon_distribution_, threshold_in_photo_electrons*gain_);
    double fraction_photon_events_passed = thresh_counts / all_counts;

    all_counts = gain_utils::sum_counts_in_histogram(full_distribution_);
    thresh_counts = gain_utils::sum_counts_in_histogram_above_threshold(
                        full_distribution_, threshold_in_photo_electrons*gain_);

    double fraction_of_all_events_passed = thresh_counts / all_counts;
//    cerr<<"Gain"<< gain_<<endl;
//    cerr<<"photon sum "<<all_counts<<endl;
//    cerr<<"Counted over thresh:"<<thresh_counts<<endl;
//    cerr<<"photon frac"<<pass_rate_for_photon_events<<endl;
//    cerr<<"Total sum"<<all_counts<<endl;
//    cerr<<"Total frac"<<pass_rate_for_all_events<<endl;


    double SNR = fraction_photon_events_passed*light_level_ /
                 sqrt(exp(fraction_of_all_events_passed) - 1.0);
    return SNR;
}

double ThresholdedSnrCalculator::find_optimum_threshold(
    const double search_min,
    const double search_max,
    const double stepsize
) const
{
    double max_SNR = calc_SNR_at_threshold(search_min);
    double best_thresh = search_min;
    for (double t = search_min; t<=search_max; t+=stepsize) {
        double s = calc_SNR_at_threshold(t);
        if (s > max_SNR) {
            best_thresh = t;
            max_SNR = s;
        }
    }
    return best_thresh;
}

double calculate_ideal_photon_counter_SNR(const double light_level)
{
    return light_level /
           sqrt(exp(light_level) - 1.0);
}

double calculate_linear_mode_SNR(const double readout_noise_in_ADU,
                                 const double mean_photon_gain_in_ADU,
                                 const double light_level,
                                 const double CICIR_event_rate_per_pixel_readout
                                )
{

    ///Calculate CICIR_signal relative to a photo-electron signal:
    double CICIR_signal =  CICIR_event_rate_per_pixel_readout /
                           log(mean_photon_gain_in_ADU);
    double normed_readout = readout_noise_in_ADU / mean_photon_gain_in_ADU;
    return light_level /
           sqrt(2.0*(light_level + CICIR_signal)  + normed_readout*normed_readout);
}

//==================================================================================
pair<int,long> get_mode(const map<int,long>& histogram_data)
{
    int peak_value(histogram_data.begin()->first);
    long peak_count(histogram_data.begin()->second);
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->second > peak_count) {
            peak_value  =it->first; peak_count = it->second;
        }
    }
    pair<int,long> peak_pair(peak_value, peak_count);
    return peak_pair;
}



float get_gaussian_HWHM(const map<int,long>& histogram_data)
{
    pair<int,long> peak = get_mode(histogram_data);
    float best_distance_to_HM=peak.second/2.0;
    int best_value = peak.first;
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->first < peak.first) {
            float current_dist_to_HM = fabs(it->second - peak.second/2.0);
            if (current_dist_to_HM < best_distance_to_HM) {
                best_distance_to_HM=current_dist_to_HM;
                best_value =it->first;
            }
        }
    }
//            cout <<"Half peak value: " <<best_value<<endl;
//            cout <<"HWHM "<<fabs(best_value -peak.first)<<endl;
    return fabs(best_value -peak.first);
}


map<int, long> get_gaussian_fitting_section(const map<int,long>& histogram_data)
{
    map<int,long> gaussian_section;
    pair<int,long> peak = get_mode(histogram_data);
    long count_threshold = max((long)(peak.second*0.005), 300l);
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->first < (peak.first+10) && it->second > count_threshold) { gaussian_section[it->first]=it->second; }
    }
    return gaussian_section;
}



map<int, long> get_tail_fitting_section(const map<int,long>& histogram_data,
                                        const double g_sigma_est, const bool output_to_screen)
{
    pair<int,long>peak =get_mode(histogram_data);
    long count_threshold=max((long)25,
                             (long)floor(peak.second/6e6)); //FIXME - this may need adjusting
    if (output_to_screen) { cout <<"\nTail count cutoff: " << count_threshold<<endl; }
    int tail_section_end_value = peak.first;
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->first > peak.first && it->second >count_threshold) {
            tail_section_end_value= max(tail_section_end_value, it->first);
        }
    }
    tail_section_end_value-=10; //jump back, otherwise may get extended by a noise spike extending the range to... a noise spike!
    if (output_to_screen) { cout <<"Resulting in cutoff value: "<<tail_section_end_value<<endl; }
    map<int,long> tail_section;
    int curve_range = tail_section_end_value - peak.first;
    int peak_offset = max(curve_range*0.25, g_sigma_est*8);
    int tail_start = (peak.first + peak_offset);
    if (output_to_screen) { cout <<"Tail_start at: "<<tail_start<<endl; }
    if (tail_start >= tail_section_end_value -5) { throw runtime_error("Not enough data to fit tail of histogram"); }
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->first >= tail_start && it->first <=tail_section_end_value) { tail_section[it->first] = it->second; }
    }
    return tail_section;
}

map<int,long> get_serial_CIC_model_fitting_section(const map<int,long>& histogram_data,
        const double bias_pedestal,
        const double readout_sigma)
{
    int low_cutoff = bias_pedestal + 4.0 * readout_sigma;
    int high_cutoff = bias_pedestal + 12*readout_sigma;
    map<int, long> fitting_region;
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->first > low_cutoff && it->first <=high_cutoff) { fitting_region[it->first] = it->second; }
    }
    return fitting_region;
}

map<int,long> get_full_model_fitting_section(const map<int,long>& histogram_data,
        const double bias_pedestal,
        const double readout_sigma)
{
//    pair<int,long>peak =get_mode(histogram_data);
//    long count_threshold=max((long)50, (long)floor(peak.second/1e6));

    map<int, long> tail_region = get_tail_fitting_section(histogram_data,
                                 readout_sigma);
    int tail_end = get_max_value_pair(tail_region).first;

    int low_cutoff = bias_pedestal - 0.5 * readout_sigma;
//    cerr<<"low cuttoff, tail end: " <<low_cutoff <<" , "<< tail_end <<endl;

    map<int, long> fitting_region;
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if (it->first > low_cutoff && it->first <=tail_end) { fitting_region[it->first] = it->second; }
    }
    return fitting_region;
}

map<int, long> get_points_at_values_greater_than(const map<int,long>& histogram_data,
        int zero_point)
{
    map<int, long> positive_section;
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if ((it->first > zero_point)) { positive_section[it->first]=it->second; }
    }
    return positive_section;
}

map<int, long> threshold_above_count(const map<int,long>& histogram_data,
                                     int count_threshold)
{
    map<int, long> threshed;
    for (map<int,long>::const_iterator it=histogram_data.begin(); it!=histogram_data.end();
            ++it) {
        if ((it->second > count_threshold)) { threshed[it->first]=it->second; }
    }
    return threshed;
}


float estimate_log_slope(const map<int,long>& tail_section)
{
    pair<int,long> tail_section_peak(get_mode(tail_section)),
         tail_end(get_max_value_pair(tail_section));
    return (log(tail_end.second) - log(tail_section_peak.second))/
           (tail_end.first - tail_section_peak.first);
}

//Returns gauss_par vector: peak value (zero value), peak count, sigma (readout sigma);
vector<double> fit_readout_gaussian(const map<int,long>& histogram_data,
                                    const bool output_to_screen)
{
    pair<int,long> peak = get_mode(histogram_data);
    if (output_to_screen) { cout <<"Data peak count at " <<peak.first<<" : "<<peak.second<<endl; }
    map<int, long> gaussian_section = get_gaussian_fitting_section(histogram_data);
    gaussian_section.erase(
        0); //remove weird spike at 0 value. fixme! - track down the origin of this.

    GaussianHistogramFitFCN hump_fitter(gaussian_section);
    ROOT::Minuit2::MnUserParameters gauss_pars;
    gauss_pars.Add("HistPeakValue", peak.first, 2);
    gauss_pars.SetLimits("HistPeakValue", -100.0, 20000.0);
    gauss_pars.Add("HistPeakCount", peak.second, peak.second*0.1);
//            gauss_pars.SetLimits("HistPeakCount", peak.first*0.5, peak.first*2.0);
    gauss_pars.SetLowerLimit("HistPeakCount",5);
    float gauss_sigma_est = get_gaussian_HWHM(histogram_data) * 2 /
                            misc_math::gaussian_fwhm_to_sigma_ratio;
    if (output_to_screen) { cout <<"Sigma init est. " << gauss_sigma_est <<endl; }

    gauss_pars.Add("HistGaussSigma", gauss_sigma_est, gauss_sigma_est*0.3);
    gauss_pars.SetLowerLimit("HistGaussSigma",0.1);

//            gauss_pars.SetLimits("HistGaussSigma",gauss_sigma_est*0.5, gauss_sigma_est*2.0);
    ROOT::Minuit2::MnMinimize min_finder(hump_fitter, gauss_pars);
    ROOT::Minuit2::FunctionMinimum min=min_finder();
    if (output_to_screen) { cout <<min<<endl; }
    vector<double> end_pars;
    end_pars=min.UserParameters().Params();
    return end_pars;
}

//Returns amplitude, gain
vector<double> fit_exponential_tail(const map<int,long>& histogram_data,
                                    const double gaussian_sigma_estimate,
                                    const bool output_to_screen)
{
    if (output_to_screen) { cout <<"Fitting tail"<<endl; }
    pair<int,long>peak =get_mode(histogram_data);
    map<int,long> tail_section = get_tail_fitting_section(histogram_data,
                                 gaussian_sigma_estimate ,output_to_screen);

    ExponentialHistogramFitFCN tail_fitter(tail_section);
    double log_slope_estimate = estimate_log_slope(tail_section);
    if (output_to_screen) { cout <<"Initial slope estimate: " <<log_slope_estimate<<endl; }
    float gain_init_est = -1.0/log_slope_estimate;
    if (output_to_screen) { cout <<" - gain therefore est. at "<<gain_init_est<<endl; }
    ROOT::Minuit2::MnUserParameters tail_par;
//            tail_par.SetPrecision(1e-7);

    double amplitude_estimate=peak.second*0.1 / exp(
                                  -1.0*peak.first/gain_init_est);     //assumes 1/10th of all pixels contain a photon event
    if (output_to_screen) { cout <<"Initial amp estimate: " <<amplitude_estimate<< endl; }
    tail_par.Add("Amp", amplitude_estimate, amplitude_estimate*0.5);
    tail_par.SetLowerLimit("Amp", 0);
    tail_par.Add("Gain", gain_init_est, gain_init_est*0.1);
    tail_par.SetLimits("Gain", 1.0,1000.0);
    tail_par.Fix("Gain"); //Gain estimate likely much more accurate than amplitude estimate, so fix this for first run
    ROOT::Minuit2::MnMinimize amp_fitter(tail_fitter, tail_par);
    ROOT::Minuit2::FunctionMinimum min = amp_fitter();
    tail_par.SetValue("Amp",min.UserParameters().Value("Amp"));
    tail_par.Release("Gain");

    ROOT::Minuit2::MnMinimize amp_gain_fitter(tail_fitter, tail_par);
    min = amp_gain_fitter();
    if (output_to_screen) { cout <<min<<endl; }
    ;
    vector<double> fitted_tail_pars = min.UserParameters().Params();
    return fitted_tail_pars;
}

vector<double> fit_exponential_tail_fixed_gain(const map<int,long>& histogram_data,
        const double gaussian_sigma_estimate, const double gain_estimate,
        const bool output_to_screen)
{
    if (output_to_screen) { cout <<"Fitting tail"<<endl; }
    pair<int,long>peak =get_mode(histogram_data);
    map<int,long> tail_section = get_tail_fitting_section(histogram_data,
                                 gaussian_sigma_estimate ,output_to_screen);

    ExponentialHistogramFitFCN tail_fitter(tail_section);

    double gain_init_est=gain_estimate;
    double amplitude_estimate=peak.second*0.1 / exp(
                                  -1.0*peak.first/gain_init_est);     //assumes 1/10th of all pixels contain a photon event
    ROOT::Minuit2::MnUserParameters tail_par;
    if (output_to_screen) { cout <<"Initial amp estimate: " <<amplitude_estimate<< ", est err: "<< peak.second*0.1<< endl; }
    tail_par.Add("Amp", amplitude_estimate, peak.second*0.1);
    tail_par.SetLowerLimit("Amp", 0);
    tail_par.Add("Gain", gain_init_est, gain_init_est*0.05);
    tail_par.SetLowerLimit("Gain", 1.0);
    tail_par.Fix("Gain");
    ROOT::Minuit2::MnMinimize amp_finder(tail_fitter, tail_par);
    ROOT::Minuit2::FunctionMinimum min = amp_finder();
    if (output_to_screen) { cout <<min<<endl; }
//
//            if (output_to_screen)cout <<"Refitting gain only"<<endl;
//            tail_par.Fix("Amp");
//            ROOT::Minuit2::MnMinimize gain_min_finder(tail_fitter, tail_par, 2);
//            min = gain_min_finder();
//            if (output_to_screen)cout <<min<<endl;
    vector<double> fitted_tail_pars = min.UserParameters().Params();
    return fitted_tail_pars;
}



GainData fit_full_CCD_model(const map<int,long>& histogram_data,
                             const GainData& approx_fit_init,
                             const bool output_to_screen
//        , const bool sqrt_serial_CIC_gain_mode
                            )
{
    //Hardcoded initial estimates (fixme?)
    GainData approx_fit(approx_fit_init);

    approx_fit.N_serial_CIC_pix= approx_fit.actual_number_pix_events_recorded -
                                 (approx_fit.N_dark_pix + approx_fit.N_light_pix);
//    approx_fit.serial_CIC_gain = sqrt(approx_fit.gain);



//    gain_utils::write_histogram_data_to_file(full_fitting_region, "full_fitting_region.txt");

    ROOT::Minuit2::MnUserParameters mnpars;
    const string key_BiasPedestal = "BiasPedestal";
    const string key_NDarkPix = "NDarkPix";
    const string key_ROSigma = "ROSigma";
    const string key_NPhotPix = "NPhotPix";
    const string key_PhotGain = "PhotGain";
    const string key_NSerialCICPix = "NSerCICPix";
//    const string key_SerialCICGain = "SerCICGain";

    //NB! Ensure ordering matches definition for FullEmccdHistogramFitFCN::operator() !!!

    mnpars.Add(key_BiasPedestal, approx_fit.bias_pedestal, 2);
    mnpars.SetLimits(key_BiasPedestal, approx_fit.bias_pedestal-2,
                     approx_fit.bias_pedestal+2);

    mnpars.Add(key_NDarkPix, approx_fit.N_dark_pix, approx_fit.N_dark_pix*0.1);
    mnpars.SetLowerLimit(key_NDarkPix,5);
//    if (output_to_screen) cout <<"Sigma init est. " << gauss_sigma_est <<endl;

    mnpars.Add(key_ROSigma, approx_fit.readout_sigma, approx_fit.readout_sigma*0.1);
    mnpars.SetLimits(key_ROSigma, approx_fit.readout_sigma*0.8, approx_fit.readout_sigma*1.2);

    mnpars.Add(key_NPhotPix, approx_fit.N_light_pix, approx_fit.N_light_pix*0.1);
    mnpars.SetLowerLimit(key_NPhotPix, 0);

    mnpars.Add(key_PhotGain, approx_fit.gain, approx_fit.gain*0.05);
    mnpars.SetLimits(key_PhotGain, 1.0, 1000);

    mnpars.Add(key_NSerialCICPix, approx_fit.N_serial_CIC_pix,
               approx_fit.actual_number_pix_events_recorded*0.2);
    mnpars.SetLowerLimit(key_NSerialCICPix, 0.0);

//    mnpars.Add(key_SerialCICGain, approx_fit.serial_CIC_gain,
//            approx_fit.serial_CIC_gain);
//    mnpars.SetLowerLimit(key_SerialCICGain, 1.0);

    //--------------------------------------------------------------------------------
    //Fit just the serial CIC

    mnpars.Fix(key_BiasPedestal);
    mnpars.Fix(key_NDarkPix);
    mnpars.Fix(key_ROSigma);
    mnpars.Fix(key_NPhotPix);
    mnpars.Fix(key_PhotGain);

//    if (sqrt_serial_CIC_gain_mode){
//        mnpars.SetValue(key_SerialCICGain, sqrt(mnpars.Value(key_PhotGain)));
//        mnpars.Fix(key_SerialCICGain);
//    }
//    mnpars.Fix(key_NSerialCICPix); //Just leave this free to start:

    map<int, long> serial_CIC_region =
        get_serial_CIC_model_fitting_section(histogram_data,
                approx_fit.bias_pedestal, approx_fit.readout_sigma);

    FullEmccdHistogramFitFCN CIC_fit_fcn(serial_CIC_region);

    if (output_to_screen) cout<<"---------------------------------\n"
                                  << "Init params:" << mnpars<<endl;
//    gain_utils::write_histogram_data_to_file(serial_CIC_region.get_model_histogram(mnpars.Params()), "initial_model.txt" );
//    ROOT::Minuit2::MnScan serial_CIC_fitter(CIC_fit_fcn, mnpars);
//    serial_CIC_fitter.Scan(5);

    ROOT::Minuit2::MnMinimize serial_CIC_fitter(CIC_fit_fcn, mnpars);

    ROOT::Minuit2::FunctionMinimum min0 = serial_CIC_fitter();
//    if (output_to_screen)cout <<min<<endl;

    ROOT::Minuit2::MnUserParameters fitted_pars = min0.UserParameters();
    if (output_to_screen) { cout<<"Params after fit 0:" << fitted_pars << endl; }

    //--------------------------------------------------------------------------------
    map<int, long> full_fitting_region =
        get_full_model_fitting_section(histogram_data,
                                       approx_fit.bias_pedestal, approx_fit.readout_sigma);
//
    FullEmccdHistogramFitFCN full_fit_fcn(full_fitting_region);

    mnpars = fitted_pars;
    mnpars.Release(key_BiasPedestal);
    mnpars.Release(key_NDarkPix);
    mnpars.Release(key_ROSigma);
    mnpars.Release(key_NPhotPix);
//    if (sqrt_serial_CIC_gain_mode == false){
    mnpars.Release(key_PhotGain);
//    }
//    mnpars.Release(key_NSerialCICPix); //Already free
//    mnpars.Release(key_SerialCICGain); //Already free
    mnpars.SetError(key_NDarkPix, mnpars.Value(key_NDarkPix)*0.1);
    mnpars.SetError(key_NPhotPix, mnpars.Value(key_NPhotPix)*0.1);
    mnpars.SetError(key_PhotGain, mnpars.Value(key_PhotGain)*0.05);
    mnpars.SetError(key_NSerialCICPix, mnpars.Value(key_NSerialCICPix)*0.1);
//    mnpars.SetError(key_SerialCICGain, mnpars.Value(key_SerialCICGain)*0.1);


    if (output_to_screen) cout<<"---------------------------------\n"
                                  << "Init params:" << mnpars<<endl;
    ROOT::Minuit2::MnMinimize free_fitter(full_fit_fcn, mnpars);
    ROOT::Minuit2::FunctionMinimum min3 = free_fitter();
    if (output_to_screen) { cout <<min3<<endl; }
    fitted_pars = min3.UserParameters();
    if (output_to_screen) { cout<<"Params after fit 3:" << fitted_pars << endl; }
    //--------------------------------------------------------------------------------
    GainData fitted_info(approx_fit);
    if (min3.IsValid()) {
        fitted_info.bias_pedestal = fitted_pars.Value(key_BiasPedestal);
        fitted_info.readout_sigma = fitted_pars.Value(key_ROSigma);
        fitted_info.N_dark_pix = fitted_pars.Value(key_NDarkPix);
        fitted_info.gain = fitted_pars.Value(key_PhotGain);
        fitted_info.N_light_pix = fitted_pars.Value(key_NPhotPix);
        fitted_info.N_serial_CIC_pix = fitted_pars.Value(key_NSerialCICPix);
//        fitted_info.serial_CIC_gain = fitted_pars.Value(key_SerialCICGain);

        const double total_calculated_hist_area = fitted_info.N_dark_pix +
                fitted_info.N_serial_CIC_pix  + fitted_info.N_light_pix;

        fitted_info.photon_event_freq = fitted_info.N_light_pix /
                                        total_calculated_hist_area;

        fitted_info.serial_CIC_rate =  fitted_info.N_serial_CIC_pix /
                                       (fitted_info.N_serial_CIC_pix +
                                        fitted_info.N_dark_pix); //NB light pixels do not show up in serial CIC rate count!
        if (output_to_screen) {
            cerr<<"Calculated hist area / actual counts " << total_calculated_hist_area
                <<" / "<<fitted_info.actual_number_pix_events_recorded<<endl;

            cerr<<"Calculated N_dark pix / light pix / n_CIC_pix"<< fitted_info.N_dark_pix
                << " / "<<fitted_info.N_light_pix
                << " / "<<fitted_info.N_serial_CIC_pix
                <<endl;
        }
    }


    return fitted_info;
}


GainData fit_CCD_histogram(const map<int,long>& histogram_bins,
                            const bool perform_advanced_fit,
                            const bool fixed_gain_mode, const double gain_value,
                            const bool output_to_screen
//        , const bool sqrt_serial_CIC_gain_mode
                           )
{
    if (fixed_gain_mode &&  gain_value==0.0) { throw logic_error("Must supply a gain value for fitting in fixed gain mode"); }
    vector<double> gauss_pars = gain_utils::fit_readout_gaussian(histogram_bins);
    vector<double> exp_pars;
    if (fixed_gain_mode) { exp_pars = gain_utils::fit_exponential_tail_fixed_gain(histogram_bins, gauss_pars.back(), gain_value, true); }
    else { exp_pars= gain_utils::fit_exponential_tail(histogram_bins, gauss_pars.back(), output_to_screen && !perform_advanced_fit); }

    GainData approx_results;
    approx_results.bias_pedestal = gauss_pars[0];
    double peak_gaussian_count = gauss_pars[1];
    approx_results.readout_sigma = gauss_pars[2];
    approx_results.gain = exp_pars[1];

    approx_results.N_dark_pix= peak_gaussian_count*approx_results.readout_sigma *sqrt(
                                   2*M_PI); //area under a gaussian

    double area_under_exponential_curve = approx_results.gain *
                                          exp_pars.front()* exp(
                                                  -1.0*approx_results.bias_pedestal/approx_results.gain);      //Integral from bias pedestal to infinity
    approx_results.N_light_pix=area_under_exponential_curve;

    approx_results.actual_number_pix_events_recorded = sum_counts_in_histogram(
                histogram_bins);

    approx_results.photon_event_freq =
        approx_results.N_light_pix /
        (approx_results.N_light_pix+approx_results.N_dark_pix);

//    approx_results.N_serial_CIC_pix =
//        approx_results.actual_number_pix_events_recorded -
//        (approx_results.N_dark_pix+approx_results.N_light_pix);

//    approx_results.serial_CIC_rate =
//            approx_results.N_serial_CIC_pix /
//            approx_results.actual_number_pix_events_recorded;


//    cerr<<"Approx fit:" << approx_results<<endl;

    if (perform_advanced_fit) {
//    cerr<<"Running advanced fit...";
        //Put in initial estimates:
        GainData full_fit = fit_full_CCD_model(histogram_bins,
                                                approx_results,
                                                output_to_screen
//                , sqrt_serial_CIC_gain_mode
                                               );

        return full_fit;
//        cerr<<"Done"<<endl;
    } else {

        return approx_results;
    }
//    vector<double> full_pars = gauss_pars; // zero_val, peak_count, sigma,
//
//    full_pars.push_back(exp_pars[0]); //photon_amplitude
//    full_pars.push_back(exp_pars[1]); //photon_gain,
//    full_pars.push_back(0.00); //serial_CIC_amplitude;
//    approx_results.photon_event_freq = area_under_exponential_curve / (area_under_gaussian_curve + area_under_exponential_curve);

}

//==================================================================================
GainData::GainData()
    : bias_pedestal(-1),readout_sigma(-1),gain(-1),
      N_dark_pix(-1), N_light_pix(-1), N_serial_CIC_pix(-1),
      photon_event_freq(-1), serial_CIC_rate(-1),
      actual_number_pix_events_recorded(0)
{}

void GainData::write_to_file(const std::string& filename)
{
    ofstream outputfile(filename.c_str());
    if (outputfile) {
        outputfile<<"#"+get_column_headers()<<endl;
        outputfile<<*this<<endl;
    } else { throw runtime_error("Cannot write to gain_inf file" + filename); }
}

GainData::GainData(const string& filename)
{
    using namespace string_utils;
    ifstream datafile(filename.c_str());
    vector<string> lines;
    string line;
    while (getline(datafile,line)) { lines.push_back(line); }
//             for (size_t i=0; i!=lines.size();++i) cout <<lines[i]<<endl;
    if (lines.size()!=2) { throw runtime_error("\nIncorrect file format for gain info: " + filename +"\n"); }
    vector<string> line_segments = tokenize_and_strip_spaces(lines[1], ",");
    if (line_segments.size()!=4) { throw runtime_error("\nIncorrect file format for gain info: " + filename +"\n"); }
    bias_pedestal=atof(line_segments[0]);
    readout_sigma=atof(line_segments[1]);
    gain=atof(line_segments[2]);
    photon_event_freq=atof(line_segments[3]);
    serial_CIC_rate=atof(line_segments[4]);
}

std::ostream& operator<<(std::ostream& os, const GainData& g_inf)
{
    os<<g_inf.bias_pedestal<<" "<<g_inf.readout_sigma<<" "
      <<g_inf.gain<<" "
      <<g_inf.photon_event_freq<<" "
      <<g_inf.serial_CIC_rate<<" "
//            <<g_inf.serial_CIC_gain<<" "
      <<g_inf.actual_number_pix_events_recorded<<" "
      ;
    return os;
}

std::string GainData::get_column_headers()
{
    return "bias_pedestal ""readout_sigma "
           "gain "
           "photon_event_frequency "
           "serial_CIC_rate "
//            "serial_CIC_gain "
           "N_pix_recorded "
           ;
}

//==================================================================================

template<typename input_datatype>
CcdImage<input_datatype>& normalise_CCD_with_uniform_gain(CcdImage<input_datatype>& input,
        const GainData& det_info)
{
    input.pix-=det_info.bias_pedestal;
    input.pix/=det_info.gain;
    return input;
}
template
CcdImage<float>& normalise_CCD_with_uniform_gain(CcdImage<float>& input,
        const GainData& det_info);

template
CcdImage<double>& normalise_CCD_with_uniform_gain(CcdImage<double>& input,
        const GainData& det_info);


//double_bitmap& normalise_CCD_with_gain_map(double_bitmap& input, const GainData& det_info, const double_bitmap& gain_map){
////    input.add_normalisation_flag();
//    input-=det_info.bias_pedestal;
//    input/=gain_map;
//    return input;
//}


//==================================================================================

CcdImage<float> threshold_bitmap(const CcdImage<float>& input,
                                 const double threshold_level)
{
    CcdImage<float> thresholded(input);
//    assert(input.key_exists("CCDNORMD"));
    for (PixelIterator i(input.pix.range()); i!=i.end; ++i) {
        input.pix(i)>threshold_level ? thresholded.pix(i)=1.0f : thresholded.pix(i)=0.0f;
    }
//    thresholded.add_keyword("THRESHED", string_utils::ftoa(threshold_level),"Datacount Level at which data was thresholded");
    return thresholded;
}


CcdImage<float> create_threshold_mask(const CcdImage<float>& input,
                                      const float threshold_count)
{
    CcdImage<float> mask(input);
    mask.pix.assign(1.0);
//     if (!input.key_exists("CCDNORMD")) throw logic_error("create_threshold_mask takes normalised input");
    for (PixelIterator pix(input.pix.range()); pix!=pix.end; ++pix) {
        if (input.pix(pix)>threshold_count) {
            //If bright area
            PixelRange mask_region(pix,pix);
            mask_region = PixelRange::pad(mask_region, 4); //FIXME - hardcoded radius
            mask_region = PixelRange::overlap(mask_region, mask.pix.range());
            for (PixelIterator i(mask_region); i!=i.end; ++i) { mask.pix(i) = 0.0; }
        }
    }
    return mask;
}



GainData analyse_histogram(const map<int, long>& full_data_histogram ,
                            const string& output_dir)
{
    gain_utils::write_histogram_data_to_file(full_data_histogram,
            output_dir+"short_hist.dat");
    vector<double> hump_pars= gain_utils::fit_readout_gaussian(full_data_histogram, true);

    map<int,long> gauss_data = gain_utils::get_gaussian_fitting_section(full_data_histogram);
    gain_utils::GaussianHistogramFitFCN gauss_fitter(gauss_data);
    map<int,double> gauss_fitted_bit = gauss_fitter.get_model_histogram(hump_pars);
    gain_utils::write_histogram_data_to_file(gauss_data, output_dir+"hump_data.dat");
    gain_utils::write_histogram_data_to_file(gauss_fitted_bit,
            output_dir+"hump_model_fit.dat");

    map<int,long> threshed_data = gain_utils::threshold_above_count(full_data_histogram, 100);

    gain_utils::GaussianHistogramFitFCN gauss_full(threshed_data);
    map<int,double> gauss_model = gauss_full.get_model_histogram(hump_pars);
    gain_utils::write_histogram_data_to_file(gauss_model, output_dir+"hump_model_full.dat");

    vector<double> tail_pars = gain_utils::fit_exponential_tail(full_data_histogram,
                               hump_pars.back(), true);

    map<int, long> fitted_tail_data = gain_utils::get_tail_fitting_section(
                                          full_data_histogram, hump_pars.back());
    gain_utils::ExponentialHistogramFitFCN tail_fitter(fitted_tail_data);
    map<int, double> fitted_tail_model = tail_fitter.get_model_histogram(tail_pars);
    gain_utils::write_histogram_data_to_file(fitted_tail_data, output_dir+"tail_data.dat");
    gain_utils::write_histogram_data_to_file(fitted_tail_model,
            output_dir+"tail_model_fit.dat");

    map<int, long> postive_section = gain_utils::get_points_at_values_greater_than(
                                         full_data_histogram, hump_pars.front()+hump_pars.back()); //
    gain_utils::ExponentialHistogramFitFCN tail_model(postive_section);
    map<int, double>full_tail_model = tail_model.get_model_histogram(tail_pars);
    gain_utils::write_histogram_data_to_file(full_tail_model,
            output_dir+"tail_model_full.dat");

    gain_utils::write_histogram_data_to_file(gain_utils::sum_histograms(full_tail_model,
            gauss_model), output_dir+"combined_model.dat");

    map<int,long> unfitted_data;
    for (map<int,long>::const_iterator it=full_data_histogram.begin();
            it!=full_data_histogram.end(); ++it) {
        int value = it->first;
        if (gauss_data.count(value)==0 && fitted_tail_data.count(value)==0) {
            unfitted_data[value]=it->second;
        }
    }

    gain_utils::write_histogram_data_to_file(unfitted_data, output_dir+"unfitted_data.dat");

    GainData results;
    results.bias_pedestal = hump_pars.front();
    results.readout_sigma = hump_pars.back();
    results.gain = tail_pars.back();

    double peak_gaussian_count = hump_pars[1];
    double deduced_proportion_dark_pix = peak_gaussian_count*results.readout_sigma * sqrt(
            2*M_PI); //area under a gaussian
    double deduced_proportion_photon_pix = results.gain * tail_pars.front()* exp(
            -1.0*results.bias_pedestal/results.gain);      //Amplitude by gain (amplitude corrected for origin shift to gaussian peak)
    results.photon_event_freq = deduced_proportion_photon_pix / (deduced_proportion_dark_pix+
                                deduced_proportion_photon_pix);
    results.N_light_pix=deduced_proportion_photon_pix;
    results.N_dark_pix=deduced_proportion_dark_pix;
    return results;
}

//=====================================================================================================================
}//end namespace
}

