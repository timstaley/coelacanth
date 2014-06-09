// Author: ts337
//
// Created on 21 May 2008, 16:50
//

#ifndef COELA_MY_STATS_H
#define COELA_MY_STATS_H

#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>

using std::vector;
using std::max;
namespace coela {

vector<double> median_filter(const vector<double>& input, const size_t neighbour_depth=1);
vector<double> box_smooth(const vector<double>& input, const size_t neighbour_depth=1);

template <class T> //NB returns double as returning int, say, would lose information
double vector_mean(const std::vector<T>& vec)
{
    double sum(
        0.0); //store the sum in a double so we don't lose any accuracy along the way (e.g. adding a big float to a little float and storing answer in a float loses information (i think?))
    for (size_t i(0); i<vec.size(); ++i) {
        sum+=vec[i];
    }
    return sum / vec.size();
};

template <class T>
double vector_squared_mean(const vector<T>& v)
{
    double sum(0.0);
    for (size_t i(0); i<v.size(); ++i) {
        sum+=v[i]*v[i];
    }
    return sum / v.size();
};

template <class T>
double vector_variance(const vector<T>& v)
{
    double sum(0.0), squares_sum(0.0);
    for (size_t i(0); i<v.size(); ++i) {
        sum+=v[i];
        squares_sum+=v[i]*v[i];
    }
    sum /= (T)v.size();  //Now equal to the mean
    squares_sum /= (T)v.size(); //Now equal to the mean of the squares
    return squares_sum - sum*sum; //mean of squares minus square of mean ==variance.
};

template <class InputIterator, class T>
double general_container_mean(InputIterator b, InputIterator e)
{
    double sum(0.0);
    int n_vals=0;
    if (b==e) { return 0.; }
    for (; b!=e; ++b) {
        sum+=(*b);
        ++n_vals;
    }
    sum /= n_vals;  //Now equal to the mean
    return sum;
};

template <class InputIterator>
double container_variance(InputIterator b, InputIterator e)
{
    double sum(0.0), squares_sum(0.0);
    int n_vals=0;
    if (b==e) { return 0.; }
    for (; b!=e; ++b) {
        sum+=(*b);
        squares_sum+=(*b) * (*b);
        ++n_vals;
    }
    sum /= n_vals;  //Now equal to the mean
    squares_sum /= n_vals; //Now equal to the mean of the squares
    return squares_sum - sum*sum; //mean of squares minus square of mean ==variance.
};

template <class T>
double vector_std_dev(const vector<T>& v)
{
    return sqrt(vector_variance(v));
};

template <class T>
T vector_rms_from_zero(const vector<T>& v)
{
    vector<T> squares(v.begin(), v.end());
    for (size_t i(0); i!=squares.size(); ++i) {
        squares[i]*=squares[i];
    }
    return sqrt(vector_mean(squares));
};

template <class T>
T vector_sigma_clipped_mean(const std::vector<T>& vec)
{

    double mean = vector_mean(vec);
    double squared_mean = vector_squared_mean(vec);
    double sd = sqrt(squared_mean - mean*mean);

    vector<T> onesigma;
    onesigma.reserve(vec.size());

    for (size_t i(0); i < vec.size(); ++i) {
        if (fabs(vec[i]- mean) < sd)  {
            onesigma.push_back(vec[i]);
        }
    }

    return vector_mean(onesigma);
};

template <class T>
double vector_3sigma_clipped_mean(const std::vector<T>& vec)
{

    double mean = vector_mean(vec);
    double squared_mean = vector_squared_mean(vec);
    double sd = sqrt(squared_mean - mean*mean);

    vector<T> threesigma;
    threesigma.reserve(vec.size());

    for (size_t i(0); i < vec.size(); ++i) {
        if (fabs(vec[i]- mean) < 3.0*sd)  {
            threesigma.push_back(vec[i]);
        }
    }

    return vector_mean(threesigma);
};

template <class T>
T vector_twice_sigma_clipped_mean(const std::vector<T>& vec)
{

    double mean = vector_mean(vec);
    double squared_mean = vector_squared_mean(vec);
    double sd = sqrt(squared_mean - mean*mean);

    vector<T> one_sigma_clip, twice_sigma_clip;
    one_sigma_clip.reserve(vec.size());
    twice_sigma_clip.reserve(vec.size());

    for (size_t i(0); i < vec.size(); ++i) {
        if (fabs(vec[i]- mean) < sd)  {
            one_sigma_clip.push_back(vec[i]);
        }
    }

    mean = vector_mean(one_sigma_clip);
    squared_mean = vector_squared_mean(one_sigma_clip);
    sd = sqrt(squared_mean - mean*mean);

    for (size_t i(0); i < one_sigma_clip.size(); ++i) {
        if (fabs(one_sigma_clip[i]- mean) < sd)  {
            twice_sigma_clip.push_back(vec[i]);
        }
    }
    return vector_mean(twice_sigma_clip);

};


template <class T>
T vector_triple_sigma_clipped_mean(const std::vector<T>& vec)
{

    double mean = vector_mean(vec);
    double squared_mean = vector_squared_mean(vec);
    double sd = sqrt(squared_mean - mean*mean);

    vector<T> one_sigma_clip, twice_sigma_clip;
    one_sigma_clip.reserve(vec.size());
    twice_sigma_clip.reserve(vec.size());

    for (size_t i(0); i < vec.size(); ++i) {
        if (fabs(vec[i]- mean) < sd)  {
            one_sigma_clip.push_back(vec[i]);
        }
    }

    mean = vector_mean(one_sigma_clip);
    squared_mean = vector_squared_mean(one_sigma_clip);
    sd = sqrt(squared_mean - mean*mean);

    for (size_t i(0); i < one_sigma_clip.size(); ++i) {
        if (fabs(one_sigma_clip[i]- mean) < sd)  {
            twice_sigma_clip.push_back(vec[i]);
        }
    }

    vector <T> triple_sigma_clip;
    mean = vector_mean(twice_sigma_clip);
    squared_mean = vector_squared_mean(twice_sigma_clip);
    sd = sqrt(squared_mean - mean*mean);

    for (size_t i(0); i < twice_sigma_clip.size(); ++i) {
        if (fabs(twice_sigma_clip[i]- mean) < sd)  {
            triple_sigma_clip.push_back(vec[i]);
        }
    }


    return vector_mean(triple_sigma_clip);

};

//NB Takes a non-const argument since the vector will be re-arranged!
template <class T>
double vector_median(std::vector<T>& v)
{
//        assert(!v.empty());
//        sort(v.begin(), v.end());
//        if (v.size()%2==0){
//            return  (v[v.size()/2] + v[v.size()/2 + 1]) /2.0;
//        }
//        else
//        return v[v.size()/2];
    assert(!v.empty());
    size_t median_index = v.size()/2 + 1
                          ; //to be exact we should use the size/2+1 element if size is odd, and avg size/2, size/2 +1 elements if even...
    //but this will do for relatively large / smooth datasets
    nth_element(v.begin(), v.begin()+median_index, v.end());
    return v[median_index];
}

template <class T>
double vector_first_quartile(std::vector<T>& v)
{
    assert(!v.empty());
    sort(v.begin(), v.end());
    return v[v.size()/4];
}

template <class T>
T vector_third_quartile(std::vector<T>& v)
{
    assert(!v.empty());
    sort(v.begin(), v.end());
    return v[v.size()*3/4];
}

template <class T>
double vector_tenth_percentile(std::vector<T>& v)
{
    assert(!v.empty());
    size_t tenth_percentile_index = v.size()/10;
    nth_element(v.begin(), v.begin()+tenth_percentile_index, v.end());
    return v[tenth_percentile_index];
}

template <class T>
T vector_ninetieth_percentile(std::vector<T>& v)
{
    assert(!v.empty());
    sort(v.begin(), v.end());
    return v[(9*v.size())/10];
}

template <class T>
double vector_percentile(std::vector<T>& v, const double percentile_as_a_proportion)
{
    assert(!v.empty());
    assert(percentile_as_a_proportion>=0.0 && percentile_as_a_proportion<=1.0);
    size_t percentile_index= v.size()*percentile_as_a_proportion;
    if (percentile_index==v.size()) { percentile_index=v.size()-1; }
    nth_element(v.begin(), v.begin()+percentile_index, v.end());
    return v[percentile_index];
}


template <class T>
vector<T> clip_by_limit_about_value(const vector<T> &vec, const T limit, const T value)
{
    vector<T> clipped;
    clipped.reserve(vec.size());
    for (size_t i(0); i < vec.size(); ++i) {
        if (fabs(vec[i]- value) < limit)  {
            clipped.push_back(vec[i]);
        }
    }
    return clipped;
};

template <class T>
float clip_until_mode_convergence(const vector<T> &vec)
{
//     std::cout <<"\n Clipping..." <<std::endl;
    float mean = vector_mean(vec);
    float median = vector_median(vec);
    float squares_mean;// =vector_squared_mean(vec);
    float SD;// = sqrt(squares_mean - mean*mean);
    float mode_new = 2.5*median - 1.5*mean;
//    std::cout <<"\nNew mode estimate" <<mode_new <<std::endl;
    float mode_old;
    vector<T> temp(vec);
    do {
        mode_old = mode_new;
        squares_mean= vector_squared_mean(temp);
        SD = sqrt(squares_mean - mean*mean);
        temp = clip_by_limit_about_value(temp, 3*SD, median);
        mean = vector_mean(temp);
        median = vector_median(temp);
        mode_new = 2.5*median - 1.5*mean;
//  std::cout <<"New mode estimate" <<mode_new <<std::endl;

    } while (fabs(mode_new - mode_old) > 0.05*fabs(mode_old));

    return mean;
};



template <class T>
vector<T> median_filter(const vector<T>& input, const size_t neighbour_depth)
{
    vector<T> filtered;

    for (size_t input_index=0; input_index!=input.size(); input_index++) {
        size_t neighbour_start=max((signed long)input_index-(signed long)neighbour_depth, 0l);
        size_t neighbour_end = min((size_t)input_index+neighbour_depth +1 , input.size());
        vector<T> neighbours(input.begin()+neighbour_start,input.begin()+neighbour_end);
        filtered.push_back(vector_median(neighbours));
    }
    return filtered;
}

template <class T>
vector<T> box_smooth(const vector<T>& input, const size_t neighbour_depth)
{
    vector<T> filtered;

    for (size_t input_index=0; input_index!=input.size(); input_index++) {
        size_t neighbour_start=max((signed long)input_index-(signed long)neighbour_depth, 0l);
        size_t neighbour_end = min((size_t)input_index+neighbour_depth +1 ,
                                   input.size()); //nb +1 since we are marking the end point, not the last data point.
        vector<T> neighbours(input.begin()+neighbour_start,input.begin()+neighbour_end);
        filtered.push_back(vector_mean(neighbours));
    }
    return filtered;
}


}//end namespace coela

#endif  /* _MY_UTILS_H */

