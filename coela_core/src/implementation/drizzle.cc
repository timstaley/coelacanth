
#include "../drizzle.h"
#include <iostream>
using namespace std;

namespace coela {
namespace drizzle {
//=====================================================================================================================================================================



double pixel_overlap_as_fraction_of_drop_region(const PixelBoxRegion& drop_pixel_outline,
        const double drop_pixel_area_in_output_pixels,
        const PixelIndex output_pixel)
{

    //std
//    PixelBoxRegion output_pixel_region(
//            output_pixel.x -1.0, output_pixel.y-1.0,
//            output_pixel.x, output_pixel.y );
//
//    PixelBoxRegion overlap = PixelBoxRegion::overlapping_region(drop_pixel_outline, output_pixel_region);
//    return overlap.area() / drop_pixel_outline.area();
    //----------------
    //optimized:
    PixelBoxRegion overlap(max(output_pixel.x-1.0, drop_pixel_outline.low.x),
                           max(output_pixel.y-1.0, drop_pixel_outline.low.y),
                           min((double)output_pixel.x, drop_pixel_outline.high.x),
                           min((double)output_pixel.y, drop_pixel_outline.high.y));
    return overlap.area() / drop_pixel_area_in_output_pixels;


}

///This must be calculated manually for a pixarray, where we would normally use conversion between CCD reference frame for image.
PixelRange calculate_valid_input_range(const PixelBoxRegion& input_outline,
                                       const PixelBoxRegion& output_outline,
                                       const PixelShift& translation_vector_at_input_pixel_scale,
                                       const double drizzle_scale_factor)
{

    PixelShift output_diagonal_in_input_space
        = (output_outline.high - PixelPosition::origin) * drizzle_scale_factor;

    PixelBoxRegion valid_input_region = PixelBoxRegion(
                                            translation_vector_at_input_pixel_scale*(-1.) + PixelPosition::origin ,
                                            (translation_vector_at_input_pixel_scale*(-1.) + output_diagonal_in_input_space) +
                                            PixelPosition::origin);


    valid_input_region = PixelBoxRegion::overlapping_region(
                             input_outline,
                             valid_input_region);

    //This line determines whether we drizzle partially overlapping input pixels;
    // shrink to pixel boundaries will only drizzle totally overlapping input pixels.
    // expand to pixel boundaries will also drizzle partially overlapping input pixels.

    //We shrink to only total overlap, so we do not have to range check the calculated output pixels per input pixel.
    return valid_input_region.shrink_to_pixel_boundaries().bounded_pixels();

}

template<typename input_array_type, typename output_array_type>
void translate_and_drizzle_frame(
    const PixelArray2d<input_array_type>& input,
    const PixelArray2d<input_array_type>& input_weights,
    PixelArray2d<output_array_type>& output, PixelArray2d<output_array_type>& output_weights,
    const PixelShift& translation_vector_at_input_pixel_scale,
    const double drizzle_scale_factor, const double drizzle_pixel_fraction)
{


    //Since we have a long parameter list, it's probably sensible to check we haven't been accidentally given duplicate args.
//    assert(&input != &output); //Doesn't compile if input type is not double - have to hope nobody tries supplying same input and output...
    assert(&input !=
           &output); //Actually, a cast should do it, since we're just comparing memory addresses anyway
    assert(&input_weights != &output_weights);

    assert(input.range()==input_weights.range());
    assert(output.range()==output_weights.range());

    assert(translation_vector_at_input_pixel_scale.is_valid());


    PixelRange valid_input_range  = calculate_valid_input_range(
                                        PixelBoxRegion::pixel_box_outline(input.range()),
                                        PixelBoxRegion::pixel_box_outline(output.range()),
                                        translation_vector_at_input_pixel_scale, drizzle_scale_factor);


    PixelShift drop_pixel_half_diagonal_in_output_pixels(drizzle_pixel_fraction*0.5 /
            drizzle_scale_factor, drizzle_pixel_fraction*0.5 / drizzle_scale_factor);

    double drop_pixel_area_as_multiple_of_output_pixel =
        drizzle_pixel_fraction*drizzle_pixel_fraction /
        (drizzle_scale_factor*drizzle_scale_factor);

    for (PixelIterator i(valid_input_range); !(i == i.end); ++i) {
        //calculate centre of shifted input pixel
        PixelShift drop_centre_offset_from_output_origin((PixelPosition::centre_of_pixel(
                    i) - PixelPosition::origin) + translation_vector_at_input_pixel_scale);

        drop_centre_offset_from_output_origin/=drizzle_scale_factor; //transform to output scale.

        //Calculate the outline of the dropped pixel
        const PixelBoxRegion drop_outline((drop_centre_offset_from_output_origin -
                                           drop_pixel_half_diagonal_in_output_pixels) + PixelPosition::origin ,
                                          (drop_centre_offset_from_output_origin + drop_pixel_half_diagonal_in_output_pixels) +
                                          PixelPosition::origin);

        //NB make sure not to modify the original drop outline here

        //original:
        //PixelRange output_pixels_covered = ( PixelBoxRegion( drop_outline).expand_to_pixel_boundaries() ).bounded_pixels();
        //---------------------------
        //Optimized - conversion from float to int is same as floor for positive values
        //(But not ceil!)
        PixelRange output_pixels_covered(drop_outline.low.x+1, drop_outline.low.y+1,
                                         ceil(drop_outline.high.x), ceil(drop_outline.high.y));

//           double input_pixel_value = input(i);
        double input_pixel_weight = input_weights(i);
        double weighted_input_val = input(i)*input_pixel_weight;

        for (PixelIterator o(output_pixels_covered); o!=o.end; ++o) {
            double drop_fraction = pixel_overlap_as_fraction_of_drop_region(drop_outline,
                                   drop_pixel_area_as_multiple_of_output_pixel,o);
            output(o) += drop_fraction * weighted_input_val;
//               output(o) += drop_fraction * input_pixel_value * input_pixel_weight;
            output_weights(o)+= drop_fraction * input_pixel_weight;
        }

    }
}

template<typename input_array_type, typename output_array_type>
void dual_translate_and_drizzle_frame(
    const PixelArray2d<input_array_type>& input,
    const PixelArray2d<input_array_type>& dual_input,
    const PixelArray2d<input_array_type>& input_weights,
    PixelArray2d<output_array_type>& output, PixelArray2d<output_array_type>& dual_output,
    PixelArray2d<output_array_type>& output_weights,
    const PixelShift& translation_vector_at_input_pixel_scale,
    const double drizzle_scale_factor, const double drizzle_pixel_fraction)
{


    assert(&input != &output);
    assert(&input_weights != &output_weights);
    assert(&dual_output != &output);

    assert(input.range()==input_weights.range());
    assert(output.range()==output_weights.range());
    assert(input.range()==dual_input.range());
    assert(output.range()==dual_output.range());

    assert(translation_vector_at_input_pixel_scale.is_valid());


    PixelRange valid_input_range  = calculate_valid_input_range(
                                        PixelBoxRegion::pixel_box_outline(input.range()),
                                        PixelBoxRegion::pixel_box_outline(output.range()),
                                        translation_vector_at_input_pixel_scale, drizzle_scale_factor);


    PixelShift drop_pixel_half_diagonal_in_output_pixels(drizzle_pixel_fraction*0.5 /
            drizzle_scale_factor, drizzle_pixel_fraction*0.5 / drizzle_scale_factor);

    double drop_pixel_area_as_multiple_of_output_pixel =
        drizzle_pixel_fraction*drizzle_pixel_fraction /
        (drizzle_scale_factor*drizzle_scale_factor);

    for (PixelIterator i(valid_input_range); i != i.end; ++i) {

        //calculate centre of shifted input pixel
        PixelShift drop_centre_offet_from_output_origin((PixelPosition::centre_of_pixel(
                    i) - PixelPosition::origin) + translation_vector_at_input_pixel_scale);

        drop_centre_offet_from_output_origin/=drizzle_scale_factor; //transform to output scale.

        //Calculate the outline of the dropped pixel
        const PixelBoxRegion drop_outline((drop_centre_offet_from_output_origin -
                                           drop_pixel_half_diagonal_in_output_pixels) + PixelPosition::origin ,
                                          (drop_centre_offet_from_output_origin + drop_pixel_half_diagonal_in_output_pixels) +
                                          PixelPosition::origin);

        //original:
        //PixelRange output_pixels_covered = ( PixelBoxRegion( drop_outline).expand_to_pixel_boundaries() ).bounded_pixels();
        //---------------------------
        //Optimized - conversion from float to int is same as floor for positive values
        //(But not ceil!)
        PixelRange output_pixels_covered(drop_outline.low.x+1, drop_outline.low.y+1,
                                         ceil(drop_outline.high.x), ceil(drop_outline.high.y));

//           double input_pixel_value = input(i);
//           double dual_input_pixel_value = dual_input(i);

        double input_pixel_weight = input_weights(i);
        double weighted_input_val = input(i)*input_pixel_weight;
        double weighted_dual_input_val = dual_input(i)*input_pixel_weight;

        for (PixelIterator o(output_pixels_covered); o!=o.end; ++o) {
            double drop_fraction = pixel_overlap_as_fraction_of_drop_region(drop_outline,
                                   drop_pixel_area_as_multiple_of_output_pixel,o);
//               output(o) +=     input_pixel_value       * drop_fraction *  input_pixel_weight;
            output(o) +=     weighted_input_val * drop_fraction ;
//               dual_output(o)+= dual_input_pixel_value * drop_fraction *  input_pixel_weight;
            dual_output(o)+= weighted_dual_input_val* drop_fraction ;
            output_weights(o)+= drop_fraction * input_pixel_weight;
        }

    }
}

template<typename array_type>
PixelArray2d<array_type> unweight_drizzle_results(const PixelArray2d<array_type>& output,
        const PixelArray2d<array_type>& output_weights)
{
    PixelArray2d<array_type> unweighted(output);
    double weight;
    for (PixelIterator i(unweighted.range()); i!=i.end; ++i) {
        weight = output_weights(i);
        if (weight != 0.0) {
            unweighted(i) /= weight ;
        } else { unweighted(i)= 0.0 ; }
    }
    return unweighted;
}

template
PixelArray2d<float> unweight_drizzle_results(const PixelArray2d<float>& output,
        const PixelArray2d<float>& output_weights);
template
PixelArray2d<double> unweight_drizzle_results(const PixelArray2d<double>& output,
        const PixelArray2d<double>& output_weights);


//=================================================================================================================
//Instantiate and compile the templated functions for various different input types:
//-----------------------------------------------------------------------------------------------------
//void translate_and_drizzle_frame

template //<double> input array type
void translate_and_drizzle_frame(
    const PixelArray2d<double>& input, const PixelArray2d<double>& input_weights,
    PixelArray2d<double>& output, PixelArray2d<double>& output_weights,
    const PixelShift& translation_vector_at_input_pixel_scale,
    const double drizzle_scale_factor, const double drizzle_pixel_fraction);

template //<float> input array type
void translate_and_drizzle_frame(
    const PixelArray2d<float>& input, const PixelArray2d<float>& input_weights,
    PixelArray2d<float>& output, PixelArray2d<float>& output_weights,
    const PixelShift& translation_vector_at_input_pixel_scale,
    const double drizzle_scale_factor, const double drizzle_pixel_fraction);

//-----------------------------------------------------------------------------------------------------
//void dual_translate_and_drizzle_frame

template //<double> input array type
void dual_translate_and_drizzle_frame(
    const PixelArray2d<double>& input, const PixelArray2d<double>& dual_input ,
    const PixelArray2d<double>& input_weights,
    PixelArray2d<double>& output, PixelArray2d<double>& dual_output,
    PixelArray2d<double>& output_weights,
    const PixelShift& translation_vector_at_input_pixel_scale,
    const double drizzle_scale_factor, const double drizzle_pixel_fraction);

template //<float> input array type
void dual_translate_and_drizzle_frame(
    const PixelArray2d<float>& input, const PixelArray2d<float>& dual_input ,
    const PixelArray2d<float>& input_weights,
    PixelArray2d<float>& output, PixelArray2d<float>& dual_output,
    PixelArray2d<float>& output_weights,
    const PixelShift& translation_vector_at_input_pixel_scale,
    const double drizzle_scale_factor, const double drizzle_pixel_fraction);

//=====================================================================================================================================================================
} //end namespace coela::drizzle
}//end namespace coela
