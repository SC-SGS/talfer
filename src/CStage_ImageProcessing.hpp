/*
 * CStage_ImageInput.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 *      Author: Marcel Schneider <marcel.schneider@studi.informatik.uni-stuttgart.de> 
 *		Author: Gregor Schwarz <gregor.schwarz@tum.de>
 */

#ifndef CSTAGE_IMAGEPROCESSING_HPP_
#define CSTAGE_IMAGEPROCESSING_HPP_

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "SDL/SDL_image.h"
#include "CDataArray2D.hpp"
#include "CDataDrawingInformation.hpp"

#include "image_parameters.h"

#include <stdint.h>
#include "imageops.c"
#ifdef USE_FFT
#include <fftw3.h>
#define ALLOC fftw_malloc
#define FREE fftw_free
#else
#define ALLOC malloc
#define FREE free
#endif
#include "hsv.h"

#define SQUARE(x) (x)*(x)
#define CUBE(x)   (x)*(x)*(x)

/**
 * class providing static image input
 *
 * the image is send to the pipeline during each main loop iteration
 */
class CStage_ImageProcessing: public CPipelineStage {
	/**
	 * global parameters
	 */
	CParameters &cParameters;

	/**
	 * input image
	 */
	CDataArray2D<unsigned char, 3> input_cDataArray2D;

	/**
	 * last input image
	 */
	CDataArray2D<unsigned char, 3> lastframe_cDataArray2D;

	/**
	 * processed image
	 */
	CDataArray2D<unsigned char, 3> output_cDataArray2D_uc3;

	/**
	 * processed flagfield
	 */
	CDataArray2D<unsigned char, 1> &output_cDataArray2D_uc1;

	/**
	 * Intermediate flagfield, before masking out areas.
	 */
	CDataArray2D<unsigned char, 1> output_flagfield;

	/**
	 * Internal image buffer 
	 */
	pixel_t *current;

	/**
	 * next can be used as a temporary buffer.
	 */
	pixel_t *next;

	/**
	 * Depthmap for determining whether image pixels are valid.
	 */
	// CDataArray2D<unsigned char,1> depth;

	/**
	 * internal buffer dimensions.
	 */
	int w, h;

	/**
	 * The internal buffer is padding pixels larger on all sides to simplify
	 * processing.
	 */
	int padding;

	/**
	 * Buffer used by hysteresis filter
	 */
	bool *old_flags;

	/**
	 * Hysteresis thresholds, automatically adapted; two for each channel
	 */
	int thresholds[6];

	/**
	 * tells whether current input image is usable (for hand detection)
	 */
	bool transmit_output;

#ifdef USE_FFT
	/**
	 * fftw_plan for the image to freq fftw
	 */
	fftw_plan fft_tofreq;

	/**
	 * fftw_plan for the freq to image fftw
	 */
	fftw_plan fft_toimg;

	/**
	 * DoG kernel in freq. domain
	 */
	fftw_complex *freq_kernel;

	/**
	 * FFT freq. domain result
	 */
	fftw_complex *freq_res;

#endif

public:
	/**
	 * constructor
	 */
	CStage_ImageProcessing(CParameters &i_cParameters) :
			CPipelineStage("ImageProcessing"), cParameters(i_cParameters), output_cDataArray2D_uc1(
					output_flagfield), current(NULL), next(NULL), old_flags(
					NULL), transmit_output(true) {
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push() {
		if (cParameters.stage_imageprocessing_output_flagfield)
			CPipelineStage::pipeline_push(
					(CPipelinePacket&) output_cDataArray2D_uc1);
		else
			CPipelineStage::pipeline_push(
					(CPipelinePacket&) output_cDataArray2D_uc3);
	}

private:

	/**
	 * determine change of two following frames and sets the parameters has_changed and transmit_output
	 */
	void check_change() {
		// special treatment for first frame of the game
		if (lastframe_cDataArray2D.isValidData() == false) {
			output_cDataArray2D_uc1.has_changed = true;
			return;
		}
		int len = input_cDataArray2D.width * input_cDataArray2D.height * 3;
		int diff = image_delta(&lastframe_cDataArray2D.getRef(0, 0),
				&input_cDataArray2D.getRef(0, 0), len, CHANGE_EPSILON);
		printf("%d ", diff);
		if (diff > (CHANGE_THRESHOLD * len)) {
			output_cDataArray2D_uc1.has_changed = true;
		} else {
			if (transmit_output == false) {
				transmit_output = true;
				return;
			}
			output_cDataArray2D_uc1.has_changed = false;
		}

		// don't update internal frame if there is too much change
		if (diff > (CHANGE_THRESHOLD_MAX * len)) {
			transmit_output = false;
		}

	}

#ifdef USE_FFT
	/**
	 * generate kernel for FFT gaussian blur
	 */
	void generate_kernel() {
		memset(current, 0, w * h * sizeof(pixel_t) * 3);
		current[h*w/2 + w/2] = 1;
		apply_filter_blur(DOG_FINE); // fine blur (just a single rect blur)
		pixel_t *orig = new pixel_t[w*h*3];
		memmove(orig, current, 3*w*h*sizeof(pixel_t));
		apply_filter_blur(DOG_COARSE);// rough blur
		for (int i = 0; i < 3*w*h; i++)
		next[i] = (orig[i] - current[i]);
		delete [] orig;
		//memmove(current, next, sizeof(pixel_t)*w*h);
		for(int y = 0; y < h; y++)
		for (int x = 0; x < w; x++) {
			int ly = (y - h/2);
			if (ly < 0) ly += h;
			current[ly * w + x] = next[y * w + ((x - w/2) % w)];
		}
		fftw_execute(fft_tofreq);
		memmove(current, freq_res, sizeof(fftw_complex)*w*(h/2+1));
	}
#endif

	/**
	 * Convert CDataArray2D representation to internal flat buffer format. 
	 * example:
	 * RRRRRRRRR
	 * BBBBBBBBB
	 * GGGGGGGGG
	 *
	 * The allocated buffers will be freed in from_internal_buffer. 
	 * Conversion to HSV Colorspace is not performed.
	 */
	void to_internal_buffer() {

		int org_w = input_cDataArray2D.width;
		int org_h = input_cDataArray2D.height;

		output_cDataArray2D_uc3.resize(org_w, org_h);

		padding = PADDING_WIDTH;
		int w = org_w + padding + padding;
		int h = org_h + padding + padding;
		this->w = w;
		this->h = h;
		if (current == NULL || next == NULL) {
			current = (pixel_t*) ALLOC(sizeof(pixel_t) * 3 * h * w);
			next = (pixel_t*) ALLOC(sizeof(pixel_t) * 3 * h * w);
#ifdef USE_FFT
			freq_kernel = (fftw_complex*) ALLOC(sizeof(fftw_complex) * (h/2+1) * w);
			freq_res = (fftw_complex*) ALLOC(sizeof(fftw_complex) * (h/2+1) * w);
			fft_tofreq = fftw_plan_dft_r2c_2d(h, w, current, freq_res, FFTW_ESTIMATE);
			fft_toimg = fftw_plan_dft_c2r_2d(h, w, freq_res, current, FFTW_ESTIMATE);
			generate_kernel();
			memmove(freq_kernel, current, sizeof(fftw_complex)*w*(h/2+1));
#endif
		}

		for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++) {
				uint8_t *inputPix = &input_cDataArray2D.getClampedRef(
						x - padding, y - padding);
				current[0 * w * h + y * w + x] = inputPix[0];
				current[1 * w * h + y * w + x] = inputPix[1];
				current[2 * w * h + y * w + x] = inputPix[2];
			}
	}

	/**
	 * Converts the internal representation back into CDataArray2D and
	 * frees the buffers.
	 */
	void from_internal_buffer() {
		for (int y = padding; y < input_cDataArray2D.height + padding; y++)
			for (int x = padding; x < input_cDataArray2D.width + padding; x++) {
				uint8_t *outputPix = &output_cDataArray2D_uc3.getRef(
						x - padding, y - padding);
				outputPix[0] = to_uint8(current[0 * w * h + y * w + x]);
				outputPix[1] = to_uint8(current[1 * w * h + y * w + x]);
				outputPix[2] = to_uint8(current[2 * w * h + y * w + x]);
			}
	}

	/**
	 * similar to from_internal_buffer but also shows the padding.
	 */
	void from_internal_buffer_all() {
		output_cDataArray2D_uc3.resize(w, h);
		for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++) {
				uint8_t *outputPix = &output_cDataArray2D_uc3.getRef(x, y);
				outputPix[0] = to_uint8(current[0 * w * h + y * w + x]);
				outputPix[1] = to_uint8(current[1 * w * h + y * w + x]);
				outputPix[2] = to_uint8(current[2 * w * h + y * w + x]);
			}
	}

	/**
	 * Exchange the two internal buffers.
	 */
	void swap() {
		pixel_t *t = current;
		current = next;
		next = t;
	}

public:
#if 0
	/**
	 * Do anisotropic diffusion to reduce noise. Similiar to blurring, but
	 * does not destroy/move edges.
	 * Probably too slow to be useful.
	 * TODO: parameters hardcoded here and in imageops.c:phi().
	 */
	void apply_filter_diffusion(int passes) {
		make_phi_tbl();
		for(int i = 0; i < passes; i++) {
			do_diff(&current[0 * w * h], &next[0 * w * h], w, h);
			do_diff(&current[1 * w * h], &next[1 * w * h], w, h);
			do_diff(&current[2 * w * h], &next[2 * w * h], w, h);
			swap();
		}
	}
#endif

	/**
	 * Approximates a gaussian blur by doing three fast box blurs (using
	 * integral images). Approximation should be within 3% [1]. 
	 * (see [1]).
	 * [1]: http://www.w3.org/TR/SVG/filters.html#feGaussianBlurElement
	 */
	void apply_filter_blur(int size, int passes = 3) {
		for (int i = 0; i < passes; i++) {
			integrate_image(&current[0 * w * h], &next[0 * w * h], w, h);
			integrate_image(&current[1 * w * h], &next[1 * w * h], w, h);
			integrate_image(&current[2 * w * h], &next[2 * w * h], w, h);
			rect_blur(&next[0 * w * h], &current[0 * w * h], size, w, h);
			rect_blur(&next[1 * w * h], &current[1 * w * h], size, w, h);
			rect_blur(&next[2 * w * h], &current[2 * w * h], size, w, h);
		}
	}

	/**
	 * applies difference-of-gauss-filter by subtracting a fine blurred from
	 * a rough blurrend image.
	 * TODO: Size parameters hardcoded here.
	 */
	void apply_filter_DoG() {
		apply_filter_blur(DOG_FINE); // fine blur (just a single rect blur)
		pixel_t *orig = new pixel_t[w * h * 3];
		memmove(orig, current, 3 * w * h * sizeof(pixel_t));
		apply_filter_blur(DOG_COARSE); // rough blur
		for (int i = 0; i < 3 * w * h; i++)
			next[i] = orig[i] - current[i] + 127;
		swap();
		delete[] orig;
	}

#ifdef USE_FFT	
	void apply_filter_fft_DoG() {
		do_fft_DoG(next);
		memmove(current, &current[w*h], sizeof(pixel_t)*w*h);
		do_fft_DoG(&next[w*h]);
		memmove(current, &current[2*w*h], sizeof(pixel_t)*w*h);
		do_fft_DoG(&next[2*w*h]);
		memmove(current, next, 3*sizeof(pixel_t)*w*h);
	}

	void do_fft_DoG(pixel_t * dest) {
		fftw_execute(fft_tofreq);
		fftw_complex *buf = freq_res;
		int len = w*(h/2+1)-1;
		for (int i = 0; i < w*(h/2+1); i++) {
			buf[i][0] = buf[i][0] * freq_kernel[i][0] - buf[i][1] * freq_kernel[i][1];
			buf[i][1] = buf[i][0] * freq_kernel[i][1] - buf[i][1] * freq_kernel[i][0];
		}
		fftw_execute(fft_toimg);
		double norm = 1.0/(w*h);
		len = w*h;
		for(int i = 0; i < w*h; i++) {
			dest[i] = current[len-i] * norm + 127;
		}
	}
#endif

	/**
	 * Very simple thresholding.
	 */
	void apply_filter_threshold() {
		int threshold = cParameters.stage_imageprocessing_threshold_value;
		for (int i = 0; i < w * h * 3; i++)
			next[i] = current[i] > threshold ? 255 : 0;
		swap();
	}

	/** 
	 * Adaptive threshold with hysteresis
	 */
	void apply_filter_hysteresis() {
		// special treatment for first frame
		if (old_flags == NULL) {
			old_flags = new bool[w * h * 3];
			apply_filter_threshold();
			for (int i = 0; i < w * h * 3; i++)
				old_flags[i] = current[i] != 0;
			for (int i = 0; i < 6; i += 2) {
				thresholds[i] = INITIAL_THRESHOLD;
				thresholds[i + 1] = INITIAL_THRESHOLD + INITIAL_HYSTERESIS;
			}
		} else { // for every other frame
			threshold_hysteresis(&current[0 * w * h], &next[0 * w * h],
					&old_flags[0 * w * h], &thresholds[0]);
			threshold_hysteresis(&current[1 * w * h], &next[1 * w * h],
					&old_flags[1 * w * h], &thresholds[2]);
			threshold_hysteresis(&current[2 * w * h], &next[2 * w * h],
					&old_flags[2 * w * h], &thresholds[4]);
			swap();
		}
	}
	/**
	 * Performs the thresholding on the given buffers.
	 */
	void threshold_hysteresis(pixel_t *current, pixel_t *next, bool *old,
			int *tresh) {
		int &threshold_low = tresh[0];
		int &threshold_high = tresh[1];
		for (int i = 0; i < w * h; i++) {
			bool newval = current[i]
					> (old[i] ? threshold_low : threshold_high);
			next[i] = newval ? 255 : 0;
			old[i] = newval;
		}
	}

	/**
	 * Generate the flagfield output values based on the filtered and thresholded values and the 
	 * color information of the original image.
	 */
	void apply_filter_flagfield() {
		output_flagfield.resize(input_cDataArray2D.width,
				input_cDataArray2D.height);

		for (int y = 0; y < input_cDataArray2D.height; y++) {
			for (int x = 0; x < input_cDataArray2D.width; x++) {
				RgbColor *d_o = (RgbColor*) &output_cDataArray2D_uc3.getRef(x,
						y);
				RgbColor *d_i = (RgbColor*) &input_cDataArray2D.getRef(x, y);
				uint8_t res = 0;
				if (d_o->r + d_o->g + d_o->b <= 255) { // at most one of three
					res = 1;
					HsvColor h = RgbToHsv(*d_i);
					if (COLOR_INPUT(h)) {
						if (h.s < INPUT_MIN_SATURATION)
							res = 1; // Wall, black
						else
							res = 2; // input, red
					} else if (COLOR_OUTPUT(h)) {
						if (h.s < OUTPUT_MIN_SATURATION)
							res = 1; // Wall, black
						else
							res = 3; // output, blue
					} else if (COLOR_TARGET(h)) {
						if (h.s < TARGET_MIN_SATURATION)
							res = 1; // Wall, black
						else
							res = 4; // output, blue
					}
				}
				output_flagfield.getRef(x, y, 0) = res;
			}
		}
	}

	/**
	 * This filter is applied after the flagfield filter. It removes
	 * the noise (black pixels) around the red, blue, and green parts of the picture.
	 * The parameter size sets the size of the area to be analyzed.
	 * Update of Pixels is done in place at the moment.
	 */
	void apply_flagfield_postprocess(int passes, int size, int redThreshold,
			int blueThreshold, int greenThreshold) {
		for (int i = 0; i < passes; i++) {
			for (int x = size; x < output_flagfield.width - size; x++) {
				for (int y = size; y < output_flagfield.height - size; y++) {
					uint8_t *flag = &output_flagfield.getRef(x, y);

					if (*flag == 1) {
						int numFlags[] = { 0, 0, 0, 0, 0 }; // white, black, red, blue, green
						count_flags(numFlags, size, x, y);

						// filter red
						if (numFlags[2] > redThreshold) {
							*flag = 2;
						}
						// filter blue
						else if (numFlags[3] > blueThreshold) {
							*flag = 3;
						}
						// filter green
						else if (numFlags[4] > greenThreshold) {
							*flag = 4;
						}
					}

				}
			}
		}
	}

	/**
	 * Returns an array indicating how many flags of a certain color
	 * could be found around the given anchor point. 
	 * The parameter size sets the size of the area to be analyzed.
	 * example: size = 2, A = anchor, x = flag
	 *   x x x x x
	 *   x x x x x
	 *   x x A x x
	 *   x x x x x
	 *   x x x x x
	 */
	void count_flags(int *numFlags, int size, int anchor_x, int anchor_y) {
		for (int x = anchor_x - size; x <= anchor_x + size; x++) {
			for (int y = anchor_y - size; y <= anchor_y + size; y++) {
				uint8_t flagColor = output_flagfield.getRef(x, y);
				numFlags[flagColor]++;
			}
		}
	}

	/**
	 * make sure that there are no single pixels sticking out or missing.
	 * This is neccessary for Navier Stokes.
	 *  XX XX       XXXXX
	 *  XXXXX  ==>  XXXXX
	 *
	 *   X
	 *  XXX    ==>  XXX
	 *  XXX		   XXX
	 */
	void apply_flagfield_pixelfix() {
		int count = 1;
		//printf("\nfixing: "); 
		while (count) {
			count = cleanup_pixels(&output_flagfield.getRef(0, 0),
					output_flagfield.width, output_flagfield.height);
			//printf(" %d", count);
		}
	}

	/** 
	 * Visualize the generated flagfield for debugging
	 */
	void apply_show_flagfield() {
		for (int y = 0; y < input_cDataArray2D.height; y++) {
			for (int x = 0; x < input_cDataArray2D.width; x++) {
				RgbColor *d_o = (RgbColor*) &output_cDataArray2D_uc3.getRef(x,
						y);
				uint8_t flag = output_cDataArray2D_uc1.getRef(x, y);

				switch (flag) {
				case 0:
					*d_o = rgb(255, 255, 255);
					break;
				case 1:
					*d_o = rgb(0, 0, 0);
					break;
				case 2:
					*d_o = rgb(255, 0, 0);
					break;
				case 3:
					*d_o = rgb(0, 0, 255);
					break;
				default:
					*d_o = rgb(0, 255, 0);
					break;
				}
			}
		}
	}

	RgbColor rgb(uint8_t r, uint8_t g, uint8_t b) {
		RgbColor t;
		t.r = r;
		t.g = g;
		t.b = b;
		return t;
	}

public:
	/**
	 * process incoming pipeline input.
	 *
	 * the only input processed so far is from the video output stage to
	 * draw something into the image.
	 */
	void pipeline_process_input(CPipelinePacket &i_cPipelinePacket) {
		printf("hallo\n");
		// we are currently only able to process "unsigned char,3" data arrays.
		if (i_cPipelinePacket.type_info_name
				!= typeid(CDataArray2D<unsigned char, 3> ).name()) {
			std::cerr
					<< "ERROR: Video Output is only able to process (char,3) arrays"
					<< std::endl;
			exit(-1);
		}

		// unpack data
		CDataArray2D<unsigned char, 3> *input = i_cPipelinePacket.getPayload<
				CDataArray2D<unsigned char, 3> >();

		// copy data to input array
		input_cDataArray2D.resize(input->width, input->height);
		input_cDataArray2D.loadData(input->data);

		// if there has been no change to the frame, push no frame
		check_change();
		if (output_cDataArray2D_uc1.has_changed == false
				&& (cParameters.stage_imageprocessing_filter_id > 4
						|| cParameters.stage_imageprocessing_filter_id == -1)) {
			//printf("nop-push\n");
			pipeline_push();
			return;
		}
		if ((!cParameters.stage_imageprocessing_freeze_output
				&& transmit_output == true)
				|| cParameters.stage_imageprocessing_filter_id < 4) {
			// convert to internal representation
			to_internal_buffer();

			// FILTER
			switch (cParameters.stage_imageprocessing_filter_id) {
			case 0:
				apply_filter_blur(20, 3);
				break;
			case 1:
				//apply_filter_threshold();
				apply_filter_DoG();
				break;
			case 2:
#ifdef USE_FFT
				apply_filter_DoG();
				//generate_kernel();
#endif
				break;
			case 3:
				//apply_filter_diffusion(5);
				break;
			case 4:
			case 5:
			case 6:
				//apply_filter_diffusion(); // does not really help
#ifdef USE_FFT
				apply_filter_fft_DoG();
#else
				apply_filter_DoG();
#endif
				apply_filter_hysteresis();
				break;

			default:	// no filter
				break;
			}
			// convert back to "normal" image
			from_internal_buffer();

			// applied in -p5 execution mode
			if (cParameters.stage_imageprocessing_output_flagfield
					|| cParameters.stage_imageprocessing_filter_id >= 5) {
				apply_filter_flagfield();
				swap(); // magic to keep FFT bufs in sync
				if (cParameters.stage_imageprocessing_filter_id == 5) {
					apply_flagfield_postprocess(POSTPROCESS_PARAMETERS); // passes, size, redThreshold, blueThreshold, greenThreshold
					apply_flagfield_pixelfix();
				}
				if (cParameters.stage_imageprocessing_filter_id >= 5)
					apply_show_flagfield();

			}
		}

		// push modifications on image to pipeline
		pipeline_push();

		// store last frame if new frame was pushed to pipeline
		if (output_cDataArray2D_uc1.has_changed == true) {
			//printf("change: %d\n", output_cDataArray2D_uc1.has_changed);
			input_cDataArray2D.swap(lastframe_cDataArray2D);
		}
	}

	void main_loop_callback() {
	}
};

#endif /* CSTAGE_IMAGEINPUT_HPP_ */
