/*
 * CParameters.hpp
 *
 *  Created on: Jul 6, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CPARAMETERS_HPP_
#define CPARAMETERS_HPP_

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include "CSDLInterface.hpp"

class CParameters {
public:
	/**
	 * this parameter is tested during each pipeline step.
	 * in case that it's true, the program exists
	 */
	bool exit;

	// a bool for pausing the game
	bool pause;

	/**
	 * verbosity level.
	 * the higher, the more output has to be generated
	 */
	int verbosity_level;

	/**
	 * image processing filter id
	 */
	int stage_imageprocessing_filter_id;

	/**
	 * image processing threshold value
	 */
	int stage_imageprocessing_threshold_value;

	/**
	 * image processing: output flagfield
	 */
	bool stage_imageprocessing_output_flagfield;

	/**
	 * image processing: freeze output
	 */
	bool stage_imageprocessing_freeze_output;

	/**
	 * visualize flag field for simulation
	 */
	bool stage_fluidsimulation_visualize_flagfield;

	/**
	 * for multithreading: how many threads are used for the pipelining
	 */
	int threading_number_of_threads_to_use;

	/**
	 * stage image input: default image to load
	 */
	std::string stage_image_input_path;

	/**
	 * pipeline id to process
	 */
	int pipeline_id;

	/**
	 * inflow factor
	 */
	float lbm_inflow_factor;

	/**
	 * outflow factor
	 */
	float lbm_outflow_factor;

	/**
	 * string to video device being used for, e.g., webcam input
	 */
	std::string video_device;

	/**
	 * request particular video width
	 */
	int input_video_width;
	/**
	 * request particular video height
	 */
	int input_video_height;

	bool mute_music;

	float rb_acceleration_scaling;
	float rb_velocity_scaling;
	int f_rep_mode;
	int initialisation_mode;
	int max_number_rb;

	/**
	 * return bool if processed
	 */
	bool key_down(char i_key) {
		switch (i_key) {
		case SDLK_j:
			stage_imageprocessing_filter_id--;
			std::cout << "Using filter id " << stage_imageprocessing_filter_id
					<< std::endl;
			return true;

		case SDLK_k:
			stage_imageprocessing_filter_id++;
			std::cout << "Using filter id " << stage_imageprocessing_filter_id
					<< std::endl;
			return true;

		case SDLK_f:
			stage_imageprocessing_freeze_output =
					!stage_imageprocessing_freeze_output;
			std::cout << "Image input frozen: "
					<< stage_imageprocessing_freeze_output << std::endl;
			return true;

		case SDLK_g:
			stage_imageprocessing_threshold_value--;
			std::cout << "Using filter threshold value "
					<< stage_imageprocessing_threshold_value << std::endl;
			return true;

		case SDLK_t:
			stage_imageprocessing_threshold_value++;
			std::cout << "Using filter threshold value "
					<< stage_imageprocessing_threshold_value << std::endl;
			return true;

		case SDLK_v:
			stage_fluidsimulation_visualize_flagfield =
					!stage_fluidsimulation_visualize_flagfield;
			std::cout << "Visualize flag field: "
					<< stage_fluidsimulation_visualize_flagfield << std::endl;
			return true;

		case SDLK_e:
			rb_acceleration_scaling += 10;
			std::cout << "Acceleration Scaling of RBs is now "
					<< rb_acceleration_scaling << std::endl;
			return true;

		case SDLK_w:
			if(rb_acceleration_scaling -10 >= 0){
				rb_acceleration_scaling -= 10;
				std::cout << "Acceleration Scaling of RBs is now "
						<< rb_acceleration_scaling << std::endl;
			}else{
				std::cout << "Acceleration Scaling of RBs is already at its lowest" << std::endl;
			}
			return true;

		case SDLK_s:

			rb_velocity_scaling += 0.1;
			std::cout << "velocity scaling of RBs is now "
					<< rb_velocity_scaling << std::endl;
			return true;

		case SDLK_a:
			if(rb_velocity_scaling - 0.1 >= 0.0){
				rb_velocity_scaling -= 0.1;
				std::cout << "velocity scaling of RBs is now "
						<< rb_velocity_scaling << std::endl;
			}else{
				std::cout << "velocity scaling of RBs is already at its lowest" << std::endl;
			}
			return true;

		case SDLK_x:

			max_number_rb  += 1;
			std::cout << "maximum number of RBs is now "
					<< max_number_rb << std::endl;
			return true;

		case SDLK_y:
			if(max_number_rb - 1 >= 0.0){
				max_number_rb -= 1;
				std::cout << "maximum number of RBs is now "
						<< max_number_rb << std::endl;
			}else{
				std::cout << "maximum number of RBs is already at its lowest" << std::endl;
			}
			return true;

		default:
			// unknown key
			break;
		}
		return false;
	}

	CParameters() {
		/*
		 * setup default parameter values
		 */
		exit = false;
		pause = false;
		verbosity_level = -1;
		threading_number_of_threads_to_use = -1;
		stage_image_input_path = "data/fa_image_2012.jpg";

		stage_imageprocessing_filter_id = 0;
		stage_imageprocessing_threshold_value = 128;
		stage_imageprocessing_output_flagfield = false;
		stage_fluidsimulation_visualize_flagfield = false;

		lbm_inflow_factor = 1.2f;
		lbm_outflow_factor = 0.8f;

		pipeline_id = 0;

		video_device = "/dev/video0";
		input_video_width = -1;
		input_video_height = -1;

		mute_music = false;

		rb_acceleration_scaling = 0.0;
		rb_velocity_scaling = 0.0;
		f_rep_mode = 1;
		initialisation_mode = 1;
		max_number_rb = 5;
	}

	/**
	 * this function prints the program parameters to stdout
	 */
	void printMainParameterInfo(int /*i_argc*/,		/// number of arguments
			char * const i_argv[]/// array of strings of arguments preprocessed by prefixed binary
			) {
		std::cout << "Parameters for " << i_argv[0] << std::endl;
		std::cout << "	-d [device]			Use video device (default: /dev/video0)"
				<< std::endl;
		std::cout
				<< "	-i [input image]	Input image (default: ../data/fa_image_2012.jpg)"
				<< std::endl;
		std::cout << "	-v [level]			Verbosity level (default: -1)" << std::endl;
		std::cout << "	-p [pipeline id]	Pipeline id (default: 0)" << std::endl;
		std::cout << "	-w [video width]	(default: -1 [automatic])" << std::endl;
		std::cout << "	-h [video height]	(default: -1 [automatic])"
				<< std::endl;
		std::cout << "  -s [inflow factor]  (default: 1.2)" << std::endl;
		std::cout << "  -o [outflow factor] (default: 0.8)" << std::endl;
		::exit(EXIT_FAILURE);
	}

	/**
	 * setup the parameters given to the program.
	 *
	 * if any (optimistic!) of those parameters are invalid,
	 * noexit_printMainParameterInfo is called.
	 */
public:
	void setup(int i_argc,		/// number of arguments
			char * const i_argv[]/// array of strings of arguments preprocessed by prefixed binary
			) {
		int optchar;
		const char *options = "v:d:i:n:p:w:h:o:s:c:b:a:x:m:";

		while ((optchar = getopt(i_argc, i_argv, options)) > 0) {
			switch (optchar) {
			/*
			 * verbosity level
			 */
			case 'v':
				verbosity_level = atoi(optarg);
				break;

				/*
				 * input video device
				 */
			case 'd':
				video_device = optarg;
				break;

				/*
				 * input image
				 */
			case 'i':
				stage_image_input_path = optarg;
				break;

				/*
				 * number of threads to use
				 */
			case 'n':
				threading_number_of_threads_to_use = atoi(optarg);
				break;

				/*
				 * LBM inflow/source factor
				 */
			case 's':
				lbm_inflow_factor = atof(optarg);
				break;

				/*
				 * LBM outflow/sink factor
				 */
			case 'o':
				lbm_outflow_factor = atof(optarg);
				break;

				/*
				 * pipeline id
				 */
			case 'p':
				pipeline_id = atoi(optarg);
				break;

			case 'w':
				input_video_width = atoi(optarg);
				break;
			case 'h':
				input_video_height = atoi(optarg);
				break;
			case 'c':
				initialisation_mode = (atoi(optarg) >= 0 && atoi(optarg) < 3) ? atoi(optarg) : 1;
				break;
			case 'b':
				f_rep_mode = (atoi(optarg) >= 0 && atoi(optarg) < 3) ? atoi(optarg) : 1;
				break;
			case 'a':
				rb_acceleration_scaling = atof(optarg) >= 0.0 ? atof(optarg) : 0.0;
				break;
			case 'x':
				rb_velocity_scaling = atof(optarg) >= 0.0 ? atof(optarg) : 0.0;
				break;
			case 'm':
				max_number_rb = atof(optarg) >= 0 ? atof(optarg) : 0;
				break;

			default:
				printMainParameterInfo(i_argc, i_argv);
				break;
			}
		}
	}

};

#endif /* CPARAMETERS_HPP_ */
