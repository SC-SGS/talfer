//============================================================================
// Name        : fa_2013.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
//============================================================================

#include <iostream>
#ifdef PARALLEL
#include <tbb/task_scheduler_init.h>
#endif

//#include "XnOpenNI.h"

#ifdef KINECT_INCLUDED
#include "CStage_Kinect.hpp"
#else
#include "CStage_KeyboardForces.hpp"
#endif

#include <stdio.h>
#include <time.h>
//#include "CVTKTools.hpp"
#include "CPipelineStage.hpp"
#include "CParameters.hpp"
#include "CStage_ImageInput.hpp"
#include "CStage_ImageProcessing.hpp"
#include "CStage_TrivialImageProcessing.hpp"
#if defined(NAVIER_STOKES)
#include "CStage_FluidSimulationNS.hpp"
#define CStage_FluidSimulation CStage_FluidSimulationNS
#else
#ifdef PARALLEL
#include "CStage_FluidSimulationLBM_TBB.hpp"
#else
#include "CStage_FluidSimulationLBM.hpp"
#endif
#define CStage_FluidSimulation CStage_FluidSimulationLBM
#endif
#ifdef __APPLE__
//#include "CStage_VideoInput_Apple.hpp"
#else
#include "CStage_VideoInput.hpp"
#endif
#include "CStage_VideoOutput.hpp"
#include "CDataParticleArray.hpp"
#include "CStage_RigidBodySimulation.hpp"
#include "CStage_RigidBodySimulation_Coupling.hpp"
#include "CStage_StaticVelocityField.hpp"
#include "CStage_DynamicVelocityField.hpp"
#include "CStage_StaticRigidBodies.hpp"
#include "CStage_Visualization.hpp"
#if defined(GAME2)
#include "Game_version2.hpp"
#else
#include "Game.hpp"
#endif

/*
 * program parameters
 */
CParameters cParameters;

/*
 * Game states and actions
 */
Game* game;

bool written;

/*
 * main_visualization
 *
 * This is a main loop used for visualization testing.
 * Pre-Visualization-Stages are replaced by static Mock-Ups.
 */
void main_visualization() {
	/*
	 * Initialize the random number generator.
	 * We use random numbers for the positions of rigid bodies and for the
	 * origin points of stream lines.
	 */
	srand(time(NULL));

	/*
	 * Velocity Field Mock Up
	 */
	CStage_StaticVelocityField cStage_StaticVelocityField(cParameters);

	/*
	 * Rigid Body Mock Up
	 */
	CStage_StaticRigidBodies cStage_StaticRigidBodies(cParameters);

	/*
	 * Visualization Stage
	 */
	CStage_Visualization cStage_Visualization(cParameters);

	/*
	 * Video Output
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * Connections
	 *
	 * Scheme:
	 * StaticVelocityField \
     * StaticRigidBodies   - Visualization - VideoOutput
	 */
	cStage_StaticVelocityField.connectOutput(cStage_Visualization);
	cStage_StaticRigidBodies.connectOutput(cStage_Visualization);
	cStage_Visualization.connectOutput(cStage_VideoOutput);

	/*
	 * Manual initial push for static stages
	 */
	cStage_StaticVelocityField.pipeline_push();
	cStage_StaticRigidBodies.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);
		cStage_VideoOutput.main_loop_callback();
		cStage_Visualization.main_loop_callback();
	}
}

void main_image_viewer() {
	/*
	 * image input: read input file from a file
	 */
	CStage_ImageInput cStage_ImageInput(cParameters);

	/*
	 * video output: repeat output of an image
	 */
	//	CStage_VideoOutput cStage_VideoOutput(cParameters);
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * here starts the connection between different pipeline stages
	 */

	/*
	 * connect video input directly to video output
	 */
	cStage_ImageInput.connectOutput(cStage_VideoOutput);

	// forward image data once
	cStage_ImageInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);
		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();
	}
}

void main_image_modifier() {
	/*
	 * image input: read input file from a file
	 */
	CStage_ImageInput cStage_ImageInput(cParameters);

	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * here starts the connection between different pipeline stages
	 */

	/*
	 * connect video input directly to video output
	 */
	cStage_ImageInput.connectOutput(cStage_VideoOutput);

	/*
	 * also connect video output (mouse & keyboard) to image input
	 */
	cStage_VideoOutput.connectOutput(cStage_ImageInput);

	/*
	 * forward image data once to have something for output
	 */
	cStage_ImageInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);
		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();
	}
}

void main_image_modifier_and_filter() {
	/*
	 * image input: read input file from a file
	 */
	CStage_ImageInput cStage_ImageInput(cParameters);

	/*
	 * image input: read input file from a file
	 */
	CStage_ImageProcessing cStage_ImageProcessing(cParameters);
	cParameters.stage_imageprocessing_filter_id = 5;
	cParameters.stage_imageprocessing_threshold_value = 120;

	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * here starts the connection between different pipeline stages
	 */

	/*
	 * image input -> image processing
	 */
	cStage_ImageInput.connectOutput(cStage_ImageProcessing);

	/*
	 * image processing -> video output
	 */
	cStage_ImageProcessing.connectOutput(cStage_VideoOutput);

	/*
	 * video output -> image input modifications
	 */
	cStage_VideoOutput.connectOutput(cStage_ImageInput);

	// forward image data once to have something for output
	cStage_ImageInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);
		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();
		//cParameters.exit = true;
	}
}

#if !defined(__APPLE__)
void main_video_viewer() {
	/*
	 * image input: read input file from a file
	 */
	CStage_VideoInput cStage_VideoInput(cParameters);

	/*
	 * image input: read input file from a file
	 */
	CStage_ImageProcessing cStage_ImageProcessing(cParameters);
	cParameters.stage_imageprocessing_filter_id = 3;	// use threshold
	cParameters.stage_imageprocessing_threshold_value = 120; // dynamically adapted

	/*
	 * video output: repeat output of an image
	 */
	//	CStage_VideoOutput cStage_VideoOutput(cParameters);
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * here starts the connection between different pipeline stages
	 */

	cStage_VideoInput.connectOutput(cStage_ImageProcessing);
	cStage_ImageProcessing.connectOutput(cStage_VideoOutput);
	cStage_VideoInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);

		// test for video input
		cStage_VideoInput.main_loop_callback();

		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();
	}
}
#endif

void main_sim_static_image() {
	/*
	 * image input: read input file from a file
	 */
	CStage_ImageInput cStage_ImageInput(cParameters);

	/*
	 * rigid body simulation: 
	 */
	CStage_RigidBodySimulation cStage_RigidBodySimulation(cParameters);

	/*
	 * image input: read input file from a file
	 */
	CStage_TrivialImageProcessing cStage_ImageProcessing(cParameters);
	cParameters.stage_imageprocessing_output_flagfield = true;

//	cParameters.stage_imageprocessing_filter_id = 5;	// use threshold
	cParameters.stage_imageprocessing_filter_id = -1; // no filter

	cParameters.stage_imageprocessing_threshold_value = 120;

	cParameters.stage_fluidsimulation_visualize_flagfield = true;

	/*
	 * Video processing to rigid body simulation
	 */
	cStage_ImageProcessing.connectOutput(cStage_RigidBodySimulation);

	/*
	 * fluid simulation test class
	 */
	CStage_FluidSimulation cStage_FluidSimulationLBM(cParameters);

	CStage_Visualization cStage_Visualization(cParameters);
	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * kinect: read gesture from kinect and give out
	 * vector of hand position
	 */
#ifdef KINECT_INCLUDED
	CStage_Kinect cStage_Kinect(cParameters);
#else
	CStage_KeyboardForces cStage_KeyboardForces(cParameters);
	/*
	 * Connect KeyboardForces to VideoOutput
	 */
	cStage_VideoOutput.setupCStageKeyboardForces(&cStage_KeyboardForces);
#endif
	/*
	 * here starts the connection between different pipeline stages
	 */
	cStage_ImageInput.connectOutput(cStage_ImageProcessing);
	cStage_ImageProcessing.connectOutput(cStage_FluidSimulationLBM);
	cStage_ImageProcessing.connectOutput(cStage_Visualization);
	cStage_FluidSimulationLBM.connectOutput(cStage_Visualization);
	cStage_Visualization.connectOutput(cStage_VideoOutput);
	cStage_RigidBodySimulation.connectOutput(cStage_Visualization);
	cStage_RigidBodySimulation.connectOutput(cStage_FluidSimulationLBM);
	cStage_FluidSimulationLBM.connectOutput(cStage_RigidBodySimulation);
	cStage_VideoOutput.connectOutput(cStage_ImageInput);
#ifdef KINECT_INCLUDED
	// kinect to rigid body and visualization
	cStage_Kinect.connectOutput(cStage_RigidBodySimulation);
	cStage_Kinect.connectOutput(cStage_Visualization);
#else
	cStage_KeyboardForces.connectOutput(cStage_RigidBodySimulation);
	cStage_KeyboardForces.connectOutput(cStage_Visualization);
#endif

	/*
	 * forward image data once to have something for output
	 */
	cStage_ImageInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	game->setCParameters(&cParameters);
	cStage_VideoOutput.connectVisualization(&cStage_Visualization);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);

		// EDIT
		game->updateVideoOutput();
		game->check_if_paused();
		// EDIT

#ifdef KINECT_INCLUDED
		// read kinect input and generate vector
#ifdef PARALLEL
		cStage_Kinect.main_loop_callback_threaded();
#else
		cStage_Kinect.main_loop_callback();
#endif

#endif
		if (!game->global_pause && !game->has_ended) {
#ifndef PARALLEL
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback();
#else
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback_threaded();
#endif
		}

		// visualize
		cStage_Visualization.main_loop_callback();

		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();

		if (!game->global_pause && !game->has_ended) {
			// simulate rigid body
			cStage_RigidBodySimulation.main_loop_callback();
		}
	}
}

void main_lbm_rbs_coupling() {
	/*
	 * image input: read input file from a file
	 */
	CStage_ImageInput cStage_ImageInput(cParameters);

	/*
	 * rigid body simulation:
	 */
	CStage_RigidBodySimulation_Coupling cStage_RigidBodySimulation(cParameters);

	/*
	 * image input: read input file from a file
	 */
	CStage_TrivialImageProcessing cStage_ImageProcessing(cParameters);
	cParameters.stage_imageprocessing_output_flagfield = true;

//	cParameters.stage_imageprocessing_filter_id = 5;	// use threshold
	cParameters.stage_imageprocessing_filter_id = -1; // no filter

	cParameters.stage_imageprocessing_threshold_value = 120;

	cParameters.stage_fluidsimulation_visualize_flagfield = true;

	/*
	 * Video processing to rigid body simulation
	 */
	cStage_ImageProcessing.connectOutput(cStage_RigidBodySimulation);

	/*
	 * fluid simulation test class
	 */
	CStage_FluidSimulation cStage_FluidSimulationLBM(cParameters);

	CStage_Visualization cStage_Visualization(cParameters);
	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * kinect: read gesture from kinect and give out
	 * vector of hand position
	 */
#ifdef KINECT_INCLUDED
	CStage_Kinect cStage_Kinect(cParameters);
#else
	CStage_KeyboardForces cStage_KeyboardForces(cParameters);
	/*
	 * Connect KeyboardForces to VideoOutput
	 */
	cStage_VideoOutput.setupCStageKeyboardForces(&cStage_KeyboardForces);
#endif
	/*
	 * here starts the connection between different pipeline stages
	 */
	cStage_ImageInput.connectOutput(cStage_ImageProcessing);
	cStage_ImageProcessing.connectOutput(cStage_FluidSimulationLBM);
	cStage_ImageProcessing.connectOutput(cStage_Visualization);
	cStage_FluidSimulationLBM.connectOutput(cStage_Visualization);
	cStage_Visualization.connectOutput(cStage_VideoOutput);
	cStage_RigidBodySimulation.connectOutput(cStage_Visualization);
	cStage_RigidBodySimulation.connectOutput(cStage_FluidSimulationLBM);
	cStage_FluidSimulationLBM.connectOutput(cStage_RigidBodySimulation);
	cStage_VideoOutput.connectOutput(cStage_ImageInput);
#ifdef KINECT_INCLUDED
	// kinect to rigid body and visualization
	cStage_Kinect.connectOutput(cStage_RigidBodySimulation);
	cStage_Kinect.connectOutput(cStage_Visualization);
#else
	cStage_KeyboardForces.connectOutput(cStage_RigidBodySimulation);
	cStage_KeyboardForces.connectOutput(cStage_Visualization);
#endif

	/*
	 * forward image data once to have something for output
	 */
	cStage_ImageInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	game->setCParameters(&cParameters);
	cStage_VideoOutput.connectVisualization(&cStage_Visualization);
	//int t_cnt = 0;
	//int step = 0;
	//float t_ave = 0.0;
	while (!cParameters.exit) {
		//struct timespec tps, tpe;
		//clock_gettime(CLOCK_REALTIME, &tps);

		game->mute_music(cParameters.mute_music);

		// EDIT
		game->updateVideoOutput();
		game->check_if_paused();

		// EDIT

#ifdef KINECT_INCLUDED
		// read kinect input and generate vector
#ifdef PARALLEL
		cStage_Kinect.main_loop_callback_threaded();
#else
		cStage_Kinect.main_loop_callback();
#endif

#endif
		if (!game->global_pause && !game->has_ended) {
#ifndef PARALLEL
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback();
#else
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback_threaded();
#endif
		}

		// visualize
		cStage_Visualization.main_loop_callback();

		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();

		if (!game->global_pause && !game->has_ended) {
			// simulate rigid body
			cStage_RigidBodySimulation.main_loop_callback();
		}
		/*
		clock_gettime(CLOCK_REALTIME, &tpe);
		float ms_time = (float)(tpe.tv_nsec-tps.tv_nsec) / 1000000.0;
		t_ave += ms_time > 0 ? ms_time : 0.0;
		t_cnt += ms_time > 0 ? 1 : 0;
		step++;
		//printf("FPS = %f with Frametime = %f ns \n", 1000.0/ms_time, ms_time);
		if(t_cnt == 100){
			std::vector<float> t_data = std::vector<float>();
			t_data.push_back(t_ave / (float)t_cnt);
			t_data.push_back(1000.0 / (t_ave / (float)t_cnt));
			writeData("timings", step, t_data);
			printf("average frametime = %f ms, average FPS = %f\n", t_ave / (float)t_cnt, 1000.0 / (t_ave / (float)t_cnt));
			t_ave = 0.0;
			t_cnt = 0;
		}*/
	}
}

#if !defined(__APPLE__)
void main_sim_video() {
	/*
	 * image input: read input file from a file
	 */
	CStage_VideoInput cStage_VideoInput(cParameters);

	/*
	 * rigid body simulation: 
	 */
	CStage_RigidBodySimulation cStage_RigidBodySimulation(cParameters);

	CStage_ImageProcessing cStage_ImageProcessing(cParameters);
	cParameters.stage_imageprocessing_output_flagfield = true;
	cParameters.stage_imageprocessing_filter_id = 5;
	cParameters.stage_imageprocessing_threshold_value = 120;

	/*
	 * fluid simulation test class
	 */
	CStage_FluidSimulation cStage_FluidSimulationLBM(cParameters);

	CStage_Visualization cStage_Visualization(cParameters);
	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);
	/*
	 * kinect: read gesture from kinect and give out
	 * vector of hand position
	 */
#ifdef KINECT_INCLUDED
	CStage_Kinect cStage_Kinect(cParameters);
#else
	CStage_KeyboardForces cStage_KeyboardForces(cParameters);
	cStage_VideoOutput.setupCStageKeyboardForces(&cStage_KeyboardForces);
#endif

	/*
	 * here starts the connection between different pipeline stages
	 */
	cStage_ImageProcessing.connectOutput(cStage_RigidBodySimulation);
	cStage_VideoInput.connectOutput(cStage_ImageProcessing);
	cStage_ImageProcessing.connectOutput(cStage_FluidSimulationLBM);
	cStage_ImageProcessing.connectOutput(cStage_Visualization);
	cStage_FluidSimulationLBM.connectOutput(cStage_Visualization);
	cStage_Visualization.connectOutput(cStage_VideoOutput);
	cStage_RigidBodySimulation.connectOutput(cStage_Visualization);
	cStage_FluidSimulationLBM.connectOutput(cStage_RigidBodySimulation);
#ifdef KINECT_INCLUDED
	// kinect to rigid body and visualization
	cStage_Kinect.connectOutput(cStage_RigidBodySimulation);
	cStage_Kinect.connectOutput(cStage_Visualization);
#else
	cStage_KeyboardForces.connectOutput(cStage_RigidBodySimulation);
	cStage_KeyboardForces.connectOutput(cStage_Visualization);
#endif
	/*
	 * forward image data once to have something for output
	 */
	cStage_VideoInput.pipeline_push();

	/*
	 * forward image data once to have something for output
	 */
	cStage_VideoInput.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);

#ifdef KINECT_INCLUDED
		// read kinect input and generate vector
		cStage_Kinect.main_loop_callback();
#endif
		// trigger image input to do something

		// test for video input

		cStage_VideoInput.main_loop_callback();

		if (!game->global_pause && !game->has_ended) {
#ifndef PARALLEL
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback();
#else
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback_threaded();
#endif
		}

		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();

		if (!game->global_pause && !game->has_ended) {
			// simulate rigid body
			cStage_RigidBodySimulation.main_loop_callback();
		}

		cStage_Visualization.main_loop_callback();
	}
}
#endif

#if !defined(__APPLE__)
void main_lbm_with_visualization_video() {
	/*
	 * image input: read input file from a file
	 */
	CStage_VideoInput cStage_VideoInput(cParameters);

	/*
	 * image input: read input file from a file
	 */
	CStage_ImageProcessing cStage_ImageProcessing(cParameters);
	cParameters.stage_imageprocessing_output_flagfield = true;
	cParameters.stage_fluidsimulation_visualize_flagfield = true;
	cParameters.stage_imageprocessing_threshold_value = 64;

	/*
	 * fluid simulation test class
	 */
	CStage_FluidSimulation cStage_FluidSimulationLBM(cParameters);

	/*
	 * Rigid Body Mock Up
	 */
	//CStage_StaticRigidBodies cStage_StaticRigidBodies(cParameters);
	/*
	 * Visualization Stage
	 */
	CStage_Visualization cStage_Visualization(cParameters);

	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * here starts the connection between different pipeline stages
	 */
	cStage_VideoInput.connectOutput(cStage_ImageProcessing);
	cStage_ImageProcessing.connectOutput(cStage_FluidSimulationLBM);
	//cStage_StaticRigidBodies.connectOutput(cStage_Visualization);
	cStage_FluidSimulationLBM.connectOutput(cStage_Visualization);
	cStage_Visualization.connectOutput(cStage_VideoOutput);

	//cStage_StaticRigidBodies.pipeline_push();

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);

		// trigger image input to do something
		cStage_VideoInput.main_loop_callback();

		if (!game->global_pause && !game->has_ended) {
#ifndef PARALLEL
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback();
#else
			// do simulation timestep
			cStage_FluidSimulationLBM.main_loop_callback_threaded();
#endif
		}

		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();

		cStage_Visualization.main_loop_callback();
	}
}
#endif

void dynamic_field_without_simulation() {

	/*
	 * Dynamic Field Computation Stage
	 */
	CStage_DynamicVelocityField cStage_DynamicVelocityField(cParameters);

	/*
	 * Visualization Stage
	 */
	CStage_Visualization cStage_Visualization(cParameters);

	/*
	 * video output: repeat output of an image
	 */
	CStage_VideoOutput cStage_VideoOutput(cParameters);

	/*
	 * here starts the connection between different pipeline stages
	 */
	cStage_DynamicVelocityField.connectOutput(cStage_Visualization);
	cStage_Visualization.connectOutput(cStage_VideoOutput);

	/*
	 * main processing queue
	 */
	game = Game::getInstance();
	game->setVideoOutput(&cStage_VideoOutput);
	while (!cParameters.exit) {
		game->mute_music(cParameters.mute_music);

		cStage_DynamicVelocityField.main_loop_callback();

		// trigger image input to do something
		cStage_VideoOutput.main_loop_callback();

		cStage_Visualization.main_loop_callback();
	}

}

/**
 * Main C entry function
 */
int main(int argc,		/// number of arguments
		char *argv[]/// array of strings of arguments preprocessed by prefixed binary
		) {

	written = false;
	/*
	 * setup program parameters
	 */
	cParameters.setup(argc, argv);

#ifdef PARALLEL
	tbb::task_scheduler_init init();
#endif

	printf("start \n");

	switch (cParameters.pipeline_id) {
	case 0:
		main_image_viewer();
		break;

	case 1:
		main_image_modifier();
		break;

	case 2:
		main_image_modifier_and_filter();
		break;

	case 3:
#if !defined(__APPLE__)
		main_video_viewer();
#endif
		break;

	case 4:
		main_sim_static_image();
		break;

	case 5:
#if !defined(__APPLE__)
		main_sim_video();
#endif
		break;

	case 6:
		main_visualization();
		break;

	case 7:
#if !defined(__APPLE__)
		main_lbm_with_visualization_video();
#endif
		break;
	case 8:
		main_lbm_rbs_coupling();
		break;
	}

	return 0;
}
