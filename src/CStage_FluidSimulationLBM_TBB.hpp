/*
 * CStage_ImageInput.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CSTAGE_FLUIDSIMULATION_LBM_HPP_
#define CSTAGE_FLUIDSIMULATION_LBM_HPP_

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "SDL/SDL_image.h"
#include "CDataArray2D.hpp"
#include "CDataDrawingInformation.hpp"
#include <math.h>
#include "tbb/blocked_range2d.h"
#include "tbb/parallel_for.h"

// D2Q9 Directions
#define C 0
#define E 1
#define N 2
#define W 3
#define S 4
#define NE 5
#define NW 6
#define SW 7
#define SE 8
#define X 0
#define Y 1

/*
 *  D2Q9 - BGK model parameters
 */
// Relaxation parameter
const float omega = 1.5; // = 1.0/tau; // Usually between 0.5 and 1.95
// D2Q9 Weights
// C, E, N, W, S, NE, NW, SW, SE
const float weights[9] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0
		/ 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
const int units[9][2] = { { 0, 0 }, { +1, 0 }, { 0, +1 }, { -1, 0 }, { 0, -1 },
		{ +1, +1 }, { -1, +1 }, { -1, -1 }, { +1, -1 } };
const int inverse[9] = { C, W, S, E, N, SW, SE, NE, NW };

/**
 * class providing static image input
 *
 * the image is send to the pipeline during each main loop iteration
 */
class CStage_FluidSimulationLBM: public CPipelineStage {
	/**
	 * global parameters
	 */
	CParameters &cParameters;

	/**
	 * input flag image
	 */
	CDataArray2D<unsigned char, 1> input_cDataArray2D;

	/**
	 * processed velocity field
	 */
	CDataArray2D<float, 2> output_velocity;

	/**
	 * processed pressure field
	 */
	CDataArray2D<float, 1> output_pressure;

	/**
	 * simulation array to run simulation on
	 */
	CDataArray2D<float, 9> simulationData_BufferSrc;
	CDataArray2D<float, 9> simulationData_BufferDst;

	int width, height;

public:
	/**
	 * constructor
	 */
	CStage_FluidSimulationLBM(CParameters &i_cParameters) :
			CPipelineStage("FluidSimulationLBM"), cParameters(i_cParameters), width(
					0), height(0) {
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push() {
//		if (cParameters.stage_fluidsimulation_visualize_flagfield)
//			CPipelineStage::pipeline_push((CPipelinePacket&)input_cDataArray2D);

		CPipelineStage::pipeline_push((CPipelinePacket&) output_velocity);
		// CPipelineStage::pipeline_push((CPipelinePacket&)output_pressure);
	}

public:
	void simulation_timestep() {
		bool updated_input = false;

		if (!input_cDataArray2D.isValidData())
			return;

		if (width != input_cDataArray2D.width
				|| height != input_cDataArray2D.height) {
			width = input_cDataArray2D.width;
			height = input_cDataArray2D.height;

			simulationData_BufferSrc.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);
			simulationData_BufferDst.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);

			output_velocity.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);
			output_pressure.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);

			// reset everything to equilibrium
			for (int y = 0; y < output_velocity.height; y++) {
				for (int x = 0; x < output_velocity.width; x++) {
					for (int i = 0; i < 9; i++) {
						simulationData_BufferSrc.getRef(x, y, i) = weights[i];
						simulationData_BufferDst.getRef(x, y, i) = weights[i];
					}
				}
			}

			updated_input = true;
		}

		if (updated_input || input_cDataArray2D.has_changed) {
			/* 
			 * process input/output/boundary flags
			 */
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					unsigned char flag = input_cDataArray2D.getRef(x, y);
					float *f_src = &simulationData_BufferSrc.getRef(x, y);
					float *f_dst = &simulationData_BufferDst.getRef(x, y);
					// float *f;

					switch (flag) {
					case 1: // obstacle
						for (int i = 0; i < 9; i++) {
							f_src[i] = weights[i];
							f_dst[i] = weights[i];
						}
						break;

					case 2: // inflow/source
						for (int i = 0; i < 9; i++) {
							f_src[i] = weights[i]
									* cParameters.lbm_inflow_factor;
							f_dst[i] = weights[i]
									* cParameters.lbm_inflow_factor;
						}
						break;

					case 3: // outflow/drain
						for (int i = 0; i < 9; i++) {
							f_src[i] = weights[i]
									* cParameters.lbm_outflow_factor;
							f_dst[i] = weights[i]
									* cParameters.lbm_outflow_factor;
						}
						break;
					}
				}
			}
		}

		tbb::parallel_for(tbb::blocked_range2d<int, int>(0, height, 0, width),
				CollisionStream(height, width, simulationData_BufferSrc,
						simulationData_BufferDst, output_velocity,
						output_pressure, input_cDataArray2D)); //cs);
		simulationData_BufferSrc.swap(simulationData_BufferDst);

		pipeline_push();
	}

public:
	void main_loop_callback() {
		simulation_timestep();
		simulation_timestep();
	}

public:
	/**
	 * process incoming pipeline input.
	 *
	 * the only input processed so far is from the video output stage to
	 * draw something into the image.
	 */
	void pipeline_process_input(CPipelinePacket &i_cPipelinePacket) {
		// we are currently only able to process "unsigned char,1" flag arrays.
		if (i_cPipelinePacket.type_info_name
				!= typeid(CDataArray2D<unsigned char, 1> ).name()) {
			std::cerr
					<< "ERROR: Simulation input is only able to process (char,1) flag arrays"
					<< std::endl;
			exit(-1);
		}

		// unpack data
		CDataArray2D<unsigned char, 1> *input = i_cPipelinePacket.getPayload<
				CDataArray2D<unsigned char, 1> >();

		// copy data to input array
		input_cDataArray2D.resize(input->width, input->height);
		input_cDataArray2D.loadData(input->data);
	}
private:
	class CollisionStream {
	public:
		CollisionStream(int height, int width,
				CDataArray2D<float, 9>& simulationData_BufferDst,
				CDataArray2D<float, 9> simulationData_BufferSrc,
				CDataArray2D<float, 2>& output_velocity,
				CDataArray2D<float, 1>& output_pressure,
				CDataArray2D<unsigned char, 1> input_cDataArray2D) :
				height(height), width(width), simulationData_BufferSrc(
						&simulationData_BufferSrc), simulationData_BufferDst(
						&simulationData_BufferDst), input_cDataArray2D(
						&input_cDataArray2D), output_velocity(&output_velocity), output_pressure(
						&output_pressure) {
		}

		const int height, width;
		const CDataArray2D<float, 9>* simulationData_BufferSrc;
		CDataArray2D<float, 9>* simulationData_BufferDst;
		const CDataArray2D<unsigned char, 1>* input_cDataArray2D;
		CDataArray2D<float, 2>* output_velocity;
		CDataArray2D<float, 1>* output_pressure;

		void operator()(const tbb::blocked_range2d<int, int>& br) const {
			for (int y_src = br.rows().begin(); y_src != br.rows().end();
					++y_src) {
				for (int x_src = br.cols().begin(); x_src != br.cols().end();
						++x_src) {
					/*
					 *  Collision step (D2Q9 model)
					 */
					const unsigned char flag_src = input_cDataArray2D->getRef(
							x_src, y_src);
					float *f_src = &simulationData_BufferSrc->getRef(x_src,
							y_src);

					if (flag_src == 1) { // obstacle
						continue; // does not stream to other cells
					}

					float* new_f_src;
					float temp_f_src[9];

					// Does not need calculation of new probability distributions
					// from in- or outflow
					if (flag_src != 2 && flag_src != 3) {
						// Calculate macroscopic quantities:
						// density
						const float rho = f_src[NW] + f_src[N] + f_src[NE]
								+ f_src[W] + f_src[C] + f_src[E] + f_src[SW]
								+ f_src[S] + f_src[SE];
						const float rho_inv = 1 / rho;

						// velocity components
						const float u_x = (-f_src[NW] + 0.0 + f_src[NE]
								- f_src[W] + 0.0 + f_src[E] - f_src[SW] + 0.0
								+ f_src[SE]) * rho_inv;
						const float u_y = (f_src[NW] + f_src[N] + f_src[NE]
								+ 0.0 + 0.0 + 0.0 - f_src[SW] - f_src[S]
								- f_src[SE]) * rho_inv;
						const float u_x_sq = u_x * u_x;
						const float u_y_sq = u_y * u_y;
						const float u_sq = u_x_sq + u_y_sq;

						// Updating velocities for visulization
						// velocity vector
						output_velocity->getRef(x_src, y_src, 0) = u_x;
						output_velocity->getRef(x_src, y_src, 1) = u_y;
						// Pressure
						//output_pressure->getRef(x_src,y_src,0) = rho_inv;

						// equilibrium function
						const float feq[9] = { weights[C] * rho
								* (1.0f - 3.0f / 2.0f * u_sq), weights[E] * rho
								* (1.0f + 3.0f * u_x + 9.0f / 2.0f * u_x_sq
										- 3.0f / 2.0f * u_sq), weights[N] * rho
								* (1.0f + 3.0f * u_y + 9.0f / 2.0f * u_y_sq
										- 3.0f / 2.0f * u_sq), weights[W] * rho
								* (1.0f - 3.0f * u_x + 9.0f / 2.0f * u_x_sq
										- 3.0f / 2.0f * u_sq), weights[S] * rho
								* (1.0f - 3.0f * u_y + 9.0f / 2.0f * u_y_sq
										- 3.0f / 2.0f * u_sq), weights[NE] * rho
								* (1.0f + 3.0f * (u_x + u_y)
										+ 9.0f / 2.0f * (u_x + u_y)
												* (u_x + u_y)
										- 3.0f / 2.0f * u_sq), weights[NW] * rho
								* (1.0f + 3.0f * (-u_x + u_y)
										+ 9.0f / 2.0f * (-u_x + u_y)
												* (-u_x + u_y)
										- 3.0f / 2.0f * u_sq), weights[SW] * rho
								* (1.0f + 3.0f * (-u_x - u_y)
										+ 9.0f / 2.0f * (-u_x - u_y)
												* (-u_x - u_y)
										- 3.0f / 2.0f * u_sq), weights[SE] * rho
								* (1.0f + 3.0f * (u_x - u_y)
										+ 9.0f / 2.0f * (u_x - u_y)
												* (u_x - u_y)
										- 3.0f / 2.0f * u_sq) };

						for (int i = 0; i < 9; i++)
							temp_f_src[i] = f_src[i]
									- omega * (f_src[i] - feq[i]);

						new_f_src = temp_f_src;
					} else {
						new_f_src = f_src;
					}

					for (int direction_src = 0; direction_src < 9;
							direction_src++) {
						// Destination coordinates
						int x_dst = x_src + units[direction_src][X];
						int y_dst = y_src + units[direction_src][Y];

						// Periodic wrapping (faster than modulo)
						x_dst = std::min(std::max(x_dst, 0), width - 1);
						y_dst = std::min(std::max(y_dst, 0), height - 1);

						// Get flag at destination
						const unsigned char flag_dst =
								input_cDataArray2D->getRef(x_dst, y_dst);

						int direction_dst = direction_src;
						switch (flag_dst) {
						case 1: // to obstacle
							// relefction
							direction_dst = inverse[direction_src];
							x_dst = x_src;
							y_dst = y_src;
							break;
						case 2: // to inflow
						case 3: // to outflow
							// Ignore, this flow will "disappear"
							continue;
							break;
						}

						// Final destination
						float *f_dst = &simulationData_BufferDst->getRef(x_dst,
								y_dst);

						f_dst[direction_dst] = new_f_src[direction_src];
					}
				}
			}
		}
	};
};

#endif /* CSTAGE_FLUIDSIMULATION_LBM_HPP_ */
