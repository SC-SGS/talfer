/*
 * CStage_ImageInput.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CSTAGE_FLUIDSIMULATION_LBM_HPP_
#define CSTAGE_FLUIDSIMULATION_LBM_HPP_
#define NO_OMP

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "SDL/SDL_image.h"
#include "CDataArray2D.hpp"
#include "CDataDrawingInformation.hpp"
#include "CFluidSimulationData.hpp"
#include "CDataParticleArray.hpp"
//#include "CVTKTools.hpp"
#include <math.h>
#include <vector>
#include <locale>
#include <stdio.h>
#include <cstring>
//#include <codecvt>

//#define M_E 2.71828182845904523536
#define PI 3.14159265359
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
const float c_2i = 3.0;
const float c_2 = 1.0 / c_2i;

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
	 * input rigid bodies
	 */
	CDataParticleArray input_particles;
	std::vector<Vector2> output_particleForces;
	CDataArray2D<unsigned char, 1> particle_cDataArray2D;
	CDataArray2D<unsigned char, 1> particle_cDataArray2D_last;
	CDataArray2D<float, 1> density_BufferSrc;
	CDataArray2D<float, 1> density_BufferDst;
	CDataArray2D<float, 2> particle_Forces;
	const bool filter = false;
	const int filter_width = 100;
	const float filter_sigma = 50.0;
	bool particle_set;
	int timestep;
	float rb_v_max;
	const float fluid_rho_max = 2.0;
	const float fluid_v_max = 0.9;

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
			CPipelineStage("FluidSimulationLBM"), cParameters(i_cParameters), particle_set(
					false), timestep(0), width(0), height(0) {
//		rb_v_max = 0.005;
		rb_v_max = (cParameters.lbm_inflow_factor / 10.0) * 1.0;
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push() {
//		if (cParameters.stage_fluidsimulation_visualize_flagfield)
//			CPipelineStage::pipeline_push((CPipelinePacket&)input_cDataArray2D);
		CFluidSimulationData fsdata(output_particleForces, output_velocity);
		CPipelineStage::pipeline_push((CPipelinePacket&) fsdata);
		// CPipelineStage::pipeline_push((CPipelinePacket&)output_pressure);
	}

private:
	void updateRBOnDataArray() {
		if (particle_set && input_particles.data->size() > 0) {
			particle_cDataArray2D.swap(particle_cDataArray2D_last);
			for (int i = 0; i < width; ++i) {
				for (int j = 0; j < height; ++j) {
					unsigned char *flag = &particle_cDataArray2D.getRef(i, j);
					*flag = 0;
				}
			}
			std::vector<CParticle*> particles = *input_particles.data;
			for (unsigned int i = 0; i < particles.size(); ++i) {
				CParticle *p = particles.at(i);

				// specifies the search area with RB position and radius
				Vector2 pos = p->position;
				int rad = p->radius;
				int startx = round(pos.x) - rad - 1;
				if (startx < 0)
					startx = 0;
				int starty = round(pos.y) - rad - 1;
				if (starty < 0)
					starty = 0;
				int endx = round(pos.x) + rad + 2;
				if (endx > width)
					startx = width;
				int endy = round(pos.y) + rad + 2;
				if (endy > height)
					starty = height;

				// for all fluid nodes in the search area
				for (int x = startx; x < endx; ++x) {
					for (int y = starty; y < endy; ++y) {
						// if node is overlapped by RB -> set flag to RB
						if (pow(particles.at(i)->position.x - (float) x, 2.0)
								+ pow(particles.at(i)->position.y - (float) y,
										2.0)
								<= particles.at(i)->radius
										* particles.at(i)->radius) {
							//printf("x = (%i, %i) -> x² = %f, y² = %f, rad² = %i\n", x, y, pow(particles.at(i)->position.x - (float)x, 2.0), pow(particles.at(i)->position.y - (float)y, 2.0), particles.at(i)->radius * particles.at(i)->radius);
							unsigned char *flag = &particle_cDataArray2D.getRef(
									x, y);
							*flag = 6;
							float *f_src = &simulationData_BufferSrc.getRef(x,
									y);
							float *f_dst = &simulationData_BufferDst.getRef(x,
									y);
							// resets the fluid node data to default to avoid errors later on
							for (int i = 0; i < 9; i++) {
								f_src[i] = weights[i];
								f_dst[i] = weights[i];
							}
						}
					}
				}
			}

		}
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

			density_BufferSrc.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);
			density_BufferDst.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);

			particle_cDataArray2D.resize(width, height);
			particle_cDataArray2D_last.resize(width, height);
			for (int i = 0; i < width; ++i) {
				for (int j = 0; j < height; ++j) {
					particle_cDataArray2D.getRef(i, j) = 0;
					particle_cDataArray2D_last.getRef(i, j) = 0;
				}
			}

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
						//printf("height = %i, y = %i\n", height, y);
						for (int i = 0; i < 9; i++) {

							//Vector2 Vec_U_y = Vector2(U_y, 0.0f);
							//float rho_0 = 1.4f;
							//f_src[i] = equilDist(rho_0, Vec_U_y, i);
							//f_dst[i] = equilDist(rho_0, Vec_U_y, i);
							f_src[i] = weights[i]
									* cParameters.lbm_inflow_factor;
							f_dst[i] = weights[i]
									* cParameters.lbm_outflow_factor;
						}
						density_BufferDst.getRef(x, y) = 1.0f;
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

		updateRBOnDataArray();

		if (particle_set && input_particles.data->size() > 0) {
			output_particleForces = std::vector<Vector2>();
			// slots for the RB forces
			for (unsigned int i = 0; i < input_particles.data->size(); ++i) {
				output_particleForces.push_back(Vector2(0.0, 0.0));
			}
			// slots for the RB torques
			for (unsigned int i = 0; i < input_particles.data->size(); ++i) {
				output_particleForces.push_back(Vector2(0.0, 0.0));
			}
		}
		collision_stream();

		simulationData_BufferDst.swap(simulationData_BufferSrc);
		density_BufferDst.swap(density_BufferSrc);

		++timestep;

		pipeline_push();
	}

private:
	/**
	 * returns the Index of the RB a given RB point belongs to
	 */
	int getParticleIndex(int x, int y) {
		for (unsigned int i = 0; i < input_particles.data->size(); ++i) {
			if (pow(input_particles.data->at(i)->position.x - (float) x, 2.0)
					+ pow(input_particles.data->at(i)->position.y - (float) y,
							2.0)
					<= input_particles.data->at(i)->radius
							* input_particles.data->at(i)->radius) {
				return i;
			}
		}
		return -1;
//		printf("ERROR: No Particle matches flagposition (%i, %i)", x, y);
//		exit(-1);
	}

private:
	/**
	 * returns the nearest RB to a given fluid point
	 */
	CParticle* getNearestParticle(int x, int y) {
		float min_dist = width * height;
		int min_dist_index = -1;
		for (unsigned int i = 0; i < input_particles.data->size(); ++i) {
			if (pow(input_particles.data->at(i)->position.x - (float) x, 2.0)
					+ pow(input_particles.data->at(i)->position.y - (float) y,
							2.0)
					- (input_particles.data->at(i)->radius
							* input_particles.data->at(i)->radius) < min_dist) {
				min_dist = pow(
						input_particles.data->at(i)->position.x - (float) x,
						2.0)
						+ pow(
								input_particles.data->at(i)->position.y
										- (float) y, 2.0)
						- (input_particles.data->at(i)->radius
								* input_particles.data->at(i)->radius);
				min_dist_index = i;
			}
		}
		if(min_dist_index < 0){
//			printf("Löschvorgang");
			return NULL;
		}
		return input_particles.data->at(min_dist_index);
	}

private:
	/**
	 * computes the weights of a gaussian kernel
	 * (used for data output)
	 */
	void computeGaussianWeights(float* res, int kwidth, float sigma) {
		float a = 1.0 / (sigma * pow(2.0 * PI, 0.5));
		for (int i = -kwidth; i <= kwidth; ++i) {
			res[i + kwidth] = a
					* pow(M_E,
							-((float) i * (float) i) / (2.0 * sigma * sigma));
		}
	}

private:
	/**
	 * computes the distance between a fluid point and the wall of a RB
	 */
	float compute_delta(float f_x, float f_y, int direction, float p_x,
			float p_y, float p_r) {

		//printf("direction = (%i, %i)", units[direction][X], units[direction][Y]);
		float l_x = p_x - f_x;
		float l_y = p_y - f_y;
		float s = l_x * units[direction][X] + l_y * units[direction][Y];
		//printf("l = (%f, %f), |s| = %f\n", l_x, l_y, s);
		float l_sq = l_x * l_x + l_y * l_y;
		float m_sq = l_sq - s * s;
		//printf("m^2 = %f, p_r = %f\n", m_sq, p_r);
		float r_minus_m = p_r * p_r - m_sq;
		if (r_minus_m < 0) {
			return 0.5;
		}
		float q = sqrt(r_minus_m);
		float t = s - q;
		//printf("q = %f, t = %f\n", q, t);
		float b_d = sqrt(
				units[direction][X] * units[direction][X]
						+ units[direction][Y] * units[direction][Y]);
		float res = t / b_d;
		//printf("delta = %f\n", res);
		if (res > 1.0) {
			res = 1.0;
			printf("x = (%f, %f), p = (%f, %f)\n", f_x, f_y, p_x, p_y);
			//exit(-1);
		}
		if (res <= 0.0) {
			res = 0.001;
		}
		return res;
	}

private:
	/**
	 * scales down the velocity of the RB that influences the LBM simulation
	 */
	Vector2 normalizeVelocity(float u_x, float u_y, float radius) {
		float abs = (2.0 * PI * radius) * cParameters.rb_velocity_scaling;
		if (abs == 0) {
			return Vector2(u_x, u_y);
		}
		Vector2 res = Vector2(u_x / abs, u_y / abs);
		if (res.magnitude() > rb_v_max) {
			res.normalize();
			res.x *= rb_v_max;
			res.y *= rb_v_max;
		}
		return res;
	}

private:
	void collision_stream() {
		/*
		 *  Collision step (D2Q9 model)
		 */
#pragma omp parallel for
		for (int y_src = 0; y_src < height; y_src++) {
			for (int x_src = 0; x_src < width; x_src++) {

				// normal flag
				const unsigned char flag_src = input_cDataArray2D.getRef(x_src,
						y_src);
				// RB flag
				const unsigned char pflag_src = particle_cDataArray2D.getRef(
						x_src, y_src);
				// incoming distributions
				float *f_src = &simulationData_BufferSrc.getRef(x_src, y_src);

				if (flag_src == 1) { // obstacle
					continue; // does not stream to other cells
				}
				if (pflag_src == 6) { // obstacle
					output_velocity.getRef(x_src, y_src, 0) = 0.0;
					output_velocity.getRef(x_src, y_src, 1) = 0.0;
					continue; // does not stream to other cells
				}

				float* new_f_src;
				float temp_f_src[9];

				float rho = 0.0;
				float u_x = 0.0;
				float u_y = 0.0;

				// Does not need calculation of new probability distributions
				// from in- or outflow
				if (flag_src != 2 && flag_src != 3) {
					// Calculate macroscopic quantities:
					// density
					float rho_inv = 0.0;
					int ex_offx = 0;
					int ex_offy = 0;

					if (particle_cDataArray2D_last.getRef(x_src, y_src) != 0
							&& (cParameters.initialisation_mode == 1
									|| cParameters.initialisation_mode == 2)) {
						int cnt = 0;
						float rho_n;
						float rho_add = 0.0;
						float off_x_add = 0.0;
						float off_y_add = 0.0;
						// densities are collected from the surrounding cells
						// and offset for extrapolation is calculated
						for (int direction_src = 0; direction_src < 9;
								direction_src++) {
							int x_dst = x_src + units[direction_src][X];
							int y_dst = y_src + units[direction_src][Y];
							// takes data only when cell is fluid cell and was last timestep too
							if (particle_cDataArray2D_last.getRef(x_dst, y_dst)
									== 0
									&& input_cDataArray2D.getRef(x_dst, y_dst)
											== 0) {
								rho_add += density_BufferSrc.getRef(x_dst,
										y_dst);
								off_x_add += (float) units[direction_src][X];
								off_y_add += (float) units[direction_src][Y];
								cnt += 1;
							}
						}
						if (cnt > 0) { // interpolation can be only used when there is data
							rho_n = rho_add / (float) cnt;
							ex_offx = roundf(off_x_add / (float) cnt);
							ex_offy = roundf(off_y_add / (float) cnt);
							density_BufferDst.getRef(x_src, y_src) = rho;
							rho_inv = 1 / rho;
						} else { // if there are no densities to interpolate use standard values
							rho_n = 1.0;
							rho_inv = 1.0;
							density_BufferDst.getRef(x_src, y_src) = rho;
						}
						// u is taken from the RB
						CParticle *np = getNearestParticle(x_src, y_src);
						// check for null pointer
						Vector2 norm_u_n;
						if (np != NULL) {
							norm_u_n = normalizeVelocity(np->velocity.x, np->velocity.y, np->radius);
						} else {
							// standard initilialization
							 norm_u_n = normalizeVelocity(0.01, 0.01, 4);
						}
						float u_x_n = norm_u_n.x; //np->velocity.x;
						float u_y_n = norm_u_n.y; //np->velocity.y;

						const float u_x_sq_n = u_x_n * u_x_n;
						const float u_y_sq_n = u_y_n * u_y_n;
						const float u_sq_n = u_x_sq_n + u_y_sq_n;

						f_src[0] = weights[C] * rho_n
								* (1.0f - 3.0f / 2.0f * u_sq_n);
						f_src[1] = weights[E] * rho_n
								* (1.0f + 3.0f * u_x_n + 9.0f / 2.0f * u_x_sq_n
										- 3.0f / 2.0f * u_sq_n);
						f_src[2] = weights[N] * rho_n
								* (1.0f + 3.0f * u_y_n + 9.0f / 2.0f * u_y_sq_n
										- 3.0f / 2.0f * u_sq_n);
						f_src[3] = weights[W] * rho_n
								* (1.0f - 3.0f * u_x_n + 9.0f / 2.0f * u_x_sq_n
										- 3.0f / 2.0f * u_sq_n);
						f_src[4] = weights[S] * rho_n
								* (1.0f - 3.0f * u_y_n + 9.0f / 2.0f * u_y_sq_n
										- 3.0f / 2.0f * u_sq_n);
						f_src[5] = weights[NE] * rho_n
								* (1.0f + 3.0f * (u_x_n + u_y_n)
										+ 9.0f / 2.0f * (u_x_n + u_y_n)
												* (u_x_n + u_y_n)
										- 3.0f / 2.0f * u_sq_n);
						f_src[6] = weights[NW] * rho_n
								* (1.0f + 3.0f * (-u_x_n + u_y_n)
										+ 9.0f / 2.0f * (-u_x_n + u_y_n)
												* (-u_x_n + u_y_n)
										- 3.0f / 2.0f * u_sq_n);
						f_src[7] = weights[SW] * rho_n
								* (1.0f + 3.0f * (-u_x_n - u_y_n)
										+ 9.0f / 2.0f * (-u_x_n - u_y_n)
												* (-u_x_n - u_y_n)
										- 3.0f / 2.0f * u_sq_n);
						f_src[8] = weights[SE] * rho_n
								* (1.0f + 3.0f * (u_x_n - u_y_n)
										+ 9.0f / 2.0f * (u_x_n - u_y_n)
												* (u_x_n - u_y_n)
										- 3.0f / 2.0f * u_sq_n);

						// The non-equilibrium initialization adds an non-equilibrium part
						if (cParameters.initialisation_mode == 2) {
							// first-order extrapolation
							int x_ex = x_src + ex_offx;
							int y_ex = y_src + ex_offy;
							for (int i = 0; i < 9; i++) {
								float f_neq = simulationData_BufferSrc.getRef(
										x_ex, y_ex, i);
								f_src[i] = (f_neq + f_src[i]) / 2.0;
							}
						}
					}

					rho = f_src[NW] + f_src[N] + f_src[NE] + f_src[W] + f_src[C]
							+ f_src[E] + f_src[SW] + f_src[S] + f_src[SE];
					// set upper limit for density
					if (rho > fluid_rho_max) {
						rho = fluid_rho_max;
					}
					density_BufferDst.getRef(x_src, y_src) = rho;
					rho_inv = 1 / rho;
					// velocity components
					u_x = (-f_src[NW] + 0.0 + f_src[NE] - f_src[W] + 0.0
							+ f_src[E] - f_src[SW] + 0.0 + f_src[SE]) * rho_inv;
					u_y = (f_src[NW] + f_src[N] + f_src[NE] + 0.0 + 0.0 + 0.0
							- f_src[SW] - f_src[S] - f_src[SE]) * rho_inv;

					// set limit for velocity
					if (u_x * u_x + u_y * u_y > fluid_v_max) {
						Vector2 u(u_x, u_y);
						u.normalize();
						u.scalarMult(fluid_v_max);
						u_x = u.x;
						u_y = u.y;
					}

					const float u_x_sq = u_x * u_x;
					const float u_y_sq = u_y * u_y;
					const float u_sq = u_x_sq + u_y_sq;

					// Updating velocities for visualization
					// velocity vector
					output_velocity.getRef(x_src, y_src, 0) = u_x;
					output_velocity.getRef(x_src, y_src, 1) = u_y;
					// Pressure
					output_pressure.getRef(x_src, y_src, 0) = rho_inv;

					// equilibrium function
					float feq[9] = { weights[C] * rho
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
									+ 9.0f / 2.0f * (u_x + u_y) * (u_x + u_y)
									- 3.0f / 2.0f * u_sq), weights[NW] * rho
							* (1.0f + 3.0f * (-u_x + u_y)
									+ 9.0f / 2.0f * (-u_x + u_y) * (-u_x + u_y)
									- 3.0f / 2.0f * u_sq), weights[SW] * rho
							* (1.0f + 3.0f * (-u_x - u_y)
									+ 9.0f / 2.0f * (-u_x - u_y) * (-u_x - u_y)
									- 3.0f / 2.0f * u_sq), weights[SE] * rho
							* (1.0f + 3.0f * (u_x - u_y)
									+ 9.0f / 2.0f * (u_x - u_y) * (u_x - u_y)
									- 3.0f / 2.0f * u_sq) };

					for (int i = 0; i < 9; i++) {
						temp_f_src[i] = f_src[i] - omega * (f_src[i] - feq[i]);
					}

					new_f_src = temp_f_src;
				} else {
					new_f_src = f_src;
				}

				// Streamingstep
				for (int direction_src = 0; direction_src < 9;
						direction_src++) {
					// Destination coordinates
					int x_dst = x_src + units[direction_src][X];
					int y_dst = y_src + units[direction_src][Y];

					// Periodic wrapping (faster than modulo)
					x_dst = std::min(std::max(x_dst, 0), width - 1);
					y_dst = std::min(std::max(y_dst, 0), height - 1);

					// Get flag at destination
					const unsigned char flag_dst = input_cDataArray2D.getRef(
							x_dst, y_dst);
					const unsigned char pflag_dst =
							particle_cDataArray2D.getRef(x_dst, y_dst);

					int direction_dst = direction_src;
					if (pflag_dst == 6) { // if x_dst, y_dst is RB
						// distribution gets reflected
						direction_dst = inverse[direction_src];
						int pIndex = getParticleIndex(x_dst, y_dst);
						bool usedTrick17 = false;
						if(pIndex < 0){
							direction_dst = inverse[direction_src];
							x_dst = x_src;
							y_dst = y_src;
							// Final destination
							float *f_dst = &simulationData_BufferDst.getRef(x_dst,
									y_dst);
							f_dst[direction_dst] = new_f_src[direction_src];
							usedTrick17 = true;
							continue;
						}
						if(usedTrick17){
							printf("Das hätte bei Trick 17 nicht passieren sollen");
						}
						CParticle *p = input_particles.data->at(pIndex);
						// delta := distance between fluid point and wall
						float delta = compute_delta(x_src, y_src, direction_src,
								p->position.x, p->position.y, p->radius);
						// second point needed for the velocity bounddary conditions
						int x_f2 = x_src + units[direction_dst][X];
						int y_f2 = y_src + units[direction_dst][Y];
						x_f2 = std::min(std::max(x_f2, 0), width - 1);
						y_f2 = std::min(std::max(y_f2, 0), height - 1);
						unsigned char flag_f2 = input_cDataArray2D.getRef(x_f2,
								y_f2);
						float f2_in;
						float f_out_b = new_f_src[direction_src];
						float f_out_f2 = new_f_src[direction_dst];
						// f2 has to be a fluid point to contain valid data
						if (flag_f2 == 1 || flag_f2 == 2 || flag_f2 == 3
								|| flag_f2 == 6) {
							f2_in = 0.0;
						} else {		// data is taken from the last timestep
							float *f2_dst = &simulationData_BufferSrc.getRef(
									x_f2, y_f2);
							f2_in = f2_dst[direction_src];
							if (f2_in == 0.0) {
								printf("keine Source-Daten\n");
								f2_in = f_out_b;
							}
						}

						Vector2 norm_u = normalizeVelocity(p->velocity.x,
								p->velocity.y, p->radius);

						//reflected distribution
						float f_tal = 0.0;

						// there are 3 different Methods
						if (cParameters.f_rep_mode == 0) {
							f_tal = f_out_b
									- (2 / c_2) * weights[direction_src] * rho
											* (units[direction_src][X]
													* norm_u.x
													+ units[direction_src][Y]
															* norm_u.y);
						} else if (cParameters.f_rep_mode == 1) {
							f_tal =
									(1.0 / (1.0 + delta))
											* ((1.0 - delta) * f2_in // falscher Ort zum Auslesen?
											+ delta * f_out_b + delta * f_out_f2
													- 2.0
															* weights[direction_src]
															* rho * (1.0 / c_2)
															* (units[direction_src][X]
																	* norm_u.x
																	+ units[direction_src][Y]
																			* norm_u.y));
						} else if (cParameters.f_rep_mode == 2) {
							if (delta < 0.5) {
								f_tal =
										2.0 * delta * f_out_b
												+ (1.0 - 2.0 * delta) * f2_in
												- (2.0 * (1.0 / c_2) * rho
														* weights[direction_src]
														* (units[direction_src][X]
																* norm_u.x
																+ units[direction_src][Y]
																		* norm_u.y));
							} else {
								f_tal =
										(1.0 / (2.0 * delta)) * f_out_b
												+ (1.0 - (1.0 / (2.0 * delta)))
														* f_out_f2
												- ((1.0 / (delta)) * (1.0 / c_2)
														* rho
														* weights[direction_src]
														* (units[direction_src][X]
																* norm_u.x
																+ units[direction_src][Y]
																		* norm_u.y));
							}
						}

						// preparations for streaming
						new_f_src[direction_src] = f_tal;
						x_dst = x_src;
						y_dst = y_src;

						// Force on RB
						Vector2 F_p;
						F_p.x = units[direction_src][X] * (f_out_b + f_tal);
						F_p.y = units[direction_src][Y] * (f_out_b + f_tal);
						output_particleForces.at(pIndex).x += F_p.x;
						output_particleForces.at(pIndex).y += F_p.y;

						// torque on RB
						Vector2 X_p;
						X_p.x = (float) x_dst - p->position.x;
						X_p.y = (float) y_dst - p->position.y;
						X_p.normalize();
						X_p.scalarMult(p->radius);

//						printf("i = %i: x_dst = %i; rho = %f; radius = %i\n", pIndex, x_dst, rho, p->radius);

						output_particleForces.at(
								pIndex + input_particles.data->size()).x +=
								X_p.x * F_p.y - X_p.y * F_p.x;
						output_particleForces.at(
								pIndex + input_particles.data->size()).y += 0.0;

					} else {
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
					}

					// Final destination
					float *f_dst = &simulationData_BufferDst.getRef(x_dst,
							y_dst);

					f_dst[direction_dst] = new_f_src[direction_src];
				}
			}
		}
	}

private:
	float get_global_pressure(CDataArray2D<float, 9> &array) {
		float p = 0.0f;

		for (int y = 0; y < array.height; y++) {
			for (int x = 0; x < array.width; x++) {
				for (int i = 0; i < 9; i++)
					p += array.getRef(x, y, i);
			}
		}

		return p;
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
				== typeid(CDataArray2D<unsigned char, 1> ).name()) {
			// unpack data
			CDataArray2D<unsigned char, 1> *input =
					i_cPipelinePacket.getPayload<CDataArray2D<unsigned char, 1> >();

			// copy data to input array
			input_cDataArray2D.resize(input->width, input->height);
			input_cDataArray2D.loadData(input->data);
		} else if (i_cPipelinePacket.type_info_name
				== typeid(CDataParticleArray).name()) {
			// unpack data
			CDataParticleArray *input = i_cPipelinePacket.getPayload<
					CDataParticleArray>();

			input_particles = *input;
			particle_set = true;

		} else {
			std::cerr
					<< "ERROR: Simulation input is only able to process (char,1) flag arrays"
					<< std::endl;
			exit(-1);
		}

	}
	/*
	 public:
	 void writeFluidData(){

	 FILE *ptr_file;
	 char buf[1000];

	 //std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converter;
	 //std::string converted_str = converter.to_bytes( fileName );

	 ptr_file =fopen("flow_config.fconf","w");
	 if (!ptr_file){
	 printf("Die Datei mit dem Namen  konnte nicht geöffnet werden\n");
	 return;
	 }

	 fprintf(ptr_file, "%i;%i;\n", width, height);
	 for(int i = 0; i < width; ++i){
	 for(int j = 0; j < height; ++j){
	 fprintf(ptr_file, "%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;\n",
	 simulationData_BufferSrc.getRef(i, j, 0),
	 simulationData_BufferSrc.getRef(i, j, 1),
	 simulationData_BufferSrc.getRef(i, j, 2),
	 simulationData_BufferSrc.getRef(i, j, 3),
	 simulationData_BufferSrc.getRef(i, j, 4),
	 simulationData_BufferSrc.getRef(i, j, 5),
	 simulationData_BufferSrc.getRef(i, j, 6),
	 simulationData_BufferSrc.getRef(i, j, 7),
	 simulationData_BufferSrc.getRef(i, j, 8),
	 simulationData_BufferDst.getRef(i, j, 0),
	 simulationData_BufferDst.getRef(i, j, 1),
	 simulationData_BufferDst.getRef(i, j, 2),
	 simulationData_BufferDst.getRef(i, j, 3),
	 simulationData_BufferDst.getRef(i, j, 4),
	 simulationData_BufferDst.getRef(i, j, 5),
	 simulationData_BufferDst.getRef(i, j, 6),
	 simulationData_BufferDst.getRef(i, j, 7),
	 simulationData_BufferDst.getRef(i, j, 8),
	 density_BufferSrc.getRef(i, j),
	 density_BufferDst.getRef(i, j));
	 }
	 }

	 fclose(ptr_file);
	 return;
	 }

	 public:
	 void readFluidData(){

	 FILE *ptr_file;
	 char buf[1000];

	 ptr_file =fopen("flow_config.fconf","r");
	 if (!ptr_file){
	 printf("Die Datei mit dem Namen konnte nicht geöffnet werden\n");
	 return;
	 }
	 fgets(buf, 1000, ptr_file);
	 std::string topline(buf);
	 int start = 0;
	 int end = topline.find(";", start);
	 int cwidth = atoi(topline.substr(start, end).c_str());
	 start = end + 1;
	 end = topline.find(";", start);
	 int cheight = atoi(topline.substr(start, end).c_str());
	 if(height != cheight || width != cwidth){
	 printf("inkompatibles Feld wurde eingelesen");
	 return;
	 }
	 for(int i = 0; i < cwidth; ++i){
	 for(int j = 0; j < cheight; ++j){
	 if(fgets(buf,1000, ptr_file)==NULL){
	 printf("Datei Fehlerhaft (zu wenig Daten)");
	 return;
	 }
	 std::string line(buf);
	 start = 0;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 0) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 1) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 2) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 3) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 4) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 5) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 6) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 7) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferSrc.getRef(i, j, 8) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 0) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 1) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 2) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 3) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 4) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 5) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 6) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 7) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 simulationData_BufferDst.getRef(i, j, 8) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 density_BufferSrc.getRef(i, j) = atof(line.substr(start, end).c_str());
	 start = end+1;
	 end = line.find(";", start);
	 density_BufferDst.getRef(i, j) = atof(line.substr(start, end).c_str());
	 }
	 }

	 fclose(ptr_file);
	 return;
	 }

	 private:
	 void writeDataArray(){
	 FILE *ptr_file;
	 char buf[1000];
	 //std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converter;
	 //std::string converted_str = converter.to_bytes( fileName );

	 ptr_file =fopen("data.op","w");
	 if (!ptr_file){
	 printf("Die Datei mit dem Namen  konnte nicht geöffnet werden\n");
	 return;
	 }
	 float* gweights = (float*)malloc((2 * filter_width + 1) * sizeof(float));
	 computeGaussianWeights(gweights, filter_width, filter_sigma);
	 for(int i = dataStart; i < dataEnd + 1; ++i){
	 float pvalue = wdata.at(i - dataStart);
	 if(filter){
	 pvalue = 0.0;
	 for(int k = -filter_width; k <= filter_width; ++k){
	 if(i + k >= dataStart && i + k < dataEnd){
	 pvalue += gweights[k+filter_width] * wdata.at(i + k - dataStart);
	 }else{
	 pvalue += gweights[k+filter_width] * wdata.at(i - k - dataStart);
	 }
	 }
	 }
	 fprintf(ptr_file, "%i %f\n",
	 i, pvalue);

	 }
	 free(gweights);
	 fclose(ptr_file);
	 return;
	 }

	 public:
	 void writeData(char* filename, std::vector<float> data){
	 if(timestep < dataStart || timestep > dataEnd){
	 return;
	 }
	 FILE *ptr_file;

	 if(timestep == dataStart){
	 ptr_file =fopen(filename,"w");
	 }else{
	 ptr_file =fopen(filename,"a");
	 }

	 if (!ptr_file){
	 printf("Die Datei mit dem Namen  konnte nicht geöffnet werden\n");
	 exit(-1);
	 //return;
	 }
	 char* strin;
	 int num = asprintf(&strin, "%i", timestep);
	 for(int i = 0; i < data.size(); ++i){
	 char* save = (char*)malloc(num * sizeof(char));
	 std::strcpy(save, strin);
	 //printf("save = %s\n", save);
	 num = asprintf(&strin, "%s %f", save, data.at(i));
	 //printf("%s, %s, %f\n", strin, save, data.at(i));
	 free(save);
	 }
	 //printf("%s\n", strin);
	 fprintf(ptr_file, "%s\n", strin);

	 fclose(ptr_file);
	 return;
	 }
	 */

};

#endif /* CSTAGE_FLUIDSIMULATION_LBM_HPP_ */
