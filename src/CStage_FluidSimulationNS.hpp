/*
 * CStage_ImageInput.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 *
 *      There's already an existing file named CStage_FluidSimulationLBM.hpp
 *		with a skeleton for the LBM simulation.  Please consider using this as
 *		your NS Skeleton.  You get a flag array from the image processing as
 *		an input which can frequently change.  The output should be a velocity
 *		field based on the simulation involving the input flag array.
 */

#ifndef CSTAGE_FLUIDSIMULATION_NS_HPP_
#define CSTAGE_FLUIDSIMULATION_NS_HPP_

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "SDL/SDL_image.h"
#include "CDataArray2D.hpp"
#include "CDataDrawingInformation.hpp"
#include <math.h>
#include <iostream>
#include <cassert>

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_ConfigDefs.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MsrMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "AztecOO.h"

enum Direction {
	CENTER = 16,
	EAST = 8,
	WEST = 4,
	SOUTH = 2,
	NORTH = 1,
	NORTHEAST = 9,
	SOUTHEAST = 10,
	SOUTHWEST = 6,
	NORTHWEST = 5,
	INFLOW = 32,
	OUTFLOW = 64
};

enum Axis {
	X, Y
};

/**
 * class providing static image input
 *
 * the image is send to the pipeline during each main loop iteration
 */
class CStage_FluidSimulationNS: public CPipelineStage {
private:
	const float density = 1000;
	const float viscosity = 1;
	const float gravity[2] = { 0.0, 0.0 };
	const float dx = 1e-1;
	const float dy = 1e-1;
	const float max_dt = 5e-4;
	float dt;
	const float inflowPressure = 5e3;
	const float fluidPressure = 1e2;
	const float outflowPressure = 1e1;
	bool initialized;
	CDataArray2D<char, 1> type;
	CDataArray2D<float, 1> pressure;
	CDataArray2D<float, 2> velocity;
	CDataArray2D<float, 2> velocity_tmp;
public:
	CDataArray2D<unsigned char, 1> input_cDataArray2D;
	CParameters &cParameters;
	CDataArray2D<float, 2>& output_cDataArray2D_f2;

	/**
	 * constructor
	 */
	CStage_FluidSimulationNS(CParameters &i_cParameters) :
			CPipelineStage("FluidSimulationNS"), dt(max_dt), initialized(false), cParameters(
					i_cParameters), output_cDataArray2D_f2(velocity) {
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push() {
		if (cParameters.stage_fluidsimulation_visualize_flagfield)
			CPipelineStage::pipeline_push(
					(CPipelinePacket&) input_cDataArray2D);

		CPipelineStage::pipeline_push(
				(CPipelinePacket&) output_cDataArray2D_f2);
	}

	// own calculators / getters / setter
	float& getPressure(int i, int j, Direction direction) {
		switch (direction) {
		case CENTER:
			return pressure.getRef(normalize(i, X), normalize(j, Y));
			break;
		case EAST:
			return pressure.getRef(normalize(i + 1, X), normalize(j, Y));
			break;
		case WEST:
			return pressure.getRef(normalize(i - 1, X), normalize(j, Y));
			break;
		case SOUTH:
			return pressure.getRef(normalize(i, X), normalize(j - 1, Y));
			break;
		case NORTH:
			return pressure.getRef(normalize(i, X), normalize(j + 1, Y));
			break;
		case NORTHEAST:
			return pressure.getRef(normalize(i + 1, X), normalize(j + 1, Y));
			break;
		case SOUTHEAST:
			return pressure.getRef(normalize(i + 1, X), normalize(j - 1, Y));
			break;
		case SOUTHWEST:
			return pressure.getRef(normalize(i - 1, X), normalize(j - 1, Y));
			break;
		case NORTHWEST:
			return pressure.getRef(normalize(i - 1, X), normalize(j + 1, Y));
			break;
		default:
			assert(false);
			break;
		}
		return pressure.getRef(i, j);
	}
	float& getVelocity(int i, int j, Direction direction, Axis axis) {
		return getVelocity(i, j, direction, axis, velocity);
	}

	float& getVelocity(int i, int j, Direction direction, Axis axis,
			CDataArray2D<float, 2>& source) {
		switch (direction) {
		case CENTER:
			return source.getRef(normalize(i, X), normalize(j, Y), axis);
			break;
		case EAST:
			return source.getRef(normalize(i + 1, X), normalize(j, Y), axis);
			break;
		case WEST:
			return source.getRef(normalize(i - 1, X), normalize(j, Y), axis);
			break;
		case SOUTH:
			return source.getRef(normalize(i, X), normalize(j - 1, Y), axis);
			break;
		case NORTH:
			return source.getRef(normalize(i, X), normalize(j + 1, Y), axis);
			break;
		case NORTHEAST:
			return source.getRef(normalize(i + 1, X), normalize(j + 1, Y), axis);
			break;
		case SOUTHEAST:
			return source.getRef(normalize(i + 1, X), normalize(j - 1, Y), axis);
			break;
		case SOUTHWEST:
			return source.getRef(normalize(i - 1, X), normalize(j - 1, Y), axis);
			break;
		case NORTHWEST:
			return source.getRef(normalize(i - 1, X), normalize(j + 1, Y), axis);
			break;
		default:
			assert(false);
			break;
		}
		return source.getRef(i, j, axis);
	}

	int getLinear(int i, int j, Direction direction) {
		switch (direction) {
		case CENTER:
			break;
		case EAST:
			i += 1;
			break;
		case WEST:
			i -= 1;
			break;
		case SOUTH:
			j -= 1;
			break;
		case NORTH:
			j += 1;
			break;
		case NORTHEAST:
			i += 1;
			j += 1;
			break;
		case SOUTHEAST:
			j -= 1;
			i += 1;
			break;
		case SOUTHWEST:
			j -= 1;
			i -= 1;
			break;
		case NORTHWEST:
			j += 1;
			i -= 1;
			break;
		default:
			assert(false);
			break;
		}
		return normalize(j, Y) * input_cDataArray2D.width + normalize(i, X);
	}

	float diffusive_x(int i, int j) {
		// diffusive term in x-dir
		float d2ud2x = (getVelocity(i, j, EAST, X)
				- 2 * getVelocity(i, j, CENTER, X) + getVelocity(i, j, WEST, X))
				/ (dx * dx);
		float d2ud2y =
				(getVelocity(i, j, NORTH, X) - 2 * getVelocity(i, j, CENTER, X)
						+ getVelocity(i, j, SOUTH, X)) / (dy * dy);
		return viscosity * (d2ud2x + d2ud2y);
	}

	float diffusive_y(int i, int j) {
		// diffusive term in y-dir
		float d2vd2x = (getVelocity(i, j, EAST, Y)
				- 2 * getVelocity(i, j, CENTER, Y) + getVelocity(i, j, WEST, Y))
				/ (dx * dx);
		float d2vd2y =
				(getVelocity(i, j, NORTH, Y) - 2 * getVelocity(i, j, CENTER, Y)
						+ getVelocity(i, j, SOUTH, Y)) / (dy * dy);
		return viscosity * (d2vd2x + d2vd2y);
	}

	float pressure_x(int i, int j) {
		//pressure term in x
		float dpdx = (getPressure(i, j, EAST) - getPressure(i, j, CENTER)) / dx;
		return -(1 / density) * dpdx;
	}

	float pressure_y(int i, int j) {
		//pressure term in y
		float dpdy = (getPressure(i, j, NORTH) - getPressure(i, j, CENTER))
				/ dy;
		return -(1 / density) * dpdy;
	}

	float gravity_x() {
		return gravity[X];
	}
	float gravity_y() {
		return gravity[Y];
	}

	float pow(float x) {
		return x * x;
	}

	float convective_x(int i, int j) {
		double du2dx = (pow(
				(getVelocity(i, j, EAST, X) + getVelocity(i, j, CENTER, X)))
				- pow(
						(getVelocity(i, j, CENTER, X)
								+ getVelocity(i, j, WEST, X)))) / (4 * dx);
		double duvdy =
				((getVelocity(i, j, CENTER, Y) + getVelocity(i, j, EAST, Y))
						* (getVelocity(i, j, CENTER, X)
								+ getVelocity(i, j, NORTH, X))
						- (getVelocity(i, j, SOUTH, Y)
								+ getVelocity(i, j, SOUTHEAST, Y))
								* (getVelocity(i, j, SOUTH, X)
										+ getVelocity(i, j, CENTER, X)))
						/ (4 * dy);
		return du2dx + duvdy;
	}

	float convective_y(int i, int j) {
		double dv2dy = (pow(
				getVelocity(i, j, CENTER, Y) + getVelocity(i, j, NORTH, Y))
				- pow(
						getVelocity(i, j, SOUTH, Y)
								+ getVelocity(i, j, CENTER, Y))) / (4 * dy);
		double duvdx = ((getVelocity(i, j, CENTER, X)
				+ getVelocity(i, j, NORTH, X))
				* (getVelocity(i, j, CENTER, Y) + getVelocity(i, j, EAST, Y))
				- (getVelocity(i, j, WEST, X) + getVelocity(i, j, NORTHWEST, X))
						* (getVelocity(i, j, WEST, Y)
								+ getVelocity(i, j, CENTER, Y))) / (4 * dx);
		return dv2dy + duvdx;
	}

	float get_maxvelocity() {
		//int max_i_x, max_j_x, max_i_y, max_j_y;
		float max = abs(getVelocity(0, 0, CENTER, X));
		for (int j = 0; j < velocity.height; j++) {
			for (int i = 0; i < velocity.width; i++) {
				if (abs(getVelocity(i, j, CENTER, X)) > max) {
					max = abs(getVelocity(i, j, CENTER, X));
					//max_i_x = i;
					//max_j_x = j;
				}
				if (abs(getVelocity(i, j, CENTER, Y)) > max) {
					max = abs(getVelocity(i, j, CENTER, Y));
					//max_i_y = i;
					//max_j_y = j;
				}
			}
		}
		return max;
	}

	int normalize(int v, Axis axis, bool quadratic) {
		int norm;
		if (quadratic)
			norm = input_cDataArray2D.width * input_cDataArray2D.height;
		else if (axis == X)
			norm = input_cDataArray2D.width;
		else
			norm = input_cDataArray2D.height;

		if (v < 0)
			v += norm;
		else if (v >= norm)
			v -= norm;

		return v;
	}

	int normalize(int v, Axis axis) {
		return normalize(v, axis, false);
	}

	void simulation_timestep() {
		// 0:fluid, 1:obstacle, 2: inflow, 3:outflow
		if (!input_cDataArray2D.isValidData())
			return;

		// initialize
		if (!initialized) {
			initialized = true;
			pressure.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);
			velocity_tmp.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);
			velocity.resize(input_cDataArray2D.width,
					input_cDataArray2D.height);
			type.resize(input_cDataArray2D.width, input_cDataArray2D.height);
			for (int j = 0; j < pressure.height; j++) {
				for (int i = 0; i < pressure.width; i++) {
					getPressure(i, j, CENTER) = fluidPressure;
					getVelocity(i, j, CENTER, X) = 0;
					getVelocity(i, j, CENTER, Y) = 0;
					getVelocity(i, j, CENTER, X, velocity_tmp) = 0;
					getVelocity(i, j, CENTER, Y, velocity_tmp) = 0;
					type.getRef(i, j) = (input_cDataArray2D.getRef(
							normalize(i, X), normalize(j, Y)) != 1) * CENTER
							+ (input_cDataArray2D.getRef(normalize(i + 1, X),
									normalize(j, Y)) != 1) * EAST
							+ (input_cDataArray2D.getRef(normalize(i - 1, X),
									normalize(j, Y)) != 1) * WEST
							+ (input_cDataArray2D.getRef(normalize(i, X),
									normalize(j - 1, Y)) != 1) * SOUTH
							+ (input_cDataArray2D.getRef(normalize(i, X),
									normalize(j + 1, Y)) != 1) * NORTH
							+ (input_cDataArray2D.getRef(normalize(i, X),
									normalize(j, Y)) == 2) * INFLOW
							+ (input_cDataArray2D.getRef(normalize(i, X),
									normalize(j, Y)) == 3) * OUTFLOW;
				}
			}
		} else {
			//calculate maximum possible timestep size
			dt = (dx < dy ? dx : dy) / (2 * get_maxvelocity()); //TODO
			dt = dt > max_dt ? max_dt : dt;
			std::cout << dt << endl;
		}

		//TODO
		//type.getRef(5,5) = INFLOW+CENTER+EAST+WEST+NORTH+SOUTH;
		//type.getRef(9,9) = OUTFLOW+CENTER+EAST+WEST+NORTH+SOUTH;

		// temporary velocities;
		for (int j = 0; j < velocity_tmp.height; j++) {
			for (int i = 0; i < velocity_tmp.width; i++) {
				getVelocity(i, j, CENTER, X, velocity_tmp) = getVelocity(i, j,
						CENTER, X)
						+ dt
								* (-convective_x(i, j) + diffusive_x(i, j)
										+ gravity_x());
				getVelocity(i, j, CENTER, Y, velocity_tmp) = getVelocity(i, j,
						CENTER, Y)
						+ dt
								* (-convective_y(i, j) + diffusive_y(i, j)
										+ gravity_y());
			}
		}

		// initialize structures
#ifdef EPETRA_MPI
		// Initialize MPI
		MPI_Init(&argc,&argv);
		Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
		Epetra_SerialComm Comm;
#endif
		int NumMyElements = input_cDataArray2D.width
				* input_cDataArray2D.height;
		Epetra_Map Map(NumMyElements, NumMyElements, 0, Comm);
		Epetra_CrsMatrix A(Copy, Map, 5);
		Epetra_Vector x(Map);
		Epetra_Vector b(Map);

		// populate matrix
		double value_A_Y = 1 / (dy * dy);
		double value_A_X = 1 / (dx * dx);
		double value_A_CENTER = -2 * value_A_Y - 2 * value_A_X;
		double value_one = 1;
		double value_minusone = -1;
		double value_minushalf = -0.5;
		for (int j = 0; j < input_cDataArray2D.height; j++) {
			for (int i = 0; i < input_cDataArray2D.width; i++) {
				int globalRow_CENTER = A.GRID(getLinear(i, j, CENTER));
				int globalRow_NORTH = A.GRID(getLinear(i, j, NORTH));
				int globalRow_SOUTH = A.GRID(getLinear(i, j, SOUTH));
				int globalRow_WEST = A.GRID(getLinear(i, j, WEST));
				int globalRow_EAST = A.GRID(getLinear(i, j, EAST));

				if ((type.getRef(i, j) & INFLOW) == INFLOW) {
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
							&globalRow_CENTER);
					b[getLinear(i, j, CENTER)] = inflowPressure;
				} else if ((type.getRef(i, j) & OUTFLOW) == OUTFLOW) {
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
							&globalRow_CENTER);
					b[getLinear(i, j, CENTER)] = outflowPressure;
				} else if ((type.getRef(i, j) & CENTER) == CENTER) { //Fluidzelle
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_A_CENTER,
							&globalRow_CENTER);
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_A_Y,
							&globalRow_NORTH);
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_A_Y,
							&globalRow_SOUTH);
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_A_X,
							&globalRow_WEST);
					A.InsertGlobalValues(globalRow_CENTER, 1, &value_A_X,
							&globalRow_EAST);
					b[getLinear(i, j, CENTER)] = density / dt
							* ((getVelocity(i, j, EAST, X, velocity_tmp)
									- getVelocity(i, j, WEST, X, velocity_tmp))
									/ (2 * dx)
									+ (getVelocity(i, j, NORTH, Y, velocity_tmp)
											- getVelocity(i, j, SOUTH, Y,
													velocity_tmp)) / (2 * dy));

				} else { //Hinderniszelle
					switch (type.getRef(i, j)) {
					case NORTH: //North
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minusone, &globalRow_NORTH);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case SOUTH: //South
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minusone, &globalRow_SOUTH);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case EAST: //East
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minusone, &globalRow_EAST);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case WEST: //West
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minusone, &globalRow_WEST);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case NORTHWEST: //Northwest
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_NORTH);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_WEST);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case NORTHEAST: //Northeast
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_NORTH);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_EAST);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case SOUTHWEST: //Southwest
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_SOUTH);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_WEST);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case SOUTHEAST: //Southeast
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_SOUTH);
						A.InsertGlobalValues(globalRow_CENTER, 1,
								&value_minushalf, &globalRow_EAST);
						b[getLinear(i, j, CENTER)] = 0;
						break;
					case 0: //komplett im Hindernis
						A.InsertGlobalValues(globalRow_CENTER, 1, &value_one,
								&globalRow_CENTER);
						b[getLinear(i, j, CENTER)] = fluidPressure;
						break;
					}
				}
			}
		}

		A.FillComplete();
		//A.Print(std::cout);
		//b.Print(std::cout);
		// run solver

		Epetra_LinearProblem Problem(&A, &x, &b);
		AztecOO Solver(Problem);
		Solver.SetAztecOption(AZ_diagnostics, AZ_none);
		Solver.SetAztecOption(AZ_output, AZ_none);
		Solver.SetAztecOption(AZ_solver, AZ_gmres);
		Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
		//Solver.SetAztecOption(AZ_solver, AZ_gmres);
		//Solver.SetAztecOption(AZ_precond, AZ_Jacobi);
		Solver.Iterate(10000, 1E-3);

		// retrieving solution

		for (int j = 0; j < pressure.height; j++) {
			for (int i = 0; i < pressure.width; i++) {
				getPressure(i, j, CENTER) = x[getLinear(i, j, CENTER)];
			}
		}

		// final velocities

		for (int j = 0; j < velocity.height; j++) {
			for (int i = 0; i < velocity.width; i++) {
				getVelocity(i, j, CENTER, X, velocity) = (getVelocity(i, j,
						CENTER, X, velocity_tmp)
						- dt / (dx * density)
								* (getPressure(i, j, EAST)
										- getPressure(i, j, CENTER)));
				getVelocity(i, j, CENTER, Y, velocity) = (getVelocity(i, j,
						CENTER, Y, velocity_tmp)
						- dt / (dy * density)
								* (getPressure(i, j, NORTH)
										- getPressure(i, j, CENTER)));
			}
		}

		// boundary conditions
		for (int j = 0; j < velocity.height; j++) {
			for (int i = 0; i < velocity.width; i++) {
				// not allowed
				if (type.getRef(i, j) <= CENTER
						&& ((type.getRef(i, j) & (NORTH + SOUTH))
								== (NORTH + SOUTH)
								|| (type.getRef(i, j) & (EAST + WEST))
										== (EAST + WEST))) {
					assert(false);
				}
				switch (type.getRef(i, j)) {
				// not interesting: 1xxxx, 00000
				case EAST:		//EAST
					getVelocity(i, j, CENTER, X) = 0;
					getVelocity(i, j, CENTER, Y) = -getVelocity(i, j, EAST, Y);
					getVelocity(i, j, WEST, X) = 0;
					getVelocity(i, j, SOUTH, Y) = -getVelocity(i, j, SOUTHEAST,
							Y);
					break;
				case WEST:		//WEST
					getVelocity(i, j, CENTER, X) = 0;
					getVelocity(i, j, CENTER, Y) = -getVelocity(i, j, WEST, Y);
					getVelocity(i, j, WEST, X) = 0;
					getVelocity(i, j, SOUTH, Y) = -getVelocity(i, j, SOUTHWEST,
							Y);
					break;
				case SOUTH:		//SOUTH
					getVelocity(i, j, CENTER, X) = -getVelocity(i, j, SOUTH, X);
					getVelocity(i, j, CENTER, Y) = 0;
					getVelocity(i, j, WEST, X) = -getVelocity(i, j, SOUTHWEST,
							X);
					getVelocity(i, j, SOUTH, Y) = 0;
					break;
				case NORTH:		//NORTH
					getVelocity(i, j, CENTER, X) = -getVelocity(i, j, NORTH, X);
					getVelocity(i, j, CENTER, Y) = 0;
					getVelocity(i, j, WEST, X) = -getVelocity(i, j, NORTHWEST,
							X);
					getVelocity(i, j, SOUTH, Y) = 0;
					break;
				case NORTHWEST:		//NORTHWEST
					getVelocity(i, j, CENTER, X) = -getVelocity(i, j, NORTH, X);
					getVelocity(i, j, CENTER, Y) = 0;
					getVelocity(i, j, WEST, X) = 0;
					getVelocity(i, j, SOUTH, Y) = -getVelocity(i, j, SOUTHWEST,
							Y);
					break;
				case SOUTHEAST:		//SOUTHEAST
					getVelocity(i, j, CENTER, X) = 0;
					getVelocity(i, j, CENTER, Y) = -getVelocity(i, j, EAST, Y);
					getVelocity(i, j, WEST, X) = -getVelocity(i, j, SOUTHWEST,
							X);
					getVelocity(i, j, SOUTH, Y) = 0;
					break;
				case NORTHEAST:		//NORTHEAST
					getVelocity(i, j, CENTER, X) = 0;
					getVelocity(i, j, CENTER, Y) = 0;
					getVelocity(i, j, WEST, X) = -getVelocity(i, j, NORTHWEST,
							X);
					getVelocity(i, j, SOUTH, Y) = -getVelocity(i, j, SOUTHEAST,
							Y);
					break;
				case SOUTHWEST:		//SOUTHWEST
					getVelocity(i, j, CENTER, X) = -getVelocity(i, j, SOUTH, X);
					getVelocity(i, j, CENTER, Y) = -getVelocity(i, j, WEST, Y);
					getVelocity(i, j, WEST, X) = 0;
					getVelocity(i, j, SOUTH, Y) = 0;
					break;
				}
			}
		}

		if (false) {
			std::cout << "***** final Types *****" << endl; // TODO
			for (int j = 0; j < input_cDataArray2D.height; j++) {
				for (int i = 0; i < input_cDataArray2D.width; i++) {
					printf("%d\t", type.getRef(i, j));
				}
				std::cout << endl;
			}
		}
		if (false) {
			std::cout << "***** final Pressures *****" << endl; // TODO
			for (int j = 0; j < input_cDataArray2D.height; j++) {
				for (int i = 0; i < input_cDataArray2D.width; i++) {
					printf("%e\t", getPressure(i, j, CENTER));
				}
				std::cout << endl;
			}
		}
		if (false) {
			std::cout << "***** final Velocities *****" << endl; // TODO
			for (int j = 0; j < input_cDataArray2D.height; j++) {
				for (int i = 0; i < input_cDataArray2D.width; i++) {
					if (getVelocity(i, j, CENTER, X) != 0
							|| getVelocity(i, j, CENTER, Y) != 0)
						printf("%e,%e\t", getVelocity(i, j, CENTER, X),
								getVelocity(i, j, CENTER, Y));
				}
				std::cout << endl;
			}
		}
		// push pipeline
		pipeline_push();

#ifdef EPETRA_MPI
		MPI_Finalize();
#endif
	}

	void main_loop_callback() {
		simulation_timestep();
	}

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
};

#endif /* CSTAGE_FLUIDSIMULATION_NS_HPP_ */
