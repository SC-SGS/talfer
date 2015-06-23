#ifndef CSTAGE_DYNAMICVELOCITYFIELD
#define CSTAGE_DYNAMICVELOCITYFIELD

#include <iostream>
#include "math.h"

#include "CPipelineStage.hpp"
#include "CParameters.hpp"
#include "CDataArray2D.hpp"

class CStage_DynamicVelocityField: public CPipelineStage {
private:

	CParameters &cParameters;

	CDataArray2D<float, 2> output_cDataArray2D_f2;

	//Time
	float time;

	//Timestep
	const float timestep;

public:

	CStage_DynamicVelocityField(CParameters &i_cParameters) :
			CPipelineStage("Dynamic Velocity Field"), cParameters(
					i_cParameters), time(0), timestep(0.1) {
		output_cDataArray2D_f2.resize(640, 480);

		for (int i = 0; i < 640; i++)
			for (int j = 0; j < 480; j++) {
				output_cDataArray2D_f2.getRef(i, j, 0) = 0;
				output_cDataArray2D_f2.getRef(i, j, 1) = 0;
			}

	}

	void pipeline_process_input(CPipelinePacket &/*i_cPipelinePacket*/) {
		std::cerr << "ERROR: Static Velocity Field can't process any input."
				<< std::endl;
		exit(-1);
	}

private:
	void setVectorField() {
		for (int i = 0; i < 640; i++)
			for (int j = 0; j < 480; j++) {
				float norm = sin(time - sqrt(i * i + j * j) / 100.0) / 1000;
				output_cDataArray2D_f2.getRef(i, j, 0) = i * norm;
				output_cDataArray2D_f2.getRef(i, j, 1) = j * norm;

				//std::cout << output_cDataArray2D_f2.getRef(i,j,1) << " "<< std::flush;	
			}

	}

public:
	void main_loop_callback() {
		//Set the vector field
		setVectorField();
		time = time + timestep;
		pipeline_push();
	}

	void pipeline_push() {
		CPipelineStage::pipeline_push(
				(CPipelinePacket&) output_cDataArray2D_f2);
	}
};

#endif
