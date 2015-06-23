#ifndef CSTAGE_STATICVELOCITYFIELD
#define CSTAGE_STATICVELOCITYFIELD

#include <iostream>
#include "math.h"

#include "CPipelineStage.hpp"
#include "CParameters.hpp"
#include "CDataArray2D.hpp"

class CStage_StaticVelocityField: public CPipelineStage {
private:

	CParameters &cParameters;

	CDataArray2D<float, 2> output_cDataArray2D_f2;

public:

	CStage_StaticVelocityField(CParameters &i_cParameters) :
			CPipelineStage("Static Velocity Field"), cParameters(i_cParameters) {
		output_cDataArray2D_f2.resize(640, 480);

		for (int i = 0; i < 640; i++)
			for (int j = 0; j < 480; j++) {
				output_cDataArray2D_f2.getRef(i, j, 0) = i / 80.0;
				output_cDataArray2D_f2.getRef(i, j, 1) = j / 60.0 - i / 80.0;
			}
	}

	void pipeline_process_input(CPipelinePacket &/*i_cPipelinePacket*/) {
		std::cerr << "ERROR: Static Velocity Field can't process any input."
				<< std::endl;
		exit(-1);
	}

	void pipeline_push() {
		CPipelineStage::pipeline_push(
				(CPipelinePacket&) output_cDataArray2D_f2);
	}

	void main_loop_callback() {
	}
};

#endif
