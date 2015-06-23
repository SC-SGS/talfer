#ifndef CSTAGE_STATICRIGIDBODIES
#define CSTAGE_STATICRIGIDBODIES

#include <iostream>
#include "math.h"

#include "CPipelineStage.hpp"
#include "CParameters.hpp"
#include "CDataArray2D.hpp"

#include "CDataParticleArray.hpp"

class CStage_StaticRigidBodies: public CPipelineStage {
private:

	CParameters &cParameters;

	CDataParticleArray output_cDataParticleArray;

public:

	CStage_StaticRigidBodies(CParameters &i_cParameters) :
			CPipelineStage("Static Rigid Bodies"), cParameters(i_cParameters) {
		output_cDataParticleArray.data = new std::vector<CParticle*>();

		for (int k = 0; k < 10; k++) {
			CParticle * particle = new CParticle();

			Vector2 pos;
			pos.x = rand() % 640;
			pos.y = rand() % 480;
			particle->position = pos;

			output_cDataParticleArray.data->push_back(particle);
		}

		output_cDataParticleArray.size = output_cDataParticleArray.data->size();
	}

	void pipeline_process_input(CPipelinePacket &/*i_cPipelinePacket*/) {
		std::cerr << "ERROR: Static Rigid Bodies can't process any input."
				<< std::endl;
		exit(-1);
	}

	void pipeline_push() {
		CPipelineStage::pipeline_push(
				(CPipelinePacket&) output_cDataParticleArray);
	}

	void main_loop_callback() {
	}
};

#endif
