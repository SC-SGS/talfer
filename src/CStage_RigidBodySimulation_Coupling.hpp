#ifndef CSTAGE_RIGIDBODYSIMULATION_COUPLING_HPP_
#define CSTAGE_RIGIDBODYSIMULATION_COUPLING_HPP_

#include <vector>
#include "array"

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "SDL/SDL_image.h"
#include "CDataArray2D.hpp"
#include "CDataVector2.hpp"
#include "CDataDrawingInformation.hpp"

#include <time.h>
#if defined(GAME2)
#include "Game_version2.hpp"
#else
#include "Game.hpp"
#endif

#define INITIAL_VELOCITY_X 5
#define INITIAL_VELOCITY_Y 5
#define TIMESTEP 10

/**
 * class providing rigid body simulation
 *
 * calcuated rigid body postions are sent over the pipeline
 */
class CStage_RigidBodySimulation_Coupling: public CPipelineStage {
	/**
	 * global parameters
	 */
	CParameters &cParameters;

	std::vector<CParticle*> particles;
	CDataArray2D<unsigned char, 1> input_cDataArray2D;
	std::vector<Vector2> forceArray;
	std::vector<float> torqueArray;
	Vector2 kinectForce;

	clock_t init;

	Vector2 startCoords;
	bool p_placed;

	Game* game;
public:
	/**
	 * constructor
	 */
	CStage_RigidBodySimulation_Coupling(CParameters &i_cParameters) :
			CPipelineStage("RigidBody"), cParameters(i_cParameters), startCoords(60, 100.0), p_placed(false) {

		init = clock();

		kinectForce.x = 0;
		kinectForce.y = 0;

		/*
		CParticle *particle = new CParticle(0, 50, clock());

		particle->position = startCoords;

		particle->velocity.x = 0.0;
		particle->velocity.y = 0.0;

		particle->acceleration.x = 0;
		particle->acceleration.y = 0;

		particle->inverseMass = 1
						/ (float) (particle->radius * particle->radius);

		particles.push_back(particle);
		*/

		game = Game::getInstance();
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push() {
		CDataParticleArray pipe_particles;
		pipe_particles.data = &particles;
		CPipelineStage::pipeline_push((CPipelinePacket&) pipe_particles);
	}

public:
	void main_loop_callback() {
		// insert and move Particles here

		simulation_timestep();

		pipeline_push();
	}

	void clear() {
		particles.clear();
	}

	void resolve_sphere_collision(CParticle* p1, CParticle* p2,
			float restitution) {
		Vector2 normal = p1->getPosition();
		normal.subVector(p2->position);
		normal.normalize();

		Vector2 relativeVelocity = p1->getVelocity();
		relativeVelocity.subVector(p2->velocity);
		float separatingVelocity = relativeVelocity.scalarProduct(normal);

		if (separatingVelocity > 0)
			return;

		float newSepVelocity = -separatingVelocity * restitution;
		float deltaVelocity = newSepVelocity - separatingVelocity;

		float totalInverseMass = p1->inverseMass + p2->inverseMass;
		if (totalInverseMass <= 0)
			return;

		float impulse = deltaVelocity / totalInverseMass;
		Vector2 impulsePerIMass = normal;
		impulsePerIMass.scalarMult(impulse);

		Vector2 p1change = impulsePerIMass;
		p1change.scalarMult(p1->inverseMass);
		p1change.addVector(p1->velocity);
		p1->setVelocity(p1change);

		Vector2 p2change = impulsePerIMass;
		p2change.scalarMult(-1);
		p2change.scalarMult(p2->inverseMass);
		p2change.addVector(p2->velocity);
		p2->setVelocity(p2change);

		p1->bumps++;
		p2->bumps++;
		game->particle_bumped();
	}

	void resolve_obstacle_collision(CParticle* p1, float restitution) {
		int weights[3][3] = { { 1, 4, 7 }, { 4, 16, 26 }, { 7, 26, 41 } };
		Vector2 p = p1->position;
		Vector2 normal;
		Vector2 v;
		int radius = p1->radius;
		Vector2 collisionPoint;
// TODO implement better collisionPoint detection
		for (int i = -radius; i <= radius; i++) {
			for (int j = -radius; j <= radius; j++) {
				if (isWall(p.x + i, p.y + j)) {
					collisionPoint = Vector2(p.x + i, p.y + j);
				}
			}
		}
		radius = 2;
		for (int i = -radius; i <= radius; i++) {
			for (int j = -radius; j <= radius; j++) {
				if (!isWall(collisionPoint.x + i, collisionPoint.y + j)) {
					v = Vector2(i, j);
					v.scalarMult((float)weights[radius - abs(i)][radius - abs(j)]);
					normal.addVector(v);
				}
			}
		}

		normal.normalize();

		Vector2 relativeVelocity = p1->getVelocity();
		// relativeVelocity.subVector(p2->velocity);
		float separatingVelocity = relativeVelocity.scalarProduct(normal);
		if (separatingVelocity > 0)
			return;

		float newSepVelocity = -separatingVelocity * restitution;
		float deltaVelocity = newSepVelocity - separatingVelocity;

		float totalInverseMass = p1->inverseMass;
		if (totalInverseMass <= 0)
			return;

		float impulse = deltaVelocity / totalInverseMass;
		Vector2 impulsePerIMass = normal;
		impulsePerIMass.scalarMult(impulse);

		Vector2 p1change = impulsePerIMass;
		p1change.scalarMult(p1->inverseMass);
		p1change.addVector(p1->velocity);
		p1->setVelocity(p1change);

		p1->bumps++;
		p1->rotationspeed *= -1;
		game->particle_bumped();

//		Vector2 p2change = impulsePerIMass;
//		p2change.scalarMult(p2->inverseMass);
//		p2change.addVector(p2->velocity);
//		p2->setVelocity(p2change);
	}

	inline bool isWall(int x, int y) {
		return input_cDataArray2D.getRef(
				CParticle::sinfiniteField(x, input_cDataArray2D.width),
				CParticle::sinfiniteField(y, input_cDataArray2D.height)) == 1;
	}

	inline bool isFluid(int x, int y) {
		return input_cDataArray2D.getRef(
				CParticle::sinfiniteField(x, input_cDataArray2D.width),
				CParticle::sinfiniteField(y, input_cDataArray2D.height)) == 0;
	}

	Vector2* checkFlagArea(Vector2 position, int flag, int radius) {
//TODO check circular area
		for (int j = -radius; j <= radius; j++) {
			for (int k = -radius; k <= radius; k++) {
				if (input_cDataArray2D.getRef(
						CParticle::sinfiniteField(position.x + j,
								input_cDataArray2D.width),
						CParticle::sinfiniteField(position.y + k,
								input_cDataArray2D.height)) == flag) {
					return new Vector2(position.x + j, position.y + k);
				}
			}
		}
		return nullptr;
	}

	void simulation_timestep() {
		if(!p_placed && ((clock()-init)/CLOCKS_PER_SEC) > 0 //&& false
				){

			CParticle *particle = new CParticle(0, 10, clock());
			particle->position = startCoords;
			particle->velocity.x = 0.0;
			particle->velocity.y = 0.0;
			particle->acceleration.x = 0;
			particle->acceleration.y = 0;
			particle->inverseMass = 1
							/ (float) (particle->radius * particle->radius);
			particles.push_back(particle);
			p_placed = true;

		}
		float frameDuration = TIMESTEP / (float) 1000;
		for (unsigned int i = 0; i < particles.size(); i++) {
			float oldmass = particles.at(i)->inverseMass;
			particles.at(i)->inverseMass = oldmass
					* cParameters.rb_acceleration_scaling;
			particles.at(i)->integrate(frameDuration);
			particles.at(i)->inverseMass = oldmass;
			particles.at(i)->forceAccum.x = 0.0;
			particles.at(i)->forceAccum.y = 0.0;
		}

	}

public:
	/**
	 * process incoming pipeline input.
	 *
	 * the only input processed so far is from the video output stage to
	 * draw something into the image.
	 */
	void pipeline_process_input(CPipelinePacket &i_cPipelinePacket) {
		// we are currently only able to process "unsigned char,3" data arrays.
		if (i_cPipelinePacket.type_info_name
				== typeid(CDataArray2D<unsigned char, 1> ).name()) {

			// unpack data
			CDataArray2D<unsigned char, 1> *input =
					i_cPipelinePacket.getPayload<CDataArray2D<unsigned char, 1> >();

			// copy data to input array
			input_cDataArray2D.resize(input->width, input->height);
			input_cDataArray2D.loadData(input->data);
		} else if (i_cPipelinePacket.type_info_name
				== typeid(CFluidSimulationData).name()) {

			CFluidSimulationData *fluidData =
					i_cPipelinePacket.getPayload<CFluidSimulationData>();

			forceArray = fluidData->getParticleForces();

			for(unsigned int i = 0; i < (forceArray.size()/2.0); ++i){
				particles.at(i)->forceAccum.x += forceArray.at(i).x;
				particles.at(i)->forceAccum.y += forceArray.at(i).y;
			}
			torqueArray = std::vector<float>();
			for(unsigned int i = (forceArray.size()/2.0); i < forceArray.size(); ++i){
				torqueArray.push_back(forceArray.at(i).x);
			}
			for(unsigned int i = 0; i < torqueArray.size(); ++i){
				particles.at(i)->rotationacc = torqueArray.at(i) / (((float)(particles.at(i)->radius * particles.at(i)->radius)) / (particles.at(i)->inverseMass * 2));
			}

		} else if (i_cPipelinePacket.type_info_name
				== typeid(CDataVector2).name()) {
			kinectForce = *i_cPipelinePacket.getPayload<Vector2>();
			//printf("kinect force 2 is %f and %f \n", kinectForce.x, kinectForce.y);
		} else {
			std::cerr
					<< "ERROR: Simulation input is only able to process (char,1) or (float, 2) flag arrays"
					<< std::endl;
			exit(-1);
		}
	}
};

#endif
