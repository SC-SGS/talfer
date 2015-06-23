/*
 * CStage_RigidBodySimulation.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */
/* comment for Dominik
 segfault behoben

 in color stehen jetzt Werte die von Game.hpp/Game_version2.hpp gesetzt werden und stehen für verschieden sprites der Visualisierung (sollte noch von der Visualisierung umgesetzt werden)

 ansonsten wüsst ich jetzt nichts was bei uns noch falschläuft
 -eventuell kinect force noch anpassen
 */

#ifndef CSTAGE_RIGIDBODYSIMULATION_HPP_
#define CSTAGE_RIGIDBODYSIMULATION_HPP_

#include "vector"
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
class CStage_RigidBodySimulation: public CPipelineStage {
	/**
	 * global parameters
	 */
	CParameters &cParameters;

	std::vector<CParticle*> particles;
	std::vector<std::array<int, 2>> inflows;
	CDataArray2D<unsigned char, 1> input_cDataArray2D;
	CDataArray2D<float, 2> velocity_cDataArray2D;
	Vector2 kinectForce;
	std::vector<Vector2> forceArray;
	std::vector<float> torqueArray;

	clock_t now;
	clock_t lastClock;
	clock_t lastSpawn;

	Game* game;
public:
	/**
	 * constructor
	 */
	CStage_RigidBodySimulation(CParameters &i_cParameters) :
			CPipelineStage("RigidBody"), cParameters(i_cParameters) {
		lastSpawn = 0;
		lastClock = clock();
		kinectForce.x = 0;
		kinectForce.y = 0;

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
		now = clock();
		float frameDuration = (now - lastClock) / (float) CLOCKS_PER_SEC;
		float leftoverTime = ((int) (frameDuration * 1000)) % TIMESTEP;
		int steps = frameDuration * 1000 / TIMESTEP;
		lastClock = now
				- (leftoverTime / (float) 1000 * (float) CLOCKS_PER_SEC);
		for (int i = 0; i < steps; i++) {
			simulation_timestep();
			check_for_collision();
		}
		if (((now - lastSpawn) / (float) CLOCKS_PER_SEC > 1)
				&& (particles.size() < cParameters.max_number_rb)) {
			create_rigid_body();
		}

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
					v.scalarMult(weights[radius - abs(i)][radius - abs(j)]);
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
//		p1change.scalarMult(0.9);        // eingefügt: randabstoßungsskalierung
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

	void check_for_collision() {
		for (unsigned int i = 0; i < particles.size(); i++) {
			// outflow collision
			Vector2 * collision = checkFlagArea(particles.at(i)->position, 3,
					particles.at(i)->radius);
			if (collision != nullptr) {
				float lifetime = (clock() - particles.at(i)->spawn_time)
						/ CLOCKS_PER_SEC;
				game->particle_destroyed(
						1 / (float) particles.at(i)->inverseMass,
						particles.at(i)->bumps, lifetime, false,
						particles.at(i)->position.x,
						particles.at(i)->position.y);
				particles.erase(particles.begin() + i);
				i--;
			} else {
				collision = checkFlagArea(particles.at(i)->position, 4,
						particles.at(i)->radius);

				if (collision != nullptr) {
					float lifetime = (clock() - particles.at(i)->spawn_time)
							/ CLOCKS_PER_SEC;
					game->particle_destroyed(
							1 / (float) particles.at(i)->inverseMass,
							particles.at(i)->bumps, lifetime, true,
							particles.at(i)->position.x,
							particles.at(i)->position.y);
					particles.erase(particles.begin() + i);
					i--;

				}
			}
			delete collision;
		}

		for (unsigned int i = 0; i < particles.size(); i++) {
			// particle collision
			for (unsigned int l = i + 1; l < particles.size(); l++) {
				if (particles.at(i)->position.distanceTo(
						particles.at(l)->position)
						< (particles.at(i)->radius + particles.at(l)->radius)
								* (particles.at(i)->radius
										+ particles.at(l)->radius)) {
					resolve_sphere_collision(particles.at(i), particles.at(l),
							1.0);
				}
			}
			// boundary collision
			Vector2 * collision = checkFlagArea(particles.at(i)->position, 1,
					particles.at(i)->radius);
			if (collision != nullptr) {
				resolve_obstacle_collision(particles.at(i), 1.0);
			}
			delete collision;
		}
	}

	Vector2* checkFlagArea(Vector2 position, int flag, int radius) {
		//TODO check circular area
//		printf("width  = %i, %f pm %i\n", input_cDataArray2D.width, position.x, radius);
//		printf("height = %i, %f pm %i\n", input_cDataArray2D.height, position.y, radius);
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
		//float density = 1000;
		float frameDuration = TIMESTEP / (float) 1000;
		for (unsigned int i = 0; i < particles.size(); i++) {
			float oldmass = particles.at(i)->inverseMass;
			particles.at(i)->inverseMass = oldmass
					* cParameters.rb_acceleration_scaling;

			Vector2 force = kinectForce;
			force.scalarMult(1000 / cParameters.rb_acceleration_scaling);
			particles.at(i)->addForce(force);

			particles.at(i)->integrate(frameDuration);
			particles.at(i)->inverseMass = oldmass;
			//particles.at(i)->forceAccum.x = forceArray.at(i).x;
			//particles.at(i)->forceAccum.y = forceArray.at(i).y;
			/*
			 particles.at(i)->integrate(frameDuration);
			 particles.at(i)->infiniteField(input_cDataArray2D.width,
			 input_cDataArray2D.height);
			 // apply fluid force on particles
			 if (velocity_cDataArray2D.isValidData()
			 && isFluid(particles.at(i)->position.x,
			 particles.at(i)->position.y)) {
			 Vector2 particleVelocity = particles.at(i)->velocity;
			 Vector2 fluidVelocity;
			 fluidVelocity.x = velocity_cDataArray2D.getRef(
			 particles.at(i)->position.x,
			 particles.at(i)->position.y, 0);
			 fluidVelocity.y = velocity_cDataArray2D.getRef(
			 particles.at(i)->position.x,
			 particles.at(i)->position.y, 1);
			 particleVelocity.scalarMult(1 / (float) 1000);
			 fluidVelocity.subVector(particleVelocity);
			 float length = fluidVelocity.magnitude();
			 fluidVelocity.normalize();
			 fluidVelocity.scalarMult(
			 (density / 2) * length * length
			 * particles.at(i)->area);
			 particles.at(i)->addForce(fluidVelocity);

			 float rotationChange = 0;
			 int currentRadius = particles.at(i)->radius;
			 if (isFluid(particles.at(i)->position.x + currentRadius,
			 particles.at(i)->position.y)) {
			 float rightDownVelocity = velocity_cDataArray2D.getRef(
			 CParticle::sinfiniteField(
			 particles.at(i)->position.x + currentRadius,
			 input_cDataArray2D.width),
			 CParticle::sinfiniteField(
			 particles.at(i)->position.y,
			 input_cDataArray2D.height), 1);
			 rotationChange += rightDownVelocity;
			 }

			 if (isFluid(particles.at(i)->position.x - currentRadius,
			 particles.at(i)->position.y)) {
			 float leftUpVelocity = -1
			 * velocity_cDataArray2D.getRef(
			 CParticle::sinfiniteField(
			 particles.at(i)->position.x
			 - currentRadius,
			 input_cDataArray2D.width),
			 CParticle::sinfiniteField(
			 particles.at(i)->position.y,
			 input_cDataArray2D.height), 1);
			 rotationChange += leftUpVelocity;
			 }
			 if (isFluid(particles.at(i)->position.x,
			 particles.at(i)->position.y + currentRadius)) {
			 float bottomToLeftVelocity = -1
			 * velocity_cDataArray2D.getRef(
			 CParticle::sinfiniteField(
			 particles.at(i)->position.x,
			 input_cDataArray2D.width),
			 CParticle::sinfiniteField(
			 particles.at(i)->position.y
			 + currentRadius,
			 input_cDataArray2D.height), 0);
			 rotationChange += bottomToLeftVelocity;
			 }
			 if (isFluid(particles.at(i)->position.x,
			 particles.at(i)->position.y - currentRadius)) {
			 float topToRightVelocity = velocity_cDataArray2D.getRef(
			 CParticle::sinfiniteField(
			 particles.at(i)->position.x,
			 input_cDataArray2D.width),
			 CParticle::sinfiniteField(
			 particles.at(i)->position.y - currentRadius,
			 input_cDataArray2D.height), 0);
			 rotationChange += topToRightVelocity;
			 }
			 particles.at(i)->addRotationSpeed(
			 ((1000 * rotationChange)
			 - particles.at(i)->rotationspeed) / 800);
			 if (i == 0) {
			 }
			 }
			 // apply kinect force on particles
			 Vector2 force = kinectForce;
			 force.scalarMult(1000);
			 particles.at(i)->addForce(force);
			 */
		}
	}

private:
	void create_rigid_body() {
		if (inflows.size() == 0)
			return;
		std::array<int, 2> point = inflows.at(rand() % inflows.size());
		CParticle *particle = new CParticle(0, 4, now);
		particle->position.x = point.at(0);
		particle->position.y = point.at(1);

		particle->velocity.x = (1 + (rand() % 2)) * (1 - 2 * (rand() % 2))
				* (cParameters.rb_velocity_scaling * 1.5 + 0.05);
		particle->velocity.y = (1 + (rand() % 2)) * (1 - 2 * (rand() % 2))
				* (cParameters.rb_velocity_scaling * 1.5 + 0.05);

		particle->acceleration.x = 0;
		particle->acceleration.y = 0;
		game->setNextParticle(particle);
		particle->inverseMass = 1.0
				/ (float) (particle->radius * particle->radius);
		particles.push_back(particle);

		lastSpawn = clock();

		game->particle_created(1 / (float) particle->inverseMass);
	}

private:
	void create_inflow_vector() {
		inflows.clear();
		for (int y = 0; y < input_cDataArray2D.height; y++) {
			for (int x = 0; x < input_cDataArray2D.width; x++) {
				if (input_cDataArray2D.getRef(x, y) == 2) {
					std::array<int, 2> point = { { x, y } };
					inflows.push_back(point);
				}
			}
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

			// search for inflows
			create_inflow_vector();
		} else if (i_cPipelinePacket.type_info_name
				== typeid(CFluidSimulationData).name()) {

			CFluidSimulationData *velocityField = i_cPipelinePacket.getPayload<
					CFluidSimulationData>();

			velocity_cDataArray2D.resize(velocityField->velocities.width,
					velocityField->velocities.height);
			velocity_cDataArray2D.loadData(velocityField->velocities.data);

			forceArray = velocityField->getParticleForces();

			for (unsigned int i = 0; i < (forceArray.size() / 2.0); ++i) {
				if (i >= particles.size()) {
					break;
				}
				particles.at(i)->forceAccum.x += forceArray.at(i).x;
				particles.at(i)->forceAccum.y += forceArray.at(i).y;
			}
			torqueArray = std::vector<float>();
			for (unsigned int i = (forceArray.size() / 2.0); i < forceArray.size(); ++i) {
				torqueArray.push_back(forceArray.at(i).x);
			}
			for (unsigned int i = 0; i < torqueArray.size(); ++i) {
				if (i >= particles.size()) {
					break;
				}
				particles.at(i)->rotationacc = torqueArray.at(i)
						/ ((particles.at(i)->radius * particles.at(i)->radius)
								/ (particles.at(i)->inverseMass
										* cParameters.rb_acceleration_scaling
										* 2.0));
//				printf("i = %i -> rotationacc = %f, torque = %f\n", i, particles.at(i)->rotationacc,
//						torqueArray.at(i));
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

#endif /* CSTAGE_RIGIDBODYSIMULATION_HPP_ */
