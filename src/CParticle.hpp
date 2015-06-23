#ifndef CPARTICLE_HPP
#define CPARTICLE_HPP
/*
 * CParticle.hpp
 *
 * Created on: Sep 10, 2013
 * Author: rigid team
 */

#include <cmath>
#include "Vector2.hpp"
#include <stdio.h>
#include <time.h>

class CParticle {

public:

	float damping;
	int radius;
	float inverseMass;
	float area;
	int color;
	int bumps;
	clock_t spawn_time;

	float rotation;
	float rotationspeed;
	float rotationacc;
	Vector2 position;
	Vector2 velocity;
	Vector2 acceleration;

	Vector2 forceAccum;

	CParticle() :
			damping(1), radius(4), inverseMass(0.100), area(400), color(0), bumps(
					0), rotation(0), rotationspeed(0), rotationacc(0) {
		spawn_time = clock();
		forceAccum.x = 0;
		forceAccum.y = 0;
	}

	CParticle(int color, int radius, clock_t spawn_time) :
			damping(1), radius(radius), inverseMass(0.100), area(400), color(
					color), bumps(0), spawn_time(spawn_time), rotation(0), rotationspeed(
					0), rotationacc(0) {
		forceAccum.x = 0;
		forceAccum.y = 0;
	}

	void integrate(float duration) {
//		printf("acc=%f, ra = (%f), rv = (%f), v = (%f, %f),  x = (%f, %f)\n", acceleration.magnitude(), rotationacc, rotationspeed, velocity.x, velocity.y, position.x, position.y);
		if (duration <= 0.0) {
			return;
		}

		position.addScaledVector(velocity, duration);

		Vector2 resultingAcc = acceleration;
		resultingAcc.addScaledVector(forceAccum, inverseMass);

		Vector2 resultingVel = velocity;
		resultingVel.addScaledVector(resultingAcc, duration);

		velocity = resultingVel;


		float correctDamping = pow(damping, duration);
		velocity.scalarMult(correctDamping);

		rotationspeed += (rotationacc * duration);
		rotation += (rotationspeed * duration);
		while (rotation < 0) {
			rotation += 360;
		}
		while (rotation > 360) {
			rotation -= 360;
		}
//		printf("acc=%f, ra = (%f), rv = (%f), v = (%f, %f),  x = (%f, %f)\n", acceleration.magnitude(), rotationacc, rotationspeed, velocity.x, velocity.y, position.x, position.y);
//		printf("Force (%f, %f), 1/m = %f, ra = (%f), rv = (%f), v = (%f, %f),  x = (%f, %f)\n", forceAccum.x, forceAccum.y, inverseMass, rotationacc, rotationspeed, velocity.x, velocity.y, position.x, position.y);

		clearAccumulator();
	}

	static int sinfiniteField(int value, int max) {
		if (value < 0) {
			return value + max;
		} else if (value > max) {
			return value - max;
		}
		return value;
	}

	static float sinfiniteField(float value, int max) {
		if (value < 0) {
			return value + max;
		} else if (value > max) {
			return value - max;
		}
		return value;
	}

	void infiniteField(int width, int height) {
		position.x = sinfiniteField(position.x, width);
		position.y = sinfiniteField(position.y, height);
	}

	void clearAccumulator() {
		forceAccum.clear();
	}

	void addForce(Vector2 newForce) {
		forceAccum.addVector(newForce);
	}

	Vector2 getVelocity() {
		return velocity;
	}

	Vector2 getPosition() {
		return position;
	}

	void setVelocity(Vector2 newVelocity) {
		velocity.x = newVelocity.x;
		velocity.y = newVelocity.y;
	}

	void addRotationSpeed(float rotationChange) {
		rotationspeed += rotationChange;
		if (rotationspeed > 360 * 3) {
			rotationspeed = 360 * 3;
		} else if (rotationspeed < -360 * 3) {
			rotationspeed = -360 * 3;
		}
	}

};
#endif
