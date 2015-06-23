/*
 * CStage_KeyboardForces.hpp
 *
 *  Created on: Oct 2, 2013
 *      Author: Patrick Buchfink
 */

#ifndef CSTAGE_KEYBOARDFORCES_HPP_
#define CSTAGE_KEYBOARDFORCES_HPP_

//(Wurzel2)/2
#define SQRT2DIV2 0.707106781

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "CDataVector2.hpp"
#include <stdlib.h>

class CStage_KeyboardForces: public CPipelineStage {
	/**
	 * global parameters
	 */
	CParameters &cParameters;
	CDataVector2* cVector;
	float applyiedForceWithKeyboard;

	//Key pressed: True
	int pressedUp;
	int pressedDown;
	int pressedLeft;
	int pressedRight;

public:
	/**
	 * constructor
	 */
	CStage_KeyboardForces(CParameters &i_cParameters) :
			CPipelineStage("KeyboardForces"), cParameters(i_cParameters) {
		applyiedForceWithKeyboard = 0.8;
		cVector = new CDataVector2();
		cVector->data->x = 0.0;
		cVector->data->y = 0.0;
		pressedUp = 0;
		pressedDown = 0;
		pressedRight = 0;
		pressedLeft = 0;
	}

	void main_loop_callback() {
	}

	void keyUp(int i_key) {
		switch (i_key) {
		case SDLK_UP:
			pressedUp = 0;
			break;
		case SDLK_DOWN:
			pressedDown = 0;
			break;
		case SDLK_RIGHT:
			pressedRight = 0;
			break;
		case SDLK_LEFT:
			pressedLeft = 0;
			break;
		default:
			return;
		}
		calculateVector();
		pipeline_push(cVector);
	}

	void keyDown(int i_key) {
		switch (i_key) {
		case SDLK_UP:
			pressedUp = 1;
			break;
		case SDLK_DOWN:
			pressedDown = 1;
			break;
		case SDLK_RIGHT:
			pressedRight = 1;
			break;
		case SDLK_LEFT:
			pressedLeft = 1;
			break;
		default:
			return;
		}
		calculateVector();
		pipeline_push(cVector);
	}

	void pipeline_process_input(CPipelinePacket &/*i_cPipelinePacket*/
	) {
	}

private:
	void calculateVector() {
		cVector->data->x = (pressedRight - pressedLeft)
				* applyiedForceWithKeyboard;
		cVector->data->y = (pressedDown - pressedUp)
				* applyiedForceWithKeyboard;
		if ((pressedRight - pressedLeft) != 0) {
			cVector->data->y *= SQRT2DIV2;
		}
		if ((pressedDown - pressedUp) != 0) {
			cVector->data->x *= SQRT2DIV2;
		}
		//printf("SQRT2DIV2: %f, X: %f, Y: %f", SQRT2DIV2, (1-(pressedDown-pressedUp)*SQRT2DIV2), (1-(pressedRight-pressedLeft)*SQRT2DIV2));
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push(CDataVector2* output_cVector) {
		CPipelineStage::pipeline_push((CPipelinePacket&) *output_cVector);
	}

};

#endif /* CSTAGE_KEYBOARDFORCES_HPP_ */
