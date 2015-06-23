/*
 * CDataParticleArray.h
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CDataParticleArray_H_
#define CDataParticleArray_H_

#include <typeinfo>
#include <string.h>
#include <vector>
#include "CParticle.hpp"
#include "CPipelinePacket.hpp"

/**
 * data storage array for type T and D components per array element
 */
class CDataParticleArray: public CPipelinePacket {
public:
	/**
	 * size of array
	 */
	int size;

	/**
	 * data storage
	 */
	std::vector<CParticle*> *data;

	/**
	 * return true if data is valid
	 */
	bool isValidData() {
		return data != nullptr;
	}

	/**
	 * return a pointer of to this instance
	 */
	void *getPayloadRaw() {
		return this;
	}

	inline CDataParticleArray() :
			size(-1), data(nullptr) {
		CPipelinePacket::setPacketTypeInfoName(typeid(*this).name());
	}

	virtual ~CDataParticleArray() {
	}
};

#endif /* CDataParticleArray_H_ */
