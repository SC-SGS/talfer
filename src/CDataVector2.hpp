/*
 * CDataVector2.h
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CDataVector2_H_
#define CDataVector2_H_

#include <typeinfo>
#include <string.h>
#include "Vector2.hpp"
#include "CPipelinePacket.hpp"

/**
 * data storage array for a 2D Vector
 */
class CDataVector2: public CPipelinePacket {
public:
	/**
	 * data storage
	 */
	Vector2* data;

	/**
	 * return true if data is valid
	 */
	bool isValidData() {
		return true;
	}

	/**
	 * return a pointer of to this instance
	 */
	void *getPayloadRaw() {
		return data;
	}

	inline CDataVector2() {
		data = new Vector2();
		CPipelinePacket::setPacketTypeInfoName(typeid(*this).name());
	}

	virtual ~CDataVector2() {
	}
};

#endif /* CDataVector2_H_ */
