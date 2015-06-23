/*
 * CFluidSimulationData.hpp
 *
 *  Created on: 07.01.2015
 *      Author: axel
 */

#ifndef CFLUIDSIMULATIONDATA_HPP_
#define CFLUIDSIMULATIONDATA_HPP_


#include "CPipelinePacket.hpp"
#include "CDataArray2D.hpp"
#include <vector>

class CFluidSimulationData: public CPipelinePacket {
public:



	CDataArray2D<float, 2> velocities;

	std::vector<Vector2> particleForces;

	/**
	 * changed flag
	 * has to be manually set and is intended for the flag image output
	 */
	bool has_changed;

	/**
	 * return true if data is valid
	 */
	bool isValidData() {
		return velocities.data != nullptr;
	}

	/**
	 * return a pointer of to this instance
	 */
	void *getPayloadRaw() {
		return this;
	}

	inline CFluidSimulationData() :
			has_changed(true) {

		particleForces = std::vector<Vector2>();
		CPipelinePacket::setPacketTypeInfoName(typeid(*this).name());
	}

	inline CFluidSimulationData(std::vector<Vector2> pforces, CDataArray2D<float, 2> vels) :
				has_changed(true) {

		velocities.resize(vels.width, vels.height);
		velocities.loadData(vels.data);
		this->setParticleForces(pforces);
		//printf("%s", typeid(*this).name());
		CPipelinePacket::setPacketTypeInfoName(typeid(*this).name());
	}

	inline CFluidSimulationData(const CFluidSimulationData& other) :
			has_changed(true) {
		CPipelinePacket::setPacketTypeInfoName(typeid(*this).name());
		resize(other.velocities.width, other.velocities.height);
		loadData(other.velocities.data);
	}

	inline CFluidSimulationData& operator=(const CFluidSimulationData& other) {
		if (this != &other) {
			resize(other.velocities.width, other.velocities.height);
			loadData(other.velocities.data);
		}

		return *this;
	}

	void swap(CFluidSimulationData& other) {
		velocities.swap(other.velocities);
		std::swap(other.particleForces, particleForces);
		std::swap(other.has_changed, has_changed);
	}

	/**
	 * fill with data
	 */
	inline void loadData(void *i_data		///< data to store to this 2d array
			) {
		velocities.loadData(i_data);
	}

	/**
	 * resize the array
	 */
	inline void resize(int i_width,	///< new width
			int i_height	///< new height
			) {
		velocities.resize(i_width, i_height);
	}

	/**
	 * return reference to array element
	 */
	inline float &getRef(int i_x, int i_y) const {
		return velocities.getRef(i_x, i_y);
	}

	/**
	 * return reference to component of array element
	 */
	inline float &getRef(int i_x, int i_y, int i_component) const {
		return velocities.getRef(i_x, i_y, i_component);
	}

	/**
	 * return reference to component of array element
	 */
	inline float &getClampedRef(int i_x, int i_y) {
		return velocities.getClampedRef(i_x, i_y);
	}

	/**
	 * return reference to component of array element
	 */
	inline float &getClampedRef(int i_x, int i_y, int i_component) {
		return velocities.getClampedRef(i_x, i_y, i_component);
	}

	inline std::vector<Vector2> getParticleForces(){
		return std::vector<Vector2>(particleForces);
	}

	inline void setParticleForces(std::vector<Vector2> forces){
		particleForces = std::vector<Vector2>(forces);
	}

	inline void cleanup() {
		velocities.cleanup();
		has_changed = true;
	}
	virtual ~CFluidSimulationData() {
		cleanup();
	}
};



#endif /* CFLUIDSIMULATIONDATA_HPP_ */
