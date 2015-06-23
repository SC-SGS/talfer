#ifndef VECTOR2_HPP_
#define VECTOR2_HPP_

#include <cmath>

class Vector2{

public:

	float x;
	float y;

	Vector2(){
		x = 0;
		y = 0;
	}

	Vector2(float x_value, float y_value){
		x = x_value;
		y = y_value;
	}

	float magnitude(){
		return sqrt(x*x + y*y);
	}
	
	void scalarMult(float scalar){
		x *= scalar;
		y *= scalar;
	}

	void addVector(Vector2 vec){
		x += vec.x;
		y += vec.y;
	}

	void subVector(Vector2 vec){
		x -= vec.x;
		y -= vec.y;
	}

	void normalize(){
		float length = magnitude();
		
		if(length > 0.0){
			scalarMult(1.0 / length);
		}
	}

	void addScaledVector(Vector2 vec, float scale){
		vec.x *= scale;
		vec.y *= scale;
		addVector(vec);
	}
	
	float scalarProduct(Vector2 vector){
		return (x * vector.x + y * vector.y);
	}
	
	void clear(){

		x = 0;
		y = 0;
	}	
	
	float distanceTo(Vector2 vec){
		return ((x - vec.x) * (x - vec.x) + (y - vec.y) * (y - vec.y));
	}
};
#endif /* VECTOR2_HPP */
