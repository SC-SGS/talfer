#ifndef MATRI2X2_HPP_
#define MATRI2X2_HPP_


class Matrix2x2{

public:

	float x_1_1;
	float x_1_2;
	float x_2_1;
	float x_2_2;

	Matrix2x2(){
		x_1_1 = 0.0;
		x_1_2 = 0.0;
		x_2_1 = 0.0;
		x_2_2 = 0.0;
	}

	Matrix2x2(float x1, float x2, float x3, float x4){
		x_1_1 = x1;
		x_1_2 = x2;
		x_2_1 = x3;
		x_2_2 = x4;
	}

	static Matrix2x2 tensorProd(Vector2 vec1, Vector2 vec2){
		return Matrix2x2(vec1.x * vec2.x, vec1.x * vec2.y, vec1.y * vec2.x, vec1.y * vec2.x);
	}

	Matrix2x2 operator+(Matrix2x2 b){
		return Matrix2x2(this->x_1_1 + b.x_1_1, this->x_1_2 + b.x_1_2, this->x_2_1 + b.x_2_1, this->x_2_2 + b.x_2_2);
	}

	Matrix2x2 operator-(Matrix2x2 b){
		return Matrix2x2(this->x_1_1 - b.x_1_1, this->x_1_2 - b.x_1_2, this->x_2_1 - b.x_2_1, this->x_2_2 - b.x_2_2);
	}

	Matrix2x2 operator*(float a){
		return Matrix2x2(a * this->x_1_1, a * this->x_1_2, a * this->x_2_1, a * this->x_2_2);
	}

	float norm(){
		return pow(x_1_1 * x_1_1 + x_1_2 * x_1_2 + x_2_1 * x_2_1 + x_2_2 * x_2_2, 0.5);
	}

};

static Matrix2x2 operator*(float a, Matrix2x2 b){
	return Matrix2x2(a * b.x_1_1, a * b.x_1_2, a * b.x_2_1, a * b.x_2_2);
}

#endif
