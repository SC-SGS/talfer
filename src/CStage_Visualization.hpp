/*
 * CStage_Visualization.hpp
 *
 * Author: Martin Schreiber
 * Author: Patrick von Steht
 * Author: Dominik Leichtle
 *
 */

#ifndef CSTAGE_VELOCITYVISUALIZATION
#define CSTAGE_VELOCITYVISUALIZATION

#include <iostream>
#include <string>
#include <list>
#include "math.h"
#include <stdlib.h>

#include "CPipelineStage.hpp"
#include "CParameters.hpp"
#include "CDataArray2D.hpp"
#include "CDataParticleArray.hpp"
#include "CFluidSimulationData.hpp"
#include "hsv.h"

class CStage_Visualization: public CPipelineStage {
private:

	CParameters &cParameters;

	//Velocity Field
	CDataArray2D<float, 2> input_cDataArray2D_f2;

	//Flag Field
	CDataArray2D<unsigned char, 1> input_cDataArray2D_uc1;

	//Velocity Norm Field
	CDataArray2D<float, 1> cDataArray2D_velocityNorm;

	//Rigid Body Data
	CDataParticleArray input_cDataParticleArray;

	//Mean Velocity Field over x updates
	CDataArray2D<float, 2> cDataArray2D_meanvelocity;

	//Output Bitmap
	CDataArray2D<unsigned char, 3> output_cDataArray2D_uc3;

	//Kinect input force
	Vector2 kinectForce;
	char kinectNumberOfZeros;

	//Kinect input flag
	bool kinectInput;

	//Width and Height
	int width;
	int height;

	//Arrow length
	static const int arrowLength = 10; //resetted, originally 7
	//Arrow width
	static const int arrowWidth = 0; //min = 0

	//Arrow head radius
	static const int arrowHeadRadius = 1;

	// IconDistances
	static const int iconDistanceX = 20; //resetted, originally 7
	static const int iconDistanceY = iconDistanceX;

	//Jitter
	static const int jitter = iconDistanceX;

	//Maximum of velocity magnitudes
	float maxVelocityValue;
	float maxVelocityValueInverse; //The inverse
	const float maxVelocityWindow = 100.0f;

	//List of points where the arrows are drawn
	CDataArray2D<char, 1> arrowPoints;

	//First input flag
	bool firstInput;

	//Threshold constant for float==0 detection
	const float magnitudeThreshold; //0=disabled (normal mode), now set to 1E-10 to handle empty scenarios
	const float arrowThreshold; //normal value = 0.001

	//Flags: "Is there already input?"
	bool velocityInput, rigidInput, flagFieldInput;

	//H-value of HSVcolor Blue
	static const int hsvBlueH = 180;

	//Max Color Value after ... of the max velocity
	static const float relativeOfVelocityMax;

	//H-value of HSVcolor that represents velocity field of magnitude zero
	static const int hsvZeroMagnitudeH = hsvBlueH;

	//Rigid Body Sprite
#define NSprites 4
	CDataArray2D<unsigned char, 3> rigidBodyIcon[NSprites];
	CDataArray2D<unsigned char, 3> zoomedRigidBodyIcon[NSprites];

	static const std::string rigidBodyIconFile;

	//Color theme
	static const int colorTheme = 0; //0=colorful (normal), 1=pink/white (marmelade)

	//Bool array for the transparent arrow in the middle
	CDataArray2D<bool, 1> transparentPixels;

public:

	CStage_Visualization(CParameters &i_cParameters) :
			CPipelineStage("Visualization"), cParameters(i_cParameters), kinectNumberOfZeros(
					0), kinectInput(false), width(640), height(480), maxVelocityValue(
					0.0), maxVelocityValueInverse(1000000), firstInput(true), magnitudeThreshold(
					0.0000000001), arrowThreshold(0.001), velocityInput(false), rigidInput(
					false), flagFieldInput(false) {
		input_cDataParticleArray.data = new std::vector<CParticle*>();

		for (int i = 0; i < NSprites; i++)
			load_Sprite(rigidBodyIcon[i],
					rigidBodyIconFile + std::to_string(i));
	}

	~CStage_Visualization() {
		delete input_cDataParticleArray.data;
	}

private:
	void zoom_Sprite(CDataArray2D<unsigned char, 3> &sprite, float zoomFactor,
			CDataArray2D<unsigned char, 3> &zoomedSprite) {
		int width = (int) (sprite.width * zoomFactor);
		int height = (int) (sprite.height * zoomFactor);

		zoomedSprite.resize(width, height);

		float zoomFactorInverse = 1.0f / zoomFactor;

		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++) {
				int x = (int) (i * zoomFactorInverse);
				int y = (int) (j * zoomFactorInverse);
				zoomedSprite.getRef(i, j, 0) = sprite.getRef(x, y, 0);
				zoomedSprite.getRef(i, j, 1) = sprite.getRef(x, y, 1);
				zoomedSprite.getRef(i, j, 2) = sprite.getRef(x, y, 2);
			}

	}

private:
	void zoom_Sprite(CDataArray2D<unsigned char, 3> &sprite, int newRadius,
			CDataArray2D<unsigned char, 3> &zoomedSprite) {
		float oldRadius = (sprite.height + sprite.width) / 4.0f;
		float zoomFactor = newRadius / oldRadius;
		zoom_Sprite(sprite, zoomFactor, zoomedSprite);
	}

private:
	void load_Sprite(CDataArray2D<unsigned char, 3> &sprite,
			const std::string &i_path) {
		SDL_Surface *image = IMG_Load(i_path.c_str());
		if (!image) {
			std::cerr << "IMG_Load: " << IMG_GetError();
			exit(-1);
		}

		if (image->format->BitsPerPixel != 24) {
			// only 24 bpp supported for image input, trying to convert...
			// See http://sdl.beuc.net/sdl.wiki/SDL_PixelFormat for documentation of format struct
			SDL_PixelFormat format = { NULL, 24, 3, 0, 0, 0, 0, 0, 8, 16, 0,
					0xFF, 0xFF00, 0xFF0000, 0, 0, 255 };
			SDL_Surface *reimage = SDL_ConvertSurface(image, &format,
			SDL_SWSURFACE);
			SDL_FreeSurface(image);
			image = reimage;
		}

		if (image->w * 3 % 4 != 0) {
			std::cerr << "Image width time three has to be  devidable by four."
					<< std::endl;
			// Otherwise pixels won't match to array indexing
			exit(1);
		}

		sprite.resize(image->w, image->h);

		SDL_LockSurface(image);
		sprite.loadData(image->pixels);
		SDL_UnlockSurface(image);

		SDL_FreeSurface(image);
	}

	void draw_Sprite(CDataArray2D<unsigned char, 3>& sprite, int x, int y,
			float alpha) //alpha in radians
			{
		float sinAlpha = sin(alpha);
		float cosAlpha = cos(alpha);

		int spriteWidth = sprite.width;
		int spriteHeight = sprite.height;
		int max_width_rotated_picture = (int) ceil(
				sqrt(spriteWidth * spriteWidth + spriteHeight * spriteHeight));

		for (int i = 0; i < max_width_rotated_picture; i++)
			for (int j = 0; j < max_width_rotated_picture; j++) {
				int v = j - max_width_rotated_picture / 2;
				int u = i - max_width_rotated_picture / 2;
				int dx = (int) (v * sinAlpha + u * cosAlpha);
				int dy = (int) (v * cosAlpha - u * sinAlpha);

				if ((abs(dy) < (spriteHeight / 2))
						&& (abs(dx) < (spriteWidth / 2))) {
					int xs = dx + spriteWidth / 2;
					int ys = dy + spriteHeight / 2;
					if (!(sprite.getRef(xs, ys, 0) < 50
							&& sprite.getRef(xs, ys, 1) > 215
							&& sprite.getRef(xs, ys, 2) < 50)) //(0,255,0) green is transparent
						color_point(x + v, y + u, sprite.getRef(xs, ys, 0),
								sprite.getRef(xs, ys, 1),
								sprite.getRef(xs, ys, 2));
				}

			}
	}

	void draw_Sprite(CDataArray2D<unsigned char, 3>& sprite, int x, int y) {
		draw_Sprite(sprite, x, y, 0); //draws the sprite without rotation
	}

	/*static float fastInverseSqrt(float x){
	 long i;
	 float x2,y;
	 const float threehalfs = 1.5f;

	 x2 = x * 0.5F;
	 y = x;
	 i = *( long * ) &y;
	 i = 0x5f3759df - ( i >> 1);
	 y = * ( float * ) &i;
	 y = y*( threehalfs - ( x2 * y * y ));
	 y = y*( threehalfs - ( x2 * y * y ));

	 return y;
	 }*/

	void set_arrowPointArray() {
		arrowPoints.resize(width, height);
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++) {
				arrowPoints.getRef(i, j) = 0;
			}
		for (int j = 0; j < height; j += iconDistanceY)
			for (int i = 0; i < width; i += iconDistanceX) {
				int a = i + rand() % jitter;
				int b = j + rand() % jitter;
				if (a >= 0 && a < width && b >= 0 && b < height)
					arrowPoints.getRef(a, b) = 1;
			}
	}

	void draw_line(int x0, int y0, int x1, int y1, unsigned char r,
			unsigned char g, unsigned char b, int xwidth, int ywidth) {

		int dx = abs(x1 - x0);
		int dy = abs(y1 - y0);
		int sx, sy;
		if (x0 < x1)
			sx = 1;
		else
			sx = -1;
		if (y0 < y1)
			sy = 1;
		else
			sy = -1;
		int err = dx - dy;
		while (true) {
			draw_point(x0, y0, r, g, b, xwidth, ywidth);
			if (x0 == x1 && y0 == y1)
				break;
			int e2 = 2 * err;
			if (e2 > -dy) {
				err = err - dy;
				x0 = x0 + sx;
			}
			if (x0 == x1 && y0 == y1) {
				draw_point(x0, y0, r, g, b, xwidth, ywidth);
				break;
			}
			if (e2 < dx) {
				err = err + dx;
				y0 = y0 + sy;
			}
		}

	}

	void draw_line_transparent(int x0, int y0, int x1, int y1, unsigned char r,
			unsigned char g, unsigned char b, int xwidth, int ywidth,
			float transparency) {

		int dx = abs(x1 - x0);
		int dy = abs(y1 - y0);
		int sx, sy;
		if (x0 < x1)
			sx = 1;
		else
			sx = -1;
		if (y0 < y1)
			sy = 1;
		else
			sy = -1;
		int err = dx - dy;
		while (true) {
			draw_point_transparent(x0, y0, r, g, b, xwidth, ywidth,
					transparency);
			if (x0 == x1 && y0 == y1)
				break;
			int e2 = 2 * err;
			if (e2 > -dy) {
				err = err - dy;
				x0 = x0 + sx;
			}
			if (x0 == x1 && y0 == y1) {
				draw_point_transparent(x0, y0, r, g, b, xwidth, ywidth,
						transparency);
				break;
			}
			if (e2 < dx) {
				err = err + dx;
				y0 = y0 + sy;
			}
		}

	}

	void set_maxVelocityValue() {
		float currentMaxVelocityValue = 0.0f;
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (currentMaxVelocityValue	< cDataArray2D_velocityNorm.getRef(i, j, 0)) {
					currentMaxVelocityValue = cDataArray2D_velocityNorm.getRef(i, j, 0);
				}
			}
		}
		// ------------------------------------
		// smoothed scaling of colors for visualization
		if (firstInput) {
			maxVelocityValue = currentMaxVelocityValue;
		} else {
			// approx. gliding mean
			maxVelocityValue += (currentMaxVelocityValue - maxVelocityValue) / maxVelocityWindow;
		}
//		printf("max velocity = %f \n", maxVelocityValue);
//		maxVelocityValue = 0.005;
		// ------------------------------------

		if (maxVelocityValue != 0)
			maxVelocityValueInverse = 1 / maxVelocityValue;
		else
			maxVelocityValueInverse = 1000000;
	}

	void draw_point(int x, int y, unsigned char r, unsigned char g,
			unsigned char b, int xwidth, int ywidth) {
		for (int dx = -floor(xwidth / 2.0); dx <= ceil(xwidth / 2.0); dx++)
			for (int dy = -floor(ywidth / 2.0); dy <= ceil(ywidth / 2.0);
					dy++) {
				color_point(x + dx, y + dy, r, g, b);
			}
	}

	void draw_point_unsafe(int x, int y, unsigned char r, unsigned char g,
			unsigned char b, int xwidth, int ywidth) {
		for (int dx = -floor(xwidth / 2.0); dx <= ceil(xwidth / 2.0); dx++)
			for (int dy = -floor(ywidth / 2.0); dy <= ceil(ywidth / 2.0);
					dy++) {
				color_point_unsafe(x + dx, y + dy, r, g, b);
			}
	}

	void draw_point_transparent(int x, int y, unsigned char r, unsigned char g,
			unsigned char b, int xwidth, int ywidth, float transparency) {
		for (int dx = -floor(xwidth / 2.0); dx <= ceil(xwidth / 2.0); dx++)
			for (int dy = -floor(ywidth / 2.0); dy <= ceil(ywidth / 2.0);
					dy++) {
				if (!transparentPixels.getRef(x + dx, y + dy)) {
					unsigned char dr =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(
									x + dx, y + dy, 0) * (1 - transparency)
									+ transparency * r));
					unsigned char dg =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(
									x + dx, y + dy, 1) * (1 - transparency)
									+ transparency * g));
					unsigned char db =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(
									x + dx, y + dy, 2) * (1 - transparency)
									+ transparency * b));
					color_point(x + dx, y + dy, dr, dg, db);

					transparentPixels.getRef(x + dx, y + dy) = true;
				}
			}
	}

	void draw_filledCircle(int x, int y, unsigned char r, unsigned char g,
			unsigned char b, int radius) {
		for (int dx = -radius; dx <= radius; dx++) {
			int ywidth = ceil(sqrt(radius * radius - dx * dx));
			for (int dy = -ywidth; dy <= ywidth; dy++)
				color_point(x + dx, y + dy, r, g, b);
		}
	}

	void draw_filledCircle_transparent(int x, int y, unsigned char r,
			unsigned char g, unsigned char b, int radius, float transparency) //transparency between 0 and 1
			{
		for (int dx = -radius; dx <= radius; dx++) {
			int ywidth = ceil(sqrt(radius * radius - dx * dx));
			for (int dy = -ywidth; dy <= ywidth; dy++) {
				if (!transparentPixels.getRef(x + dx, y + dy)) {
					unsigned char dr =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(
									x + dx, y + dy, 0) * (1 - transparency)
									+ transparency * r));
					unsigned char dg =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(
									x + dx, y + dy, 1) * (1 - transparency)
									+ transparency * g));
					unsigned char db =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(
									x + dx, y + dy, 2) * (1 - transparency)
									+ transparency * b));
					color_point(x + dx, y + dy, dr, dg, db);
					transparentPixels.getRef(x + dx, y + dy) = true;
				}
			}
		}
	}

	void draw_filledCircle(int x, int y, RgbColor c, int radius) {
		draw_filledCircle(x, y, c.r, c.g, c.b, radius);
	}

	void draw_filledCircle(int x, int y, HsvColor c, int radius) {
		draw_filledCircle(x, y, HsvToRgb(c), radius);
	}

	inline void color_point(int x, int y, unsigned char r, unsigned char g,
			unsigned char b) {
		if (x >= 0 && x < width && y >= 0 && y < height) {
			output_cDataArray2D_uc3.getRef(x, y, 0) = r;
			output_cDataArray2D_uc3.getRef(x, y, 1) = g;
			output_cDataArray2D_uc3.getRef(x, y, 2) = b;
		}
	}

	inline void color_point_unsafe(int x, int y, unsigned char r,
			unsigned char g, unsigned char b) {
		output_cDataArray2D_uc3.getRef(x, y, 0) = r;
		output_cDataArray2D_uc3.getRef(x, y, 1) = g;
		output_cDataArray2D_uc3.getRef(x, y, 2) = b;
	}

	void color_point(int x, int y, RgbColor c) {
		color_point(x, y, c.r, c.g, c.b);
	}

	void color_point(int x, int y, HsvColor c) {
		color_point(x, y, HsvToRgb(c));
	}

	void clearOutput() {
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++)
				for (int k = 0; k < 3; k++)
					output_cDataArray2D_uc3.getRef(i, j, k) = 255;
	}

	void drawIconVisualization() {

		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++) {
				float v1 = 2 * (cDataArray2D_meanvelocity.getRef(i, j, 0));
				float v2 = 2 * (cDataArray2D_meanvelocity.getRef(i, j, 1));
				v1 = v1 * arrowPoints.getRef(i, j);
				v2 = v2 * arrowPoints.getRef(i, j);
				float sqr_norm = v1 * v1 + v2 * v2;
				if (sqr_norm > 0.00001) {
					float norm = sqrt(v1 * v1 + v2 * v2);
					if (norm > arrowThreshold) {
						v1 = (int) ((arrowLength * v1) / norm);
						v2 = (int) ((arrowLength * v2) / norm);
						draw_line(i, j, i + v1, j + v2, 255, 255, 255,
								arrowWidth, arrowWidth);
						draw_filledCircle(i + v1, j + v2, 255, 255, 255,
								arrowHeadRadius);
					}
				}
			}

	}

	void drawParticleTrace(int i, int j) {
		//Initial value is completely random, just dont choose "0" if you want a trace
		int v1 = 1234;
		int v2 = 1234;
		while (i >= 0 && i < width && j >= 0 && j < height
				&& (v1 != 0 || v2 != 0)) {
			v1 = (int) ceil(input_cDataArray2D_f2.getRef(i, j, 0));
			v2 = (int) ceil(input_cDataArray2D_f2.getRef(i, j, 1));
			draw_line(i, j, i + v1, j + v2, 0, 0, 200, 1, 1);
			i += v1;
			j += v2;
		}
	}

	void drawParticleTraceVisualization() {
		// We want to draw 50 lines
		for (int k = 0; k < 50; k++) {
			int i, j;
			// Choose a random starting point for the stream line
			i = rand() % width;
			j = rand() % height;

			drawParticleTrace(i, j);
		}
	}

	void drawRigidBodies() {
		const float rotation_enforcement_factor = 10;
		for (std::vector<CParticle*>::iterator it =
				input_cDataParticleArray.data->begin();
				it != input_cDataParticleArray.data->end(); ++it) {
			unsigned int c = (*it)->color;
			if (!(c == 0 || c == 1 || c == 2))
				c = 2;
			if (c < NSprites) {
				zoom_Sprite(rigidBodyIcon[c], (*it)->radius,
						zoomedRigidBodyIcon[c]);
				//if ((*it)->rotation != (*it)->rotation){
				//(*it)->rotation = 0;
				//std::cout << "nan entdeckt" << std::endl;
				//}else
				draw_Sprite(zoomedRigidBodyIcon[c],
						(int) floor((*it)->position.x),
						(int) floor((*it)->position.y),
						((*it)->rotation / 360.0f * -2.0f * M_PI
								* rotation_enforcement_factor));
			}
//			{
//				HsvColor c;
//				c.h = 210; c.s = 255; c.v = 255;
//            	draw_filledCircle((int) floor((*it)->position.x), (int) floor((*it)->position.y), c, (*it)->radius);  				
//			}

			//std::cout << (*it)->inverseMass << std::endl;

		}
	}

	void drawRigidBodyTraces() {
		for (std::vector<CParticle*>::iterator it =
				input_cDataParticleArray.data->begin();
				it != input_cDataParticleArray.data->end(); ++it) {
			drawParticleTrace((*it)->position.x, (*it)->position.y);
		}
	}

	void drawVelocityNormField() {
		// Array is indexed by a part of the float bits. This part contains the exponent bits to get a log range.
		// TODO: Suppress type-pun pointer warning.
		static RgbColor colors[1 << 13];

		manipulate_velocityNormField();

		if (maxVelocityValue >= magnitudeThreshold) {
			for (int i = 0; i < (1 << 13); i++) {
				HsvColor c;
				c.v = 200;
				c.s = 250;
				int iv = i << 17;
				float *vp = (float*) &iv;
				float v = *vp;
//			    float v = *((float*) &iv);
				c.h = computeHsvHvalue(v, 7);
				colors[i] = HsvToRgb(c);
			}
		} else {
			for (int i = 0; i < (1 << 13); i++) {
				HsvColor c;
				c.v = 200;
				c.s = 250;
				;
				c.h = hsvZeroMagnitudeH;
				colors[i] = HsvToRgb(c);
			}

		}
		float *velocitynorm = &cDataArray2D_velocityNorm.getRef(0, 0, 0);
		RgbColor *out = (RgbColor*) &output_cDataArray2D_uc3.getRef(0, 0);
		for (int i = 0; i < (width * height); i++) {
			float v = *velocitynorm++;
			int *idp = (int*) &v;
			int idx = (*idp & 0x3FFE0000) >> 17;
//		    int idx = (*((int*) &v) & 0x3FFE0000) >> 17;
			*out++ = colors[idx];
		}

	}

	int computeHsvHvalue(float velocity, char /*key*/) {
// 		static double min, max;
// 		if (velocity < min && velocity > 0) {
// 			min = velocity;
// 			printf("min: %f\n", min);	
// 		}
// 		if (velocity > max) {
// 			max = velocity;
// 			printf("max: %f\n", max);	
// 		}

		float a;
		/*		switch(key)
		 {
		 case 0: //easy linear
		 return (int)(hsvZeroMagnitudeH*(1-velocity*maxVelocityValueInverse));
		 case 1: //exponential - not working
		 return (int)(hsvZeroMagnitudeH*exp(-velocity));
		 case 2: //better exponential - well working - slow
		 a = log2(hsvZeroMagnitudeH+1)*maxVelocityValueInverse;
		 return (int)((hsvZeroMagnitudeH+1)*exp(-velocity*a)-1);
		 case 3:
		 a = hsvZeroMagnitudeH/log2(1+maxVelocityValue);
		 return (int)(a*log2(1+maxVelocityValue-velocity));
		 case 4:
		 if(velocity<maxVelocityValue*relativeOfVelocityMax)
		 return (int)(hsvZeroMagnitudeH*(1-velocity*maxVelocityValueInverse/relativeOfVelocityMax));
		 else
		 return 0;
		 case 5: //Quadratic
		 a = hsvZeroMagnitudeH*maxVelocityValueInverse*maxVelocityValueInverse;
		 return a*(velocity-maxVelocityValue)*(velocity-maxVelocityValue);
		 case 6: //to pow 3
		 a = hsvZeroMagnitudeH*maxVelocityValueInverse*maxVelocityValueInverse*maxVelocityValueInverse;
		 return -a*(velocity-maxVelocityValue)*(velocity-maxVelocityValue)*(velocity-maxVelocityValue);*/
//			case 7: //to pow 5
		a = hsvZeroMagnitudeH * maxVelocityValueInverse
				* maxVelocityValueInverse * maxVelocityValueInverse
				* maxVelocityValueInverse * maxVelocityValueInverse;
		return -a * (velocity - maxVelocityValue)
				* (velocity - maxVelocityValue) * (velocity - maxVelocityValue)
				* (velocity - maxVelocityValue) * (velocity - maxVelocityValue);
//			default:
//				return hsvZeroMagnitudeH;
//		}
	}

	/*inline float expApprox(float x)
	 {
	 return 1+x+x*x*0.5f+x*x*x*0.166666f;
	 }*/

	void manipulate_velocityNormField() {
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				float& a = cDataArray2D_velocityNorm.getRef(i, j, 0);
				a = std::min(a, 255.0f);
			}
		}
	}

	void compute_velocityNormField() {
		cDataArray2D_velocityNorm.resize(width, height);

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				float *a = &input_cDataArray2D_f2.getRef(i, j, 0);
				cDataArray2D_velocityNorm.getRef(i, j, 0) = (a[0] * a[0]
						+ a[1] * a[1]);
			}
		}
	}

	void drawFlagField() {
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++) {
				if (input_cDataArray2D_uc1.getRef(i, j) == 1) {
					//half-transparent
					/*unsigned char r = (unsigned char)((output_cDataArray2D_uc3.getRef(i,j,0)+3*0)/4.0f);
					 unsigned char g = (unsigned char)((output_cDataArray2D_uc3.getRef(i,j,1)+3*0)/4.0f);
					 unsigned char b = (unsigned char)((output_cDataArray2D_uc3.getRef(i,j,2)+3*0)/4.0f);
					 color_point(i,j,r,g,b);*/
					color_point(i, j, 0, 0, 0);
				} else if (input_cDataArray2D_uc1.getRef(i, j) == 2) {
					//half-transparent
					unsigned char r =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 0) + 3 * 255) / 4.0f);
					unsigned char g =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 1) + 3 * 0) / 4.0f);
					unsigned char b =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 2) + 3 * 0) / 4.0f);
					color_point(i, j, r, g, b);
				} else if (input_cDataArray2D_uc1.getRef(i, j) == 3) {
					//half-transparent
					unsigned char r =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 0) + 3 * 0) / 4.0f);
					unsigned char g =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 1) + 3 * 0) / 4.0f);
					unsigned char b =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 2) + 3 * 255) / 4.0f);
					color_point(i, j, r, g, b);
				} else if (input_cDataArray2D_uc1.getRef(i, j) == 4) {
					//half-transparent
					unsigned char r =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 0) + 3 * 0) / 4.0f);
					unsigned char g =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 1) + 3 * 255) / 4.0f);
					unsigned char b =
							(unsigned char) ((output_cDataArray2D_uc3.getRef(i,
									j, 2) + 3 * 0) / 4.0f);
					color_point(i, j, r, g, b);
				} else if (input_cDataArray2D_uc1.getRef(i, j) == 4)
					color_point(i, j, 0, 255, 0);
			}
	}

	void updateMeanVelocity() {
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++) {
				cDataArray2D_meanvelocity.getRef(i, j, 0) = 0.9f
						* cDataArray2D_meanvelocity.getRef(i, j, 0)
						+ 0.1f * input_cDataArray2D_f2.getRef(i, j, 0);
				cDataArray2D_meanvelocity.getRef(i, j, 1) = 0.9f
						* cDataArray2D_meanvelocity.getRef(i, j, 1)
						+ 0.1f * input_cDataArray2D_f2.getRef(i, j, 1);
			}
	}

	void drawKinectForce() {
		static const float transparency = 0.7;
		float norm = sqrt(
				kinectForce.x * kinectForce.x + kinectForce.y * kinectForce.y);
		if (norm == 0)
			norm = 1;
		int fx = (int) (kinectForce.x / norm * 100);
		int fy = (int) (kinectForce.y / norm * 100);
		int i = width / 2;
		int j = height / 2;

		transparentPixels.resize(width, height);
		memset(transparentPixels.data, 0, width * height * sizeof(bool));

		draw_line_transparent(i, j, i + fx, j + fy, 255, 255, 255, 10, 10,
				transparency);
		draw_filledCircle_transparent(i + fx, j + fy, 255, 255, 255, 15,
				transparency);
	}

	void drawKinectForceBorder() {
		float norm = sqrt(
				kinectForce.x * kinectForce.x + kinectForce.y * kinectForce.y);
		float x = kinectForce.x / norm;
		float y = kinectForce.y / norm;

		int xpos = (int) width * x;
		int ypos = (int) height * y;

		float whratio = width / (float) height;

		const int size = 30;

		//Kante 1
		if (y > 0 && abs(x) < whratio * y) {
			ypos = height - size;
			xpos = width / 2 + (ypos - height / 2) / y * x;
		}
		//Kante 3
		else if (y < 0 && abs(x) < whratio * abs(y)) {
			ypos = size;
			xpos = width / 2 - (height / 2 - ypos) / y * x;
		}
		//Kante 2
		else if (x > 0 && x > abs(y)) {
			xpos = width - size;
			ypos = height / 2 + (xpos - width / 2) / x * y;
		}
		//Kante 4
		else if (x < 0 && abs(x) > abs(y)) {
			xpos = size;
			ypos = height / 2 - (width / 2 - xpos) / x * y;
		}

		draw_point(xpos, ypos, 255, 255, 255, size, size);
		draw_point(xpos, ypos, 255, 0, 0, size * norm, size * norm);
	}

	void drawColorChoice() {
		for (int x = 1; x <= width; x++) {
			for (int y = height; y <= height + 30; y++) {
				color_point_unsafe(x, y, 50, 50, 50);
			}
		} // Grey Background
		/*draw_point_unsafe(360,(int)(height+15),255,0,0,20,20); // Red
		 draw_point_unsafe(390,(int)(height+15),0,0,255,20,20); // Blue
		 draw_point_unsafe(420,(int)(height+15),0,255,0,20,20); // Green
		 draw_point_unsafe(450,(int)(height+15),0,0,0,20,20); // Black
		 draw_point_unsafe(480,(int)(height+15),255,255,255,20,20); // White (for the Eraser)*/

		draw_point_unsafe(20, (int) (height + 15), 255, 0, 0, 20, 20); // Red
		draw_point_unsafe(50, (int) (height + 15), 0, 0, 255, 20, 20); // Blue
		draw_point_unsafe(80, (int) (height + 15), 0, 255, 0, 20, 20); // Green
		draw_point_unsafe(110, (int) (height + 15), 0, 0, 0, 20, 20); // Black
		draw_point_unsafe(140, (int) (height + 15), 255, 255, 255, 20, 20); // White (for the Eraser)

		// draw the pause button
		draw_point_unsafe(500, (int) (height + 15), 255, 255, 0, 60, 20);

		// draw the image save button
		draw_point_unsafe(580, (int) (height + 15), 255, 0, 255, 60, 20);
	}

public:

	void pipeline_process_input(CPipelinePacket &i_cPipelinePacket) {
		//Velocity Field
		if (i_cPipelinePacket.type_info_name
				== typeid(CFluidSimulationData).name()) {

			CFluidSimulationData* input_fsd = i_cPipelinePacket.getPayload<
					CFluidSimulationData>();

			CDataArray2D<float, 2> *input = &(input_fsd->velocities);

			width = input->width;
			height = input->height;

			input_cDataArray2D_f2.resize(width, height);
			input_cDataArray2D_f2.loadData(input->data);

			compute_velocityNormField();
			set_maxVelocityValue();
//			 std::cout << "Max Velocity Value: " << maxVelocityValue << std::endl << std::flush;

			if (firstInput == true) {
				set_arrowPointArray();
				cDataArray2D_meanvelocity.resize(width, height);
				cDataArray2D_meanvelocity.loadData(input->data);
			} else {
				updateMeanVelocity();
			}

			//Set first input flag
			firstInput = false;

			//Set "theres already an input" flag
			velocityInput = true;
		} else if (i_cPipelinePacket.type_info_name
				== typeid(CDataArray2D<unsigned char, 1> ).name()) {
			CDataArray2D<unsigned char, 1>* input =
					i_cPipelinePacket.getPayload<CDataArray2D<unsigned char, 1>>();

			width = input->width;
			height = input->height;

			input_cDataArray2D_uc1.resize(width, height);
			input_cDataArray2D_uc1.loadData(input->data);

			flagFieldInput = true;

		}
		//Rigid Bodies
		else if (i_cPipelinePacket.type_info_name
				== typeid(CDataParticleArray).name()) {
			CDataParticleArray* input = i_cPipelinePacket.getPayload<
					CDataParticleArray>();

			input_cDataParticleArray.data->assign(input->data->begin(),
					input->data->end());

			rigidInput = true;
		}
		//Kinect Input
		else if (i_cPipelinePacket.type_info_name
				== typeid(CDataVector2).name()) {

#ifdef KINECT_INPUT
			Vector2 newKinectForce = *i_cPipelinePacket.getPayload<Vector2>();

			//Kinect force inertia
			if(newKinectForce.x == 0 && newKinectForce.y == 0)
			{
				kinectNumberOfZeros++;
				if(kinectNumberOfZeros >= 3)
				kinectForce = newKinectForce;
			} else {
				kinectForce = newKinectForce;
				kinectInput = true;
				kinectNumberOfZeros = 0;
			}
#else
			kinectForce = *i_cPipelinePacket.getPayload<Vector2>();
			kinectInput = true;
#endif //KINECT_INPUT

			//printf("Visu: kinect force is %f and %f \n", kinectForce.x, kinectForce.y);

		}
		//Something we don't know
		else {
			std::cerr
					<< "ERROR: Velocity Visualization is only able to process (float,2) arrays"
					<< std::endl;
			//exit(-1);
		}

		output_cDataArray2D_uc3.resize(width, height + 30);
	}

	void main_loop_callback() {
		update_and_draw();
	}

	void update_and_draw() {

		clearOutput();

		/*
		 * Visualize velocity magnitude field
		 */
		if ((firstInput == false) && (velocityInput == true))
			drawVelocityNormField();

		/*
		 * Visualize velocity vector field through small lines (arrows would be nice later on)
		 */
		if ((firstInput == false) && (velocityInput == true))
			drawIconVisualization();

		/*
		 * Visualize flag field
		 */
		if (flagFieldInput == true)
			drawFlagField();

		/*
		 * Visualize velocity vector field through particle traces originating from some random points.
		 */
		//drawParticleTraceVisualization();
		/*
		 * Visualize rigid bodies by big red points.
		 */
		if (rigidInput == true)
			drawRigidBodies();

		/*
		 * Draw particle traces originating at the rigid bodies
		 */
		//drawRigidBodyTraces();
		/*
		 * Draws the kinect force
		 */
		if (kinectInput == true) {
			drawKinectForce();
			drawKinectForceBorder();
		}

		drawColorChoice();

		CPipelineStage::pipeline_push(
				(CPipelinePacket&) output_cDataArray2D_uc3);
	}

	void pipeline_push() {
		CPipelineStage::pipeline_push(
				(CPipelinePacket&) output_cDataArray2D_uc3);
	}

	void save_flag_field_as_image(const char* path) {
		CDataArray2D<unsigned char, 3> cDataArray2D_uc3;
		cDataArray2D_uc3.resize(input_cDataArray2D_uc1.width,
				input_cDataArray2D_uc1.height);

		for (int y = 0; y < cDataArray2D_uc3.height; y++) {
			for (int x = 0; x < cDataArray2D_uc3.width; x++) {
				unsigned char *d = &cDataArray2D_uc3.getRef(x, y, 0);
				switch (input_cDataArray2D_uc1.getRef(x, y)) {
				case 0:	// fluid
					d[0] = 255;
					d[1] = 255;
					d[2] = 255;
					break;

				case 1:	// obstacle
					d[0] = 0;
					d[1] = 0;
					d[2] = 0;
					break;

				case 2:	// inflow
					d[0] = 255;
					d[1] = 0;
					d[2] = 0;
					break;

				case 3:	// outflow
					d[0] = 0;
					d[1] = 0;
					d[2] = 255;
					break;

				case 4:	// outflow
					d[0] = 0;
					d[1] = 255;
					d[2] = 0;
					break;
				}
			}
		}

		char* pixels = reinterpret_cast<char*>(cDataArray2D_uc3.data);
		SDL_Surface *surface = SDL_CreateRGBSurfaceFrom(pixels,
				cDataArray2D_uc3.width, cDataArray2D_uc3.height, 3 * 8, // channels*8, the bits per pixel
				cDataArray2D_uc3.width * 3, // pitch
				0x0000FF, // red mask
				0x00FF00, // green mask
				0xFF0000, // blue mask
				0); // alpha mask (none)
		SDL_SaveBMP(surface, path);
	}

};

const float CStage_Visualization::relativeOfVelocityMax = 1 / 2.0;

//const std::string CStage_Visualization::rigidBodyIconFile = "data/rigidBodyIcon.png";
//const std::string CStage_Visualization::rigidBodyIconFile = "data/rigidBodyIcon_tuxPenguin.jpg";
const std::string CStage_Visualization::rigidBodyIconFile =
		"data/rigGraphics/rigIcon.png"; //"rigidBodyIcon_tuxPenguin_supersmall_green.png";
//const std::string CStage_Visualization::rigidBodyIconFile = "data/rigGraphics/kaiserschmarrn_gruen.png";

#endif
