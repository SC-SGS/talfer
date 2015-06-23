/*output_cDataArray2D_uc
 * CStage_ImageInput.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CSTAGE_VIDEOOUTPUT_HPP_
#define CSTAGE_VIDEOOUTPUT_HPP_

#include <iostream>
#include <list>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>

#include "CPipelineStage.hpp"
#include "CParameters.hpp"
#include "CDataArray2D.hpp"
#include "CSDLInterface.hpp"
#include "CGlTexture.hpp"
#include "CFont.hpp"
#include "CDataDrawingInformation.hpp"
#include "CStage_Visualization.hpp"

#include "hsv.h"

#ifndef KINECT_INCLUDED
#include "CStage_KeyboardForces.hpp"
#endif

typedef struct _Text {
	std::string msg;
	int ttl;
	int x;
	int y;
	int size;
} Text;

class CStage_VideoOutput: public CPipelineStage, public CSDLInterface {
	CParameters &cParameters;
	//time_t last_button_pressed;
	struct timeval last_button_pressed;
	CStage_Visualization *cStage_Visualization;

	CDataArray2D<unsigned char, 3> output_cDataArray2D_uc3;
	CDataArray2D<float, 2> output_cDataArray2D_f2;
	CDataArray2D<float, 3> output_cDataArray2D_f3;

	CDataArray2D<unsigned char, 4> output_cDataArray2D_uc4_flags;

	CGlTexture cGlTextureFlags;

#ifndef KINECT_INCLUDED
	CStage_KeyboardForces* cStage_KeyboardForces;
#endif

	//Font Drawing Class
	CFont *F;

	//Texts to draw onto screen
	std::list<Text*> *Texts;
	std::string *Points;

	std::string *Rest_time;

	/**
	 * OpenGL handler to texture
	 */
	CGlTexture cGlTexture;

	/**
	 * update the texture
	 */
	bool updateTexture;

	// The active color: the color that is drawn in the window with the mouse
	RgbColor activeColor;

public:
	CStage_VideoOutput(CParameters &i_cParameters) :
			CPipelineStage("ImageInput"), cParameters(i_cParameters),
			//last_button_pressed(0),
			Points(0), Rest_time(0), cGlTexture(GL_TEXTURE_2D, GL_RGB, GL_RGB,
			GL_UNSIGNED_BYTE), updateTexture(false) {
		F = new CFont();
		if (!F->init) {
			std::cerr << "Failed to load font" << std::endl;
			exit(-1);
		}
		F->setFontSize(20);
		Texts = new std::list<Text*>();

		gettimeofday(&last_button_pressed, NULL);
		last_button_pressed.tv_sec = 0;
		last_button_pressed.tv_usec = 0;
	}

	virtual ~CStage_VideoOutput() {
	}

	void connectVisualization(CStage_Visualization* visu) {
		cStage_Visualization = visu;
	}

	/**
	 * incoming data to process from pipeline
	 */
	void pipeline_process_input(CPipelinePacket &i_cPipelinePacket) {

		output_cDataArray2D_uc3.cleanup();
		output_cDataArray2D_f2.cleanup();
		output_cDataArray2D_f3.cleanup();

		/**
		 * RGB image
		 */
		if (i_cPipelinePacket.type_info_name
				== typeid(CDataArray2D<unsigned char, 3> ).name()) {
			CDataArray2D<unsigned char, 3> *input =
					i_cPipelinePacket.getPayload<CDataArray2D<unsigned char, 3> >();

			output_cDataArray2D_uc3.resize(input->width, input->height);
			output_cDataArray2D_uc3.loadData(input->data);

			updateTexture = true;
			return;
		}

		/**
		 * velocity field
		 */
		if (i_cPipelinePacket.type_info_name
				== typeid(CDataArray2D<float, 2> ).name()) {
			CDataArray2D<float, 2> *input = i_cPipelinePacket.getPayload<
					CDataArray2D<float, 2> >();

			output_cDataArray2D_f2.resize(input->width, input->height);
			output_cDataArray2D_f2.loadData(input->data);

			updateTexture = true;
			return;
		}

		/**
		 * velocity field
		 */
		if (i_cPipelinePacket.type_info_name
				== typeid(CDataArray2D<float, 3> ).name()) {
			CDataArray2D<float, 3> *input = i_cPipelinePacket.getPayload<
					CDataArray2D<float, 3> >();

			output_cDataArray2D_f3.resize(input->width, input->height);
			output_cDataArray2D_f3.loadData(input->data);

			updateTexture = true;
			return;
		}

		/**
		 * flag field
		 */
		if (i_cPipelinePacket.type_info_name
				== typeid(CDataArray2D<unsigned char, 1> ).name()) {
			CDataArray2D<unsigned char, 1> *input =
					i_cPipelinePacket.getPayload<CDataArray2D<unsigned char, 1> >();

			setupFlagField(input);
			return;
		}

		std::cerr
				<< "ERROR: Video Output is only able to process (char,3) and (float,2) arrays, got "
				<< i_cPipelinePacket.type_info_name << std::endl;
		exit(-1);
	}

	/**
	 * callback from mainloop
	 */
	void main_loop_callback() {
		sdl_process_events();

		draw_scene();

		sdl_swap_buffers();
	}

	/**
	 * setup the flag field for visualization output based on flag input field
	 */
	void setupFlagField(CDataArray2D<unsigned char, 1> *i_input) {
		//printf("Flag field was updated...\n");
		output_cDataArray2D_uc4_flags.resize(i_input->width, i_input->height);

		for (int y = 0; y < output_cDataArray2D_uc4_flags.height; y++) {
			for (int x = 0; x < output_cDataArray2D_uc4_flags.width; x++) {
				unsigned char *d = &output_cDataArray2D_uc4_flags.getRef(x, y,
						0);

				switch (i_input->getRef(x, y)) {
				case 0:	// fluid
					d[0] = 255;
					d[1] = 255;
					d[2] = 255;
					d[3] = 0;
					break;

				case 1:	// obstacle
					d[0] = 0;
					d[1] = 0;
					d[2] = 0;
					d[3] = 255;
					break;

				case 2:	// inflow
					d[0] = 255;
					d[1] = 0;
					d[2] = 0;
					d[3] = 255;
					break;

				case 3:	// outflow
					d[0] = 0;
					d[1] = 0;
					d[2] = 255;
					d[3] = 255;
					break;

				default:
					std::cout << "Unknown flag " << (int) i_input->getRef(x, y)
							<< std::endl;
					break;
				}
			}
		}

		cGlTextureFlags.bind();

		if (cGlTextureFlags.width != i_input->width
				|| cGlTextureFlags.height != i_input->height)
			cGlTextureFlags.resize(i_input->width, i_input->height);

		cGlTextureFlags.setTextureParameters(
		GL_TEXTURE_2D,
		GL_RGBA,		// internal format
				GL_RGBA,
				GL_UNSIGNED_BYTE);

		cGlTextureFlags.setData(output_cDataArray2D_uc4_flags.data);

		cGlTextureFlags.unbind();
	}

public:
	void postprocess() {
		if (!cParameters.stage_fluidsimulation_visualize_flagfield)
			return;

		if (!output_cDataArray2D_uc4_flags.isValidData())
			return;

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		cGlTextureFlags.bind();

		glBegin(GL_QUADS);

		glTexCoord2f(0, 1);
		glVertex2f(-1, -1);
		glTexCoord2f(1, 1);
		glVertex2f(1, -1);
		glTexCoord2f(1, 0);
		glVertex2f(1, 1);
		glTexCoord2f(0, 0);
		glVertex2f(-1, 1);

		glEnd();

		cGlTextureFlags.unbind();
		glDisable(GL_BLEND);
	}

private:
	/**
	 * draw the scene to screen
	 */
	void draw_scene() {
		/*
		 // draw the color choice buttons
		 drawColorChoice();*/

		glClear(GL_COLOR_BUFFER_BIT);

		/*
		 * setup projection matrix to orthogonal
		 */
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-1, 1, -1, 1, -1, 1);

		/*
		 * setup model-view matrix to identity
		 */
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glEnable(GL_TEXTURE_2D);

		cGlTexture.bind();

		if (updateTexture) {
			if (output_cDataArray2D_uc3.isValidData()) {
				cGlTexture.setTextureParameters(
				GL_TEXTURE_2D,
				GL_RGB,		// internal format
						GL_RGB,
						GL_UNSIGNED_BYTE);

				if (cGlTexture.width != output_cDataArray2D_uc3.width
						|| cGlTexture.height != output_cDataArray2D_uc3.height)
					cGlTexture.resize(output_cDataArray2D_uc3.width,
							output_cDataArray2D_uc3.height);

				cGlTexture.setData(output_cDataArray2D_uc3.data);
			} else if (output_cDataArray2D_f2.isValidData()) {
				cGlTexture.setTextureParameters(
				GL_TEXTURE_2D,
				GL_RG32F,		// internal format
						GL_RG,
						GL_FLOAT);

				if (cGlTexture.width != output_cDataArray2D_f2.width
						|| cGlTexture.height != output_cDataArray2D_f2.height)
					cGlTexture.resize(output_cDataArray2D_f2.width,
							output_cDataArray2D_f2.height);

				cGlTexture.setData(output_cDataArray2D_f2.data);
			} else if (output_cDataArray2D_f3.isValidData()) {
				cGlTexture.setTextureParameters(
				GL_TEXTURE_2D,
				GL_RGB32F,		// internal format
						GL_RGB,
						GL_FLOAT);

				if (cGlTexture.width != output_cDataArray2D_f3.width
						|| cGlTexture.height != output_cDataArray2D_f3.height)
					cGlTexture.resize(output_cDataArray2D_f3.width,
							output_cDataArray2D_f3.height);

				cGlTexture.setData(output_cDataArray2D_f3.data);
			} else {
				std::cout << "nothing to draw" << std::endl;
				cGlTexture.unbind();
				return;
			}

			updateTexture = false;
		}

		glBegin(GL_QUADS);

		glTexCoord2f(0, 1);
		glVertex2f(-1, -1);
		glTexCoord2f(1, 1);
		glVertex2f(1, -1);
		glTexCoord2f(1, 0);
		glVertex2f(1, 1);
		glTexCoord2f(0, 0);
		glVertex2f(-1, 1);

		glEnd();

		cGlTexture.unbind();

		overlayText();

		postprocess();

	}

private:
	void overlayText() {
		for (std::list<Text*>::iterator it = Texts->begin(); it != Texts->end();
				++it) {
			if ((*it)->ttl == 0) {
				continue;
			}
			float old_size = F->getFontSize();
			F->setFontSize((*it)->size);
			F->printfxy(((*it)->x) * width / output_cDataArray2D_uc3.width,
					((*it)->y) * height / output_cDataArray2D_uc3.height,
					(*it)->msg.c_str());
			F->setFontSize(old_size);
			(*it)->ttl--;
		}

		//Print Points
		if (Points)
			F->printfxy(width * 0.3, height * 0.93, Points->c_str());

		//Print Remaining Time
		if (Rest_time)
			F->printfxy(width * 0.5, height * 0.93, Rest_time->c_str());

		/*//Print the "Choose Color" text to the bottom of the window
		 std::string const test_text = "Choose color:";
		 F->printfxy(210,height*0.93,test_text.c_str());*/
	}

public:
	/**
	 * callback from SDL for key down event
	 */
	void sdl_handle_key_down(SDL_keysym &i_keysym) {
		char key;
		switch (i_keysym.sym) {
		case SDLK_q:
			cParameters.exit = true;
			break;

		default:
			key = i_keysym.sym;
			cParameters.key_down(key);
#ifndef KINECT_INCLUDED
			if (cStage_KeyboardForces) { //Checking for NULL-Pointer
				cStage_KeyboardForces->keyDown(i_keysym.sym);
			}
#endif
			break;
		}
	}

	void sdl_handle_key_up(SDL_keysym &i_keysym) {
		switch (i_keysym.sym) {
		default:
#ifndef KINECT_INCLUDED
			if (cStage_KeyboardForces) { //Checking for NULL-Pointer
				cStage_KeyboardForces->keyUp(i_keysym.sym);
			}
#endif
			break;
		}
	}

	inline void pipeline_push_drawing(int i_mouse_x, int i_mouse_y,
			unsigned char r, unsigned char g, unsigned char b) {
		CDataDrawingInformation d;
		d.x = (double) i_mouse_x / (double) width;
		d.y =
				(double) i_mouse_y
						/ (double) (height
								- 30
										* (height
												/ (double) output_cDataArray2D_uc3.height));
		d.size = 0.1;
		d.color.r = r;
		d.color.g = g;
		d.color.b = b;

		pipeline_push(d);
	}

	void interactiveStuff(int i_mouse_buttons, int i_mouse_x, int i_mouse_y) {
		// check if the left mouse button was pressed
		if (i_mouse_buttons & 2) {
			// Check if the mouse is in zone at the bottom of the window
			if ((double) i_mouse_y
					/ (double) (height
							- 30
									* (height
											/ (double) output_cDataArray2D_uc3.height))
					> 1) {
				struct timeval now;
				gettimeofday(&now, NULL);
				float difftime_exact = (now.tv_sec - last_button_pressed.tv_sec)
						* 1000.0
						+ (now.tv_usec - last_button_pressed.tv_usec) / 1000.0;

				std::stringstream sst;
				sst.flush();
				sst << difftime_exact;
				std::string s = sst.str();

				if (difftime_exact < 500) {
					// if the last button press is less than 0.5 sec ago, then do nothing
					printf(
							"You cannot press the button because it is to early again. Elapsed time: ");
					printf(s.c_str());
					printf("\n");
					return;
				} else {
					// else, update the time when the button was pressed for the last time
					// and go on with the action, that the button envoked
					gettimeofday(&last_button_pressed, NULL);
				}

				// If it is, check, if a new color was chosen
				double scaledX = (double) output_cDataArray2D_uc3.width
						* (double) i_mouse_x / (double) width; //calculate the normal positions because the window may have been scaled (e.g. fullscreen)
				double scaledY =
						(double) output_cDataArray2D_uc3.height
								* (double) i_mouse_y
								/ (double) (height
										- 30
												* (height
														/ (double) output_cDataArray2D_uc3.height));
				if ((output_cDataArray2D_uc3.height + 5 <= scaledY)
						&& (output_cDataArray2D_uc3.height + 25 >= scaledY)) {
					// red
					if ((scaledX >= 10) && (scaledX <= 30)) {
						activeColor.r = 255;
						activeColor.g = 0;
						activeColor.b = 0;
					}
					// blue
					else if ((scaledX >= 40) && (scaledX <= 60)) {
						activeColor.r = 0;
						activeColor.g = 0;
						activeColor.b = 255;
					}
					// green
					else if ((scaledX >= 70) && (scaledX <= 90)) {
						activeColor.r = 0;
						activeColor.g = 255;
						activeColor.b = 0;
					}
					// black
					else if ((scaledX >= 100) && (scaledX <= 120)) {
						activeColor.r = 0;
						activeColor.g = 0;
						activeColor.b = 0;
					}
					// white (the eraser)
					else if ((scaledX >= 130) && (scaledX <= 150)) {
						activeColor.r = 255;
						activeColor.g = 255;
						activeColor.b = 255;
					}
					// start button
					else if ((scaledX >= 470) && (scaledX <= 530)) {
						if (!cParameters.pause) {
							printf(
									"The start button was pressed, but the game is already running...\n");
						} else {
							printf(
									"The start button was pressed and the game will continue.\n");
						}
						//cParameters.pause = !cParameters.pause;
						cParameters.pause = !cParameters.pause;
					}
					// save the image
					else if ((scaledX >= 550) && (scaledX <= 610)) {
						std::stringstream ss;
						ss.flush();
						ss << "data/img-saved_";
						time_t rawtime;
						struct tm * timeinfo;
						time(&rawtime);
						timeinfo = localtime(&rawtime);
						ss << timeinfo->tm_year + 1900 << "-"
								<< timeinfo->tm_mon + 1 << "-"
								<< timeinfo->tm_mday << "_";
						ss << timeinfo->tm_hour << "-" << timeinfo->tm_min
								<< "-" << timeinfo->tm_sec;
						ss << ".bmp";
						std::string os_path = ss.str();
						/*char* pixels = reinterpret_cast<char*>(output_cDataArray2D_uc3.data);
						 SDL_Surface *surface = SDL_CreateRGBSurfaceFrom(pixels,
						 output_cDataArray2D_uc3.width,
						 output_cDataArray2D_uc3.height,
						 3*8, // channels*8, the bits per pixel
						 output_cDataArray2D_uc3.width*3, // pitch
						 0x0000FF, // red mask
						 0x00FF00, // green mask
						 0xFF0000, // blue mask
						 0); // alpha mask (none)
						 printf(reinterpret_cast<const char*>(output_cDataArray2D_uc3_flags.width));
						 printf(reinterpret_cast<const char*>(output_cDataArray2D_uc3_flags.height));
						 char* pixels = reinterpret_cast<char*>(output_cDataArray2D_uc3_flags.data);
						 SDL_Surface *surface = SDL_CreateRGBSurfaceFrom(pixels,
						 output_cDataArray2D_uc3_flags.width,
						 output_cDataArray2D_uc3_flags.height,
						 3*8, // channels*8, the bits per pixel
						 output_cDataArray2D_uc3_flags.width*3, // pitch = width*channel
						 0x0000FF, // red mask
						 0x00FF00, // green mask
						 0xFF0000, // blue mask
						 0); // alpha mask (none)
						 SDL_SaveBMP(surface, os_path);*/
						cStage_Visualization->save_flag_field_as_image(
								os_path.c_str());
						printf("The image was saved to the path '");
						printf("%s", os_path.c_str());
						printf("'.\n");
					}

				}

			} else { // If its not, draw in the window

				/*
				 // left+right
				 if (i_mouse_buttons == 10)
				 {
				 // clear
				 pipeline_push_drawing(i_mouse_x, i_mouse_y, 255, 255, 255);
				 return;
				 }

				 // left
				 if (i_mouse_buttons & 2)
				 pipeline_push_drawing(i_mouse_x, i_mouse_y, 255, 0, 0);

				 // right
				 if (i_mouse_buttons & 8)
				 pipeline_push_drawing(i_mouse_x, i_mouse_y, 0, 0, 255);

				 // middle
				 if (i_mouse_buttons & 4)
				 pipeline_push_drawing(i_mouse_x, i_mouse_y, 0, 0, 0);
				 */

				// you can draw just with the left mouse button,
				// the color is defined by "activeColor",
				if (i_mouse_buttons & 2)
					pipeline_push_drawing(i_mouse_x, i_mouse_y, activeColor.r,
							activeColor.g, activeColor.b);

			}
		}
	}

	void sdl_handle_mouse_button_down(int /*i_button_id*/, int i_mouse_buttons,
			int i_mouse_x, int i_mouse_y) {
		interactiveStuff(i_mouse_buttons, i_mouse_x, i_mouse_y);
	}

	void sdl_handle_mouse_button_up(int /*i_button_id*/, int i_mouse_buttons,
			int i_mouse_x, int i_mouse_y) {
		interactiveStuff(i_mouse_buttons, i_mouse_x, i_mouse_y);
	}

	void sdl_handle_mouse_motion(int i_mouse_buttons, int i_mouse_x,
			int i_mouse_y) {
		interactiveStuff(i_mouse_buttons, i_mouse_x, i_mouse_y);
	}

	virtual void sdl_handle_quit() {
		cParameters.exit = true;
	}

	void push_Text(std::string msg, int ttl, int x, int y, int size) {
		Text* txt = new Text();
		txt->msg = msg;
		txt->ttl = ttl;
		txt->x = x;
		txt->y = y;
		txt->size = size;
		Texts->push_back(txt);
	}

	void change_Points(std::string *newPoints) {
		Points = newPoints;
	}

	void change_Time(std::string *newTime) {
		Rest_time = newTime;
	}

	/*inline void color_point_unsafe(int x, int y, unsigned char r, unsigned char g, unsigned char b)
	 {
	 output_cDataArray2D_uc3.getRef(x, y, 0) = r;
	 output_cDataArray2D_uc3.getRef(x, y, 1) = g;
	 output_cDataArray2D_uc3.getRef(x, y, 2) = b;
	 }

	 void draw_point_unsafe(int x, int y, unsigned char r, unsigned char g, unsigned char b, int xwidth, int ywidth)
	 {
	 for (int dx=-floor(xwidth/2.0); dx<=ceil(xwidth/2.0); dx++)
	 for (int dy=-floor(ywidth/2.0); dy<=ceil(ywidth/2.0); dy++)
	 {
	 color_point_unsafe(x+dx, y+dy, r, g, b);
	 }
	 }

	 void drawColorChoice()
	 {
	 for (int x=1; x<=width; x++){
	 for(int y=height-30; y<=height; y++){
	 color_point_unsafe(x,y,50,50,50);
	 }
	 } // Grey Background
	 draw_point_unsafe(360,(int)(height-15),255,0,0,20,20); // Red
	 draw_point_unsafe(390,(int)(height-15),0,0,255,20,20); // Blue
	 draw_point_unsafe(420,(int)(height-15),0,255,0,20,20); // Green
	 draw_point_unsafe(450,(int)(height-15),0,0,0,20,20); // Black
	 draw_point_unsafe(480,(int)(height-15),255,255,255,20,20); // White (for the Eraser)
	 }*/

#ifndef KINECT_INCLUDED
	void setupCStageKeyboardForces(CStage_KeyboardForces* cKeyForce) {
		cStage_KeyboardForces = cKeyForce;
	}
#endif
};

#endif /* CSTAGE_VIDEOOUTPUT_HPP_ */
