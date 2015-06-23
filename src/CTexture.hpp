/*
 * Copyright 2010 Martin Schreiber
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TEXTURE_H
#define TEXTURE_H
/**
 * version 0.2
 *
 * changes:
 * 2007-09-27:
 * 	added proxy texture to detect memory overrun
 */

#include <SDL_image.h>
#include "CGLTools.hpp"
#include <iostream>

using namespace std;

class CTexture {
public:
	bool valid;
	int width, height;

	GLuint textureid;

	GLint int_format;
	GLenum ext_format;
	GLenum ext_type;

	GLenum target;	// usually GL_TEXTURE_2D

	/**
	 * create empty texture object
	 */
	inline CTexture() {
		valid = false;
	}

	/**
	 * create empty texture but create gl texture handle with given texture size
	 */
	inline CTexture(GLint p_width, GLint p_height, GLint p_int_format = GL_RGBA,
			GLenum p_ext_format = GL_RGBA, GLenum p_ext_type = GL_UNSIGNED_BYTE,
			GLenum p_target = GL_TEXTURE_2D) :
			width(p_width), height(p_height), int_format(p_int_format), ext_format(
					p_ext_format), ext_type(p_ext_type), target(p_target) {
		// create one texture name
		glGenTextures(1, &textureid);

		bind();

		if (p_target == GL_TEXTURE_2D)	// only works with GL_TEXTURE_2D yet
		{
			// try to create proxy texture
			glTexImage2D(GL_PROXY_TEXTURE_2D, 0, int_format, width, height, 0,
					ext_format, ext_type, 0);

			GLint x, y;
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,
					&x);
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,
					&y);

			if (x == 0 && y == 0) {
				cerr
						<< "ERROR: failed to create texture (maybe not enough memory)"
						<< endl;
				glDeleteTextures(1, &textureid);
				valid = false;
				return;
			}

		}

		// this reads from the sdl surface and puts it into an opengl texture
		glTexImage2D(target, 0, int_format, width, height, 0, ext_format,
				ext_type, 0);
		GL_checkError();

		if (p_target == GL_TEXTURE_2D) {
			setParamInt(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			setParamInt(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			setParamInt(GL_TEXTURE_WRAP_S, GL_CLAMP);
			setParamInt(GL_TEXTURE_WRAP_T, GL_CLAMP);
		} else {
			setParamInt(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			setParamInt(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			setParamInt(GL_TEXTURE_WRAP_S, GL_CLAMP);
			setParamInt(GL_TEXTURE_WRAP_T, GL_CLAMP);
		}
		GL_checkError();

		valid = true;
	}

	inline void loadFromFile(const char *filename) {
		cleanup();

		SDL_Surface *surface = IMG_Load(filename);
		target = GL_TEXTURE_2D;

		// could not load filename
		if (!surface) {
			valid = false;
			cerr << "cannot load image" << endl;
			return;
		}

		// work out what format to tell glTexImage2D to use...
		if (surface->format->BytesPerPixel == 3) { // RGB 24bit
			int_format = GL_RGB;
		} else if (surface->format->BytesPerPixel == 4) { // RGBA 32bit
			int_format = GL_RGBA;
		} else {
			SDL_FreeSurface(surface);
			cerr << "image has to be 3 or 4 bytes per pixel" << endl;
			valid = false;
			return;
		}

		width = surface->w;
		height = surface->h;

		// create one texture name
		glGenTextures(1, &textureid);

		bind();

		ext_format = int_format;
		ext_type = GL_UNSIGNED_BYTE;
		// this reads from the sdl surface and puts it into an opengl texture
		glTexImage2D(target, 0, int_format, surface->w, surface->h, 0,
				ext_format, ext_type, surface->pixels);
		SDL_FreeSurface(surface);

		setParamInt(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		setParamInt(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		setParamInt(GL_TEXTURE_WRAP_S, GL_REPEAT);
		setParamInt(GL_TEXTURE_WRAP_T, GL_REPEAT);

		valid = true;
	}

	/**
	 * create texture from file
	 */
	inline CTexture(const char *filename) {
		valid = false;
		loadFromFile(filename);
	}

	/**
	 * bind texture to texture unit
	 */
	inline void bind() {
		glBindTexture(target, textureid);
	}

	inline void unbind() {
		glBindTexture(target, 0);
	}

	inline void setParamFloat(GLenum pname, GLfloat param) {
		glTexParameterf(target, pname, param);
	}

	inline void setParamInt(GLenum pname, GLint param) {
		glTexParameteri(target, pname, param);
	}

	inline GLfloat getParami(GLenum pname) {
		GLfloat param;
		glGetTexParameterfv(target, pname, &param);
		return param;
	}

	inline GLint getParam(GLenum pname) {
		GLint param;
		glGetTexParameteriv(target, pname, &param);
		return param;
	}

	inline GLuint getTexture() {
		return textureid;
	}

	inline GLenum getFormat() {
		return int_format;
	}

	/**
	 * copy the data pointed by data to the texture
	 *
	 * \param data	pointer to new texture data
	 */
	inline void setData(void *data) {
		bind();
		glTexSubImage2D(target, 0, 0, 0, width, height, ext_format, ext_type,
				data);
		unbind();
	}

	/**
	 * copy the texture data to the memory location referred by data
	 * 
	 * \param data	pointer to store the texture data
	 */
	inline void getData(void *data) {
		bind();
		glGetTexImage(target, 0, ext_format, ext_type, data);
		unbind();
	}

	/**
	 * calculate the difference of the red channel from this texture to texture tex
	 * \param tex		texture to calculate difference
	 * \param channel	offset of color channel (r=0, g=1, b=2, a=3)
	 */
	inline float getDiff(CTexture *tex, int channel = 0) {
		if (this->width != tex->width || this->height != tex->height) {
			fprintf(stderr, "width / height doesnt match\n");
			return -1.0;
		}

		float *data1 = new float[width * height * 4];
		float *data2 = new float[width * height * 4];
		this->getData(data1);
		this->getData(data2);

		float diff = 0;
		float *datap1 = data1 + channel;
		float *datap2 = data2 + channel;
		for (int i = 0; i < width * height; i++) {
			diff = (*data1 - *data2) * (*data1 - *data2);
			datap1 += 4;
			datap2 += 4;
		}

		delete data1;
		delete data2;

		return diff;
	}

	/**
	 * setup quickly nearest neighbor filtering
	 */
	inline void setupNearestNeighbor() {
		bind();
		setParamInt(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		setParamInt(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		unbind();
	}

	inline void cleanup() {
		if (valid) {
			glDeleteTextures(1, &textureid);
			valid = false;
		}
	}

	inline ~CTexture() {
		cleanup();
	}
};

#endif
