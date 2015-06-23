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

#ifndef __FONT_H
#define __FONT_H

#include "CTexture.hpp"
#include <SDL.h>
#include <GL/glu.h>
#include <string.h>
#include <stdarg.h>
#include "CFont.hpp"
#include <iostream>
#include <string>

class CFont {
private:
	CTexture *fontTexture;
	GLuint fontDisplayList;
	float charWidth, charHeight;
	float scaleX, scaleY;

public:
	/**
	 * init successful
	 */
	int init;

	CFont(const string &file = "data/fonts/fontmap.png") :
			fontTexture(NULL) {
		setup(file);
	}

	void setup(const string &file = "data/fonts/fontmap.png") {
		if (fontTexture != NULL)
			delete fontTexture;

		fontTexture = new CTexture(file.c_str());

		if (fontTexture->valid == false) {
			std::cerr << "failed to load fontfile '" << file << "'"
					<< std::endl;
			init = false;
			return;
		}

		charWidth = (float) fontTexture->width / 16.0f;
		charHeight = (float) fontTexture->height / 16.0f;

		setFontSize(charHeight);

		fontDisplayList = glGenLists(256);
		if (fontDisplayList == 0) {
			std::cout << "fontdl failed" << std::endl;
			init = false;
			return;
		}

		init = true;

		fontTexture->bind();

		float x, y;
		for (int pos = 0; pos < 256; pos++) {
			x = ((float) (pos & 0x0f)) / 16.0f;
			y = ((float) (pos >> 4)) / 16.0f;

			glNewList(fontDisplayList + pos, GL_COMPILE);

			glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);
			glEnable(GL_BLEND);
			glBegin(GL_QUADS);
			glTexCoord2f(x, y);
			glVertex2d(0, 1);	// top left
			glTexCoord2f(x, y + (1.0f / 16.0f));
			glVertex2d(0, 0);	// bottom left
			glTexCoord2f(x + (1.0f / 16.0f), y + (1.0f / 16.0f));
			glVertex2d(-1, 0);	// bottom right
			glTexCoord2f(x + (1.0f / 16.0f), y);
			glVertex2d(-1, 1);	// top right
			glEnd();
			//		glTranslatef(-charWidth, 0, 0);
			glTranslatef(-1, 0, 0);

			glEndList();
		}
	}

	void setFontSize(float size) {
		scaleY = size;
		scaleX = charWidth * (scaleY / charHeight);
	}

	float getFontSize() {
		return scaleY;
	}

	void printfxy(int x, int y, const char *format, ...) {
		va_list list;
		static char text[2048];
		float viewport[4];
		GLuint EDepth, EBlend, ETexture, ELight;

		if (!init)
			return;

		va_start(list, format);
		vsnprintf(text, 2048, format, list);

		EDepth = glIsEnabled(GL_DEPTH_TEST);
		EBlend = glIsEnabled(GL_BLEND);
		ETexture = glIsEnabled(GL_TEXTURE_2D);
		ELight = glIsEnabled(GL_LIGHTING);

		glDisable(GL_DEPTH_TEST);
		glEnable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);

		//	glEnable(GL_BLEND);
		//	glBlendEquation(GL_FUNC_ADD);
		//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		glGetFloatv(GL_VIEWPORT, viewport);
		glOrtho(viewport[0], -viewport[2], viewport[1], viewport[3], -1.0, 1.0);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		fontTexture->bind();

		glTranslatef(-x, viewport[3] - y - charHeight, 0.0);
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		glMatrixMode(GL_MODELVIEW);	// for glTranslate in displaylist
		glScalef(scaleX, scaleY, 1);

		glListBase(fontDisplayList - 32);
		glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
		glListBase(0);

		fontTexture->unbind();

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		if (EDepth)
			glEnable(GL_DEPTH_TEST);
		if (EBlend)
			glEnable(GL_BLEND);
		if (!ETexture)
			glDisable(GL_TEXTURE_2D);
		if (ELight)
			glEnable(GL_LIGHTING);
	}
};

#endif
