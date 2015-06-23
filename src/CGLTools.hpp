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

/**
 * version 0.2
 *
 * history:
 *
 * 2009-05-23:
 * - update to class
 *
 * 2007-09-02:
 * - removed loadTexture because this is done by CTexture.hpp
 * - renamed WaitFPS to limit_fps
 */

#ifndef CGL_TOOLS_H__
#define CGL_TOOLS_H__

#include <GL/gl.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * check for gl error and output the line and function where the error occured
 */

class CGLTools {
public:
	/**
	 * restore the frustum from a matrix and return the values to parameters
	 */
	static bool unFrustum(float *matrix, float *left, float *right,
			float *bottom, float *top, float *near, float *far);

	/**
	 * check the opengl error state variable and output some information if error exists
	 */
	static void checkGLError_(const char *file, const char *function, int line);

};

#if CHECK_GL_ERROR_OFF
#define checkGlError()
#else
#define GL_checkError()	CGLTools::checkGLError_(__FILE__, __FUNCTION__, __LINE__)
#endif

#endif
