/*
 * CStage_ImageInput.hpp
 *
 *  Created on: Jul 5, 2013
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CSTAGE_VIDEOINPUT_APPLE_HPP_
#define CSTAGE_VIDEOINPUT_APPLE_HPP_

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "SDL/SDL_image.h"
#include "CDataArray2D.hpp"
#include "CDataDrawingInformation.hpp"
#include <cv.h>
#include <cassert>

/**
 * class providing static image input
 *
 * the image is send to the pipeline during each main loop iteration
 */
class CStage_VideoInput	:	public
	CPipelineStage,
	ImageReader
{
	/**
	 * global parameters
	 */
	CParameters &cParameters;

	/**
	 * storage of image and packet for pipelining
	 */
	CDataArray2D<unsigned char,3> cDataArray2D;
	
	/**
	 * OpenCV capture device and frame
	 */
	cv::VideoCapture* capture;
	cv::Mat frame;
	
	int image_width, image_height


public:
	/**
	 * constructor
	 */
	CStage_VideoInput(CParameters &i_cParameters)	:
		CPipelineStage("ImageInput"),
		image_width(-1),
		image_height(-1),
		cParameters(i_cParameters)
	{
		capture = new cv::VideoCapture(-1);
		
		image_height = cParameters.input_video_height;
		image_width = cParameters.input_video_width;
		capture->set(CV_CAP_PROP_FRAME_HEIGHT, image_height);
		capture->set(CV_CAP_PROP_FRAME_WIDTH, image_width);

		if (cParameters.verbosity_level > 5)
			std::cout << "Using input video size of " << image_width << " x " << image_height << std::endl;

		// dummy data: blue color
		cDataArray2D.resize(image_width, image_height);
		unsigned char *d = cDataArray2D.data;
		for (int i = 0; i < image_width*image_height*3; i++)
			d[i] = 100;

	}

	void main_loop_callback()
	{
		capture.read(frame);
		memcpy(cDataArray2D.data, frame.data, image_width*image_height*3);
		
		pipeline_push();
	}

	virtual ~CStage_VideoInput()
	{
		std::cout << "CV: Releasing capturing-device" << std::endl;
		capture->release();

		std::cout << "CV: Killing capturing-device" << std::endl;
	    delete capture;
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push()
	{
		CPipelineStage::pipeline_push((CPipelinePacket&)cDataArray2D);
	}



private:
	/**
	 * draw a circle into an image
	 */
	void draw_circle(
			CDataDrawingInformation &i_d
	)
	{
		int i_x = i_d.x*(double)cDataArray2D.width;
		int i_y = i_d.y*(double)cDataArray2D.height;

		int i_size = i_d.size*(double)cDataArray2D.width;

		for (int ry = -i_size; ry < i_size; ry++)
		{
			int y = i_y + ry;

			if (y < 0 || y >= cDataArray2D.height)
					continue;

			for (int rx = -i_size; rx < i_size; rx++)
			{
				if (rx*rx + ry*ry > i_size)
					continue;

				int x = i_x + rx;

				if (x < 0 || x >= cDataArray2D.width)
						continue;

				unsigned char *d = &cDataArray2D.getRef(x, y);
				d[0] = i_d.color.r;
				d[1] = i_d.color.g;
				d[2] = i_d.color.b;
			}
		}
	}

public:
	/**
	 * process incoming pipeline input.
	 *
	 * the only input processed so far is from the video output stage to
	 * draw something into the image.
	 */
	void pipeline_process_input(
			CPipelinePacket &i_cPipelinePacket
	)
	{
		if (i_cPipelinePacket.type_info_name != typeid(CDataDrawingInformation).name())
		{
			std::cerr << "ERROR: Only data drawing information allowed for image input stage" << std::endl;
			exit(-1);
		}

		// unpack
		CDataDrawingInformation *d = i_cPipelinePacket.getPayload<CDataDrawingInformation>();

		// handle package by drawing a circle
		draw_circle(*d);

		// push modifications on image to pipeline
		pipeline_push();
	}
};

#endif /* CSTAGE_IMAGEINPUT_APPLE_HPP_ */
