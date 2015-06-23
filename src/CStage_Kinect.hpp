#ifndef CSTAGE_KINECT_HPP_
#define CSTAGE_KINECT_HPP_

#include "CParameters.hpp"
#include "CPipelineStage.hpp"
#include "CDataVector2.hpp"
#include <XnOpenNI.h>
#include <XnCppWrapper.h>
#include <exception>
#include <XnCyclicStackT.h>
#include <XnHashT.h>

#if defined(GAME2)
#include "Game_version2.hpp"
#else
#include "Game.hpp"
#endif

#define SAMPLE_XML_FILE_LOCAL "Kinect.xml"

#define CHECK_RC(rc, what)							\
    if (rc != XN_STATUS_OK)							\
{									\
    printf("%s failed: %s\n", what, xnGetStatusString(rc));	\
    return rc;							\
}

class CStage_Kinect: public CPipelineStage {
	/**
	 * global parameters
	 */
	CParameters &cParameters;

	xn::Context context; // OpenNI.
	xn::EnumerationErrors errors;	// errors
	xn::ScriptNode context_scriptNode;

	xn::ScriptNode scriptNode;
	static xn::HandsGenerator g_HandsGenerator;
	xn::GestureGenerator g_GestureGenerator;
	CDataVector2* cVector;
	CDataVector2* old_output_vector;
	CDataVector2* new_output_vector;
	//Declaration ZeroPoint
	static XnPoint3D zero;
	static XnPoint3D bufferedPoint;

	static int flag_hand;
	static int flag_valid_data;

	static Game* game;

public:
	/**
	 * constructor
	 */
	CStage_Kinect(CParameters &i_cParameters) :
			CPipelineStage("Kinect"), cParameters(i_cParameters) {
		cVector = new CDataVector2();
		old_output_vector = new CDataVector2();
		new_output_vector = new CDataVector2();
		XnStatus rc = context.InitFromXmlFile(SAMPLE_XML_FILE_LOCAL,
				context_scriptNode, &errors);
		if (rc == XN_STATUS_NO_NODE_PRESENT) {
			XnChar strError[1024];
			errors.ToString(strError, 1024);
			printf("ERROR: InitFromXmlFile - %s\n", strError);
			exit(-1);
		} else if (rc != XN_STATUS_OK) {
			printf("Open failed: %s\n", xnGetStatusString(rc));
			exit(-1);
		}

		rc = context.FindExistingNode(XN_NODE_TYPE_HANDS, g_HandsGenerator);

		if (rc != XN_STATUS_OK) {
			std::cerr << "FindExistingNode" << std::endl;
			exit(-1);
		}

		//dummy = CHECK_RC(rc, "Find hands generator");
		rc = context.FindExistingNode(XN_NODE_TYPE_GESTURE, g_GestureGenerator);

		if (rc != XN_STATUS_OK) {
			std::cerr << "FindExistingNode" << std::endl;
			exit(-1);
		}

		//dummy = CHECK_RC(rc, "Find gest generator");
		rc = context.StartGeneratingAll();

		if (rc != XN_STATUS_OK) {
			std::cerr << "StartGeneratingAll" << std::endl;
			exit(-1);
		}

		//Gestures to recognize

		rc = g_GestureGenerator.AddGesture("Click", NULL);

		if (rc != XN_STATUS_OK) {
			std::cerr << "AddGesture" << std::endl;
			exit(-1);
		}
		/*
		 rc = g_GestureGenerator.AddGesture("HandRaised", NULL);

		 if (rc != XN_STATUS_OK)
		 {
		 std::cerr << "AddGesture Swipe" << std::endl;
		 exit(-1);
		 }
		 */
		XnCallbackHandle chandle;

		rc = g_HandsGenerator.RegisterHandCallbacks(HandCreate, Hand_Update,
				Hand_Destroy, this, chandle);
		//CHECK_RC(rc, "Register Hand Callback");

		rc = g_GestureGenerator.RegisterGestureCallbacks(Gesture_Recognized,
				Gesture_Process, this, chandle);
		//CHECK_RC(rc, "Register Gest Callback");

		// CStage_Kinect::game->global_pause = true;
	}

	void main_loop_callback() {
		context.WaitAndUpdateAll();
		new_output_vector = calculateVector();
		if (CStage_Kinect::flag_valid_data == 1) {
			pipeline_push(new_output_vector);
			old_output_vector = new_output_vector;
			//printf("valid vector is %f  and %f  \n", new_output_vector->data->x,new_output_vector->data->y);
		} else {
			pipeline_push(old_output_vector);
			//printf("invalid vector is %f  and %f  \n", old_output_vector->data->x,old_output_vector->data->y);
		}

	}

	//---------------------------------------------------------------------------------------------------------------------------------------------
	//CALLBACKS
	//---------------------------------------------------------------------------------------------------------------------------------------------

	static void XN_CALLBACK_TYPE SessionProgress(const XnChar* /*strFocus*/,
			const XnPoint3D& /*ptFocusPoint*/, XnFloat /*fProgress*/,
			void* /*UserCxt*/) {
		//printf("Session progress (%6.2f,%6.2f,%6.2f) - %6.2f [%s]\n", ptFocusPoint.X, ptFocusPoint.Y, ptFocusPoint.Z, fProgress,  strFocus);
	}
	// Callback for session start
	static void XN_CALLBACK_TYPE SessionStart(const XnPoint3D& /*ptFocusPoint*/,
			void* /*UserCxt*/) {
		//printf("Session started. Please push or swipe (%6.2f,%6.2f,%6.2f)...\n", ptFocusPoint.X, ptFocusPoint.Y, ptFocusPoint.Z);
	}
	// Callback for session end
	static void XN_CALLBACK_TYPE SessionEnd(void* /*UserCxt*/) {
		//printf("Session ended. Please perform focus gesture to start session\n");
	}
	// Callback for hand created
	static void XN_CALLBACK_TYPE HandCreate(xn::HandsGenerator& /*generator*/,
			XnUserID /*user*/, const XnPoint3D* /*pPosition*/,
			XnFloat /*fTime*/, void* /*pCookie*/) {
		//printf("New Hand: %d @ (%f,%f,%f)\n", user, pPosition->X, pPosition->Y, pPosition->Z);
	}
	// Callback for hand updated
	static void XN_CALLBACK_TYPE Hand_Update(xn::HandsGenerator& /*generator*/,
			XnUserID /*nId*/, const XnPoint3D* pPosition, XnFloat /*fTime*/,
			void* /*pCookie*/) {
		//printing the position of the hand
		//printf("Hand Update: %d @ (%f,%f,%f)\n", nId, pPosition->X, pPosition->Y, pPosition->Z);
		/*
		 float xdiff, ydiff, scale;
		 //value = sqrt((pPosition->X - zero.X)*(pPosition->X - zero.X) + (pPosition->Y - zero.Y)*(pPosition->Y - zero.Y));
		 scale = 850;
		 xdiff = (zero.X - pPosition->X)/scale;
		 ydiff = (pPosition->Y - zero.Y)/scale;
		 value = sqrt(ydiff*ydiff + xdiff*xdiff);
		 printf("magnitude is %f  \n", value);
		 */
		int checkVal = 60;
		if ((pPosition->X - bufferedPoint.X) > checkVal
				|| (pPosition->Y - bufferedPoint.Y) > checkVal
				|| (pPosition->Z - bufferedPoint.Z) > checkVal) {
			printf("Buffered:  X=%f, Y=%f, Z=%f \n", bufferedPoint.X,
					bufferedPoint.Y, bufferedPoint.Z);
			printf("pPosition: X=%f, Y=%f, Z=%f)\n\n", pPosition->X,
					pPosition->Y, pPosition->Z);
		}
		CStage_Kinect::bufferedPoint = *pPosition;

	}
	//Callback for hand lost
	static void XN_CALLBACK_TYPE Hand_Destroy(xn::HandsGenerator& /*generator*/,
			XnUserID nId, XnFloat /*fTime*/, void* /*pCookie*/) {
		printf("Lost Hand: %d\n", nId);
		//CStage_Kinect::game->global_pause = true;
		CStage_Kinect::flag_hand = 0;
	}
	static void XN_CALLBACK_TYPE Gesture_Process(
			xn::GestureGenerator& /*generator*/, const XnChar* /*strGesture*/,
			const XnPoint3D* /*pPosition*/, XnFloat /*fProgress*/,
			void* /*pCookie*/) {
	}
	//Callback for gesture recognizing
	static void XN_CALLBACK_TYPE Gesture_Recognized(
			xn::GestureGenerator& /*generator*/, const XnChar* strGesture,
			const XnPoint3D* /*pIDPosition*/, const XnPoint3D* pEndPosition,
			void* /*pCookie*/) {
		//std::cout << "GESTURE RECOGNIZED" << std::endl;
		printf("Geste: %s", strGesture);
		if (CStage_Kinect::flag_hand == 0) { 				//just track 1 hand
			//printf("Click gesture recognized kinect \n");
			//printf("value of flag is %d \n", flag_hand);
			CStage_Kinect::flag_hand = 1;
			//CStage_Kinect::game->global_pause = false;

			//set ZeroPoint
			CStage_Kinect::zero = *pEndPosition;
			//printf("Zero changed to: x=%f y=%f\n", CStage_Kinect::zero.X, CStage_Kinect::zero.Y);
			bufferedPoint = *pEndPosition;
			//printf("ZeroPoint: (%f, %f, %f) \n", zero->X, zero->Y, zero->Z);
			g_HandsGenerator.StartTracking(*pEndPosition);
		} else {
			//printf("value of flag is %d \n", flag_hand);
		}

	}

	void pipeline_process_input(CPipelinePacket &i_cPipelinePacket) {
	}

public:
	/**
	 * manually triggered pushing of next image to the pipeline
	 */
	void pipeline_push(CDataVector2* output_cVector) {
		CPipelineStage::pipeline_push((CPipelinePacket&) *output_cVector);
		//printf("Pushed \n");
	}

public:
	CDataVector2* calculateVector() {
		//CDataVector2 cVector;
		double scale = 650.0;
		CStage_Kinect::flag_valid_data = 1;
		if (CStage_Kinect::flag_hand == 1) {
			cVector->data->x = (-1) * (bufferedPoint.X - zero.X) / scale;
			//printf("ZERO: x=%f, y=%f",zero.X, zero.Y);
			//printf("x val %f, bufX %f, y val %f, bufY %f", ((-1)*(bufferedPoint.X - zero.X)/scale) , bufferedPoint.X, ((-1)*(bufferedPoint.Y - zero.Y)/scale),bufferedPoint.Y);
			cVector->data->y = (-1) * (bufferedPoint.Y - zero.Y) / scale;
			if (cVector->data->x < 0.05 && cVector->data->x > -0.05) {
				cVector->data->x = 0;
			}
			if (cVector->data->y < 0.1 && cVector->data->y > -0.1) {
				cVector->data->y = 0;
			}

		} else {
			cVector->data->x = 0;
			cVector->data->y = 0;

		}
		if (cVector->data->x > 1 || cVector->data->x < -1) {
			cVector->data->x = 0;
			cVector->data->y = 0;
			CStage_Kinect::flag_valid_data = 0;

		}
		if (cVector->data->y > 1 || cVector->data->y < -1) {
			cVector->data->x = 0;
			cVector->data->y = 0;
			CStage_Kinect::flag_valid_data = 0;
		}
		if (cVector->data->x == 0 && cVector->data->y == 0) {
			//printf("buffered point is %f   and %f  \n", bufferedPoint.X, bufferedPoint.Y);
			//printf("calculated point is %f  and %f  \n", (-1)*(bufferedPoint.X - zero->X)/scale, (-1)*(bufferedPoint.Y - zero->Y)/scale);
		}
		//printf("vector values are %f and %f \n",cVector->data->x, cVector->data->y );
		return cVector;
	}

};

int CStage_Kinect::flag_hand = 0;
int CStage_Kinect::flag_valid_data = 1;
XnPoint3D CStage_Kinect::zero;
XnPoint3D CStage_Kinect::bufferedPoint;
xn::HandsGenerator CStage_Kinect::g_HandsGenerator;
Game* CStage_Kinect::game = Game::getInstance();

#endif /* CSTAGE_KINECT_HPP_ */
