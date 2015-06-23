// if difference of two pixel values is higher then CHANGE_EPSILON, pixels are marked as different
#define CHANGE_EPSILON		10

// define change threshold so that a new frame is committed to the pipeline. The global variable output_cDataArray2D_uc1.has_changed is relying on this threshold.
#define CHANGE_THRESHOLD	0.005

// change_threshold_max defines the percentage of pixels that have to be detected in order to freeze screen for hand detection. Threshold is responsible for variable transmit_output.
#define CHANGE_THRESHOLD_MAX	0.01

// padding around every frame to avoid boundary problems with Gaussian filter
#define PADDING_WIDTH 		20

// size of fine Gaussian Blur and coarse (rough) gaussian blur
#define DOG_FINE 		2
#define DOG_COARSE 		10

// initial threshold for object detection
#define INITIAL_THRESHOLD 	120

// initial hysteresis threshold for object detection
#define INITIAL_HYSTERESIS	5

// define 'hue' and 'value' and 'saturation' of hsv color to detect RED color
#define COLOR_INPUT(c)		(c .h < 20 || c .h > 250) && c .v > 120
#define INPUT_MIN_SATURATION 	100
// define 'hue' and 'saturation' of hsv color to detect BLUE color
#define COLOR_OUTPUT(c)		c .h > 140 && c .h < 170 && c. v > 110
#define OUTPUT_MIN_SATURATION 	150
// define 'hue' and 'value' and 'saturation' of hsv color to detect GREEN color
#define COLOR_TARGET(c)		c .h > 90 && c .h < 140 && c .v > 90
#define TARGET_MIN_SATURATION 	170

// passes, size, redThreshold, blueThreshold, greenThreshold
#define POSTPROCESS_PARAMETERS 	5, 2, 7, 3, 3 
