#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "hsv.h"

/**
 * A smallprogram to read out colorvalues of pixels on the screen.
 * Left-Click: get color values.
 * Right-Click: exit.
 * 
 * @autor: Marcel Schneider; mostly based on Stackoverflow snippets [1] and [2]
 * [1]: http://stackoverflow.com/questions/16122196/getting-mouseclick-coordinates-with-xlib
 * [2]: http://stackoverflow.com/questions/17518610/how-to-get-a-screen-pixels-color-in-x11
 */
int main() {
	int x = -1, y = -1;
	XEvent event;
	int button = 0;
	Display *display = XOpenDisplay(NULL);
	if (display == NULL) {
		fprintf(stderr, "Cannot connect to X server!\n");
		exit(EXIT_FAILURE);
	}
	Window root = XDefaultRootWindow(display);
	XGrabPointer(display, root, False, ButtonReleaseMask, GrabModeAsync,
	GrabModeAsync, None, None, CurrentTime);

	XSelectInput(display, root, ButtonReleaseMask);
	printf("Right-click to exit.\n");
	while (button != Button3) {
		while (1) {
			XNextEvent(display, &event);
			switch (event.type) {
			case ButtonRelease:
				switch (event.xbutton.button) {
				case Button1:
					x = event.xbutton.x;
					y = event.xbutton.y;
					button = Button1;
					break;

				case Button3:
					x = event.xbutton.x;
					y = event.xbutton.y;
					button = Button3;
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
			if (x >= 0 && y >= 0)
				break;
		}

		XColor c;
		XImage *image;
		image = XGetImage(display, RootWindow(display, DefaultScreen (display)),
				x, y, 1, 1, AllPlanes, XYPixmap);
		c.pixel = XGetPixel(image, 0, 0);
		XFree(image);
		XQueryColor(display, DefaultColormap(display, DefaultScreen (display)),
				&c);
		c.red /= 256;
		c.green /= 256;
		c.blue /= 256;
		RgbColor rgb;
		rgb.r = c.red;
		rgb.g = c.green;
		rgb.b = c.blue;
		HsvColor hsv = RgbToHsv(rgb);
		printf("x,y: %d,%d\n\trgb(%3d,%3d,%3d)\n", x, y, c.red, c.green,
				c.blue);
		printf("\thsv(%3d,%3d,%3d)\n", hsv.h, hsv.s, hsv.v);
	}
	XCloseDisplay(display);
	return 0;
}
