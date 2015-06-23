/**
 * Basic image operations that need to be fast. 
 * 
 * Mainly anisotropic diffusion and blur code. The code is not highly 
 * optimized, but designed with speed in mind.
 * 
 * All operations expect one channel (grayscale) pixelmaps as line-major arrays.
 * 
 * Author: Marcel Schneider <marcel.schneider@studi.informatik.uni-stuttgart.de> 
 * 
 */

/**
 * Representation of a pixel. Integer types do not work for anisotropic 
 * diffusion, float causes rounding errors for the integral image 
 * (8 + log(n_pixels) bits of precision required).
 * Performance difference is hardly measurable. 
 * 
 * Range of pixel values is always 0..255.
 */
#ifdef USE_FFT
typedef double pixel_t;
#else
typedef int pixel_t;
#endif

static uint8_t to_uint8(pixel_t x) {
	if (x > 255)
		x = 255;
	if (x < 0)
		x = 0;
	return (uint8_t) x;
}

/*
 * Anisotropic Diffusion
 */

static pixel_t phi(pixel_t x) {
	float d = (float) x;
	float p = d / 10.0;
	//return x;
	//return (pixel_t) (exp(-(p*p)) * d);
	float r = (d / (1.0 + p * p));
	//r += r > 0.0 ? 0.9 : -0.9;
	return (pixel_t) r;
}

static pixel_t phi_tbl[256];

static void make_phi_tbl() {
	for (int i = 0; i < 256; i++) {
		phi_tbl[i] = phi(i);
	}
}

static pixel_t phi_tbl_lookup(pixel_t x) {
	int i = (x > 0 ? x : -x);
	return (x > 0 ? 1 : -1) * phi_tbl[i];
}

void do_diff(pixel_t *src, pixel_t *dest, int w, int h) {
	for (int y = 1; y < h - 1; y++) {
		for (int x = 1; x < w - 1; x++) {
			pixel_t current_i = src[y * w + x];
			pixel_t ci[4];
			ci[0] = src[(y - 1) * w + x];
			ci[1] = src[y * w + (x + 1)];
			ci[2] = src[(y + 1) * w + x];
			ci[3] = src[y * w + (x - 1)];

			pixel_t delta = 0;
			for (int r = 0; r < 4; r++) {
				pixel_t nabla = ci[r] - current_i;
				delta += phi_tbl_lookup(nabla);
				//delta += phi(nabla);
			}
			dest[y * w + x] = current_i + (delta / 4);
		}
	}
}

/*
 * Integral Image
 */

void integrate_image(pixel_t *__restrict src, pixel_t *__restrict dest, int w,
		int h) {
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			pixel_t left = x > 0 ? dest[y * w + (x - 1)] : 0;
			pixel_t top = y > 0 ? dest[(y - 1) * w + x] : 0;
			pixel_t topleft = x > 0 && y > 0 ? dest[(y - 1) * w + (x - 1)] : 0;
			pixel_t i = src[y * w + x];
			dest[y * w + x] = left + top - topleft + i;
		}
	}
}

pixel_t integral(int from_x, int from_y, int to_x, int to_y,
		pixel_t *integral_image, int w) {
	pixel_t topleft =
			(from_y > 0 && from_x > 0) ?
					integral_image[(from_y - 1) * w + (from_x - 1)] : 0;
	pixel_t topright =
			(from_y > 0 && to_x > 0) ?
					integral_image[(from_y - 1) * w + (to_x - 1)] : 0;
	pixel_t botleft =
			(to_y > 0 && from_x > 0) ?
					integral_image[(to_y - 1) * w + (from_x - 1)] : 0;
	pixel_t botright =
			(to_y > 0 && to_x > 0) ?
					integral_image[(to_y - 1) * w + (to_x - 1)] : 0;
	return botright - topright - botleft + topleft;
}

void rect_blur(pixel_t *integral_image, pixel_t *dest, int s, int w, int h) {
	float inv_s = 1.0 / ((s + 1 + s) * (s + 1 + s));
	for (int y = s; y < h - s; y++) {
		for (int x = s; x < w - s; x++) {
			pixel_t sum = integral(x - s, y - s, x + s + 1, y + s + 1,
					integral_image, w);
			dest[y * w + x] = sum * inv_s;
		}
	}
}

void rect_blur_unnormalized(pixel_t *integral_image, pixel_t *dest, int s,
		int w, int h) {
	for (int y = s; y < h - s; y++) {
		for (int x = s; x < w - s; x++) {
			pixel_t sum = integral(x - s, y - s, x + s + 1, y + s + 1,
					integral_image, w);
			dest[y * w + x] = sum;
		}
	}
}

void corner_detect(pixel_t *integral_image, pixel_t *dest, int s, int w, int h) {
	int inv_s = (s + 1 + s) * (s + 1 + s);
	for (int y = s; y < h - 3 * s; y++) {
		for (int x = s; x < w - 3 * s; x++) {
			pixel_t sum1 = integral(x - s, y - s, x + s + 1, y + s + 1,
					integral_image, w);
			pixel_t sum2 = integral(x + s, y - s, x + 3 * s + 1, y + s + 1,
					integral_image, w);
			pixel_t sum3 = integral(x - s, y + s, x + s + 1, y + 3 * s + 1,
					integral_image, w);
			pixel_t sum4 = integral(x + s, y + s, x + 3 * s + 1, y + 3 * s + 1,
					integral_image, w);
			dest[y * w + x] = (sum1 + sum2 + sum3 - 3 * sum4) / 6 / inv_s;
		}
	}
}

// get rid of pixels missing or sticking out. This is neccessary for Navier Stokes equations.
int cleanup_pixels(uint8_t *img, int w, int h) {
	uint8_t *top = img - w;
	uint8_t *left = img - 1;
	uint8_t *right = img + 1;
	uint8_t *bot = img + w;
	int flip_count = 0;
	for (int y = 1; y < h - 1; y++) {
		for (int x = 1; x < w - 1; x++) {
			int idx = y * w + x;
			int count_non_white = (top[idx] != 0) + (left[idx] != 0)
					+ (right[idx] != 0) + (bot[idx] != 0);
			if (img[idx] != 0 && count_non_white == 2 && top[idx] == bot[idx]
					&& left[idx] == right[idx]) {
				img[idx] = 0;
				flip_count++;
			} else if (img[idx] == 0 && count_non_white >= 3) {
				flip_count++;
				img[idx] = 1;
			} else if (img[idx] != 0 && count_non_white <= 1) {
				flip_count++;
				img[idx] = 0;
			}
		}
	}
	return flip_count;
}

int image_delta(uint8_t *a, uint8_t *b, int len, uint8_t epsilon, int step = 2) {
	int diff = 0;
	for (int i = 0; i < len; i += step) {
		if ((a[i] - b[i]) > epsilon || (b[i] - a[i]) > epsilon) {
			diff++;
		}
	}
	return diff;
}

