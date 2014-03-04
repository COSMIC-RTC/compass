#include <carma_utils.h>

void getNumBlocksAndThreads(int device, int n, int &blocks, int &threads) {

	struct cudaDeviceProp deviceProperties;
	cudaGetDeviceProperties(&deviceProperties, device);

	int maxThreads = deviceProperties.maxThreadsPerBlock;
	blocks = deviceProperties.multiProcessorCount * 8;
	threads = (n + blocks - 1) / blocks;

	if (threads > maxThreads) {
		threads = maxThreads;
		blocks =
				MIN( deviceProperties.maxGridSize[0],(n + threads - 1) / threads);
	}
}

float ran1() {
	float norm;
	norm = 2147483647.f;

	//srand(time(NULL)); // initialisation de rand
	//return rand()/norm;
	return random() / norm;
}

extern "C" {
int _dist(float *d, long dimx, long dimy, float xc, float yc) {
	/* Declarations */
	long i, j;

	/* Loop and fill d with distance values */
	for (i = 0; i < dimx; ++i) {
		for (j = 0; j < dimy; ++j) {
			d[i + j * dimx] = (float) sqrt(
					(xc - (float) i) * (xc - (float) i)
							+ (yc - (float) j) * (yc - (float) j));
		}
	}
	return EXIT_SUCCESS;
}

void _poidev(float *xmv, long n)
/* */
{
	float gammln(float xx);
	/*  float ran1(long *idum);*/
	static float sq, alxm, g, oldm = (-1.0);
	float xm, em, t, y;
	long i;

	for (i = 0; i < n; i++) {
		xm = (float) xmv[i];
		if (xm == 0.0f)
			continue;
		if (xm < 20.0) { /* Use direct method. */
			if (xm != oldm) {
				oldm = xm;
				g = exp(-xm); /* If xm is new, compute the exponential. */
			}
			em = -1;
			t = 1.0;
			do {
				++em;
				t *= ran1();
			} while (t > g);
		} else { /* Use rejection method. */
			if (xm != oldm) {
				oldm = xm;
				sq = sqrt(2.0 * xm);
				alxm = log(xm);
				g = xm * alxm - gammln(xm + 1.0);
			}
			do {
				do {
					y = tan(3.1415926535897932384626433832 * ran1());
					em = sq * y + xm;
				} while (em < 0.0);
				em = floor(em);
				t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - g);
			} while (ran1() > t);
		}
		xmv[i] = (float) em;
	}
}

float gammln(float xx) {
	/* Returns the value ln[?(xx)] for xx>0. */
	float x, y, tmp, ser;
	static float cof[6] = { 76.18009172947146, -86.50532032941677,
			24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
			-0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++)
		ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);
}

}
