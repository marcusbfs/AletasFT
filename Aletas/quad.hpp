#ifndef QUAD_HPP
#define QUAD_HPP

#include <functional>

#pragma region DATA

// Gauss
#define c_p0 0.000000000000000
#define c_w0 0.417959183673469
#define c_p1 0.949107912342759
#define c_w1 0.129484966168870
#define c_p2 0.741531185599394
#define c_w2 0.279705391489277
#define c_p3 0.405845151377397
#define c_w3 0.381830050505119

// Kronrod
#define kc_p0 0.000000000000000
#define kc_w0 0.209482141084728
#define kc_p1 0.991455371120813
#define kc_w1 0.022935322010529
#define kc_p2 0.949107912342759
#define kc_w2 0.063092092629979
#define kc_p3 0.864864423359769
#define kc_w3 0.104790010322250
#define kc_p4 0.741531185599394
#define kc_w4 0.140653259715525
#define kc_p5 0.586087235467691
#define kc_w5 0.169004726639267
#define kc_p6 0.405845151377397
#define kc_w6 0.190350578064785
#define kc_p7 0.207784955007898
#define kc_w7 0.204432940075298

// Gauss arrays
double GAUSS_NODES[7] = { c_p0,
						c_p1, -c_p1,
						c_p2, -c_p2,
						c_p3, -c_p3 };
double GAUSS_WEIGHTS[7] = { c_w0,
						c_w1, c_w1,
						c_w2, c_w2,
						c_w3, c_w3 };
// Kronrod arrays
double KRONROD_NODES[15] = { kc_p0,
						kc_p1, -kc_p1,
						kc_p2, -kc_p2,
						kc_p3, -kc_p3,
						kc_p4, -kc_p4,
						kc_p5, -kc_p5,
						kc_p6, -kc_p6,
						kc_p7, -kc_p7 };
double KRONROD_WEIGHTS[15] = { kc_w0,
						kc_w1, kc_w1,
						kc_w2, kc_w2,
						kc_w3, kc_w3,
						kc_w4, kc_w4,
						kc_w5, kc_w5,
						kc_w6, kc_w6,
						kc_w7, kc_w7 };

#pragma endregion

double quad15points(double (*fun)(const double&), const double& ll, const double& hl) {
	double integral = 0.0;
	const double c1 = .5 * (hl - ll);
	const double c2 = .5 * (hl + ll);
	for (int i = 0; i < 15; i++) {
		integral += KRONROD_WEIGHTS[i] * fun(c1 * KRONROD_NODES[i] + c2);
	}
	return integral * c1;
}

double quad15points(std::function<double(const double&)>& fun, const double& ll, const double& hl) {
	double integral = 0.0;
	const double c1 = .5 * (hl - ll);
	const double c2 = .5 * (hl + ll);
	for (int i = 0; i < 15; i++) {
		integral += KRONROD_WEIGHTS[i] * fun(c1 * KRONROD_NODES[i] + c2);
	}
	return integral * c1;
}


#endif // QUAD_HPP
