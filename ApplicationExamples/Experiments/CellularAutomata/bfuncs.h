// The 8 b functinos as defined in 
// http://fab.cba.mit.edu/classes/862.16/notes/computation/Omohundro-1984.pdf

#include <cmath>

constexpr double mu = .01/2.0;

constexpr double SQ(const double z) {
	return z*z;
}

constexpr double b1(const double z) {
	return ((-.01 < z && z < .01) ?
		// using an alternative
		std::exp(-SQ(z)/SQ(mu)) : 0.0);
}

constexpr double b1int(const double z) {
	// Undetermined integral of b1
	return (z<-0.01 ? 0. : (
		z> 0.01 ? 1. : 
		0.5 * mu * std::sqrt(M_PI) * std::erf(z/mu)
	));
}

constexpr double b2(const double z) {
	return b1(z/20.)/b1(0.);
}

constexpr double b3(const double z) {
	return b1int(z)-b1int(-1.0); //)/(b1int(.01)-b1int(-.01)); // is 1
}

constexpr double b4(const double z) {
	return b3(z-.1)*b3(.1-z);
}

constexpr double b5(const double z) {
	return b3(z-.2)*b3(.2-z);
}

constexpr double around0(const double z) {
	// move z around the origin, so we don't have to evaluate the sum
	return (std::abs(z)-std::round(std::abs(z)));
}

constexpr double b6(const double z) {
	return b4(around0(z));
}

constexpr double b7(const double z) {
	return b5(around0(z));
}

constexpr double b8(const double z) {
	return b2(around0(z));
}

constexpr double db8(const double z) {
	// d b8(z) / dz
	return b1(around0(z));
}