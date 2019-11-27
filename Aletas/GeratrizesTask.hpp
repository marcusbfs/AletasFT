#ifndef GERATRIZES_TASK_HPP
#define GERATRIZES_TASK_HPP

#include "Aletas.hpp"
#include "AletaTask.hpp"
#include "FunDer.hpp"
#include "GeratrizesTask.hpp"

// ========== Altere aqui! ==========

// F(z) = D/2
class GeratrizA : public AletaTaskInput {
	// Returns F(z)
	virtual double F(const double& z) {
		return  .5 * ft_D;
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return 0.0;
	}

	// Returns dAs(z)/dz
	virtual double dAsdz(const double& z) {
		// As = PI * D * z
		return PI * ft_D;
	}

	// Return ID
	virtual std::string ID() {
		return "A";
	}
};

// F(z) = a + b*z
class GeratrizB : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = ft_D / (2 * ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return  a + b * z;
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return b;
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2.0 * PI * std::sqrt(1.0 + std::pow(b, 2)) * (a + b * z);
	//}

	// Return ID
	virtual std::string ID() {
		return "B";
	}
};

// F(z) = a - b*z
class GeratrizC : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = ft_D / (2 * ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return  a - b * z;
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return -b;
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2.0 * PI * std::sqrt(1.0 + std::pow(b, 2)) * (a - b * z);
	//}

	// Return ID
	virtual std::string ID() {
		return "C";
	}
};

// F(z) = a + b*z*z
class GeratrizD : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (1.0 / 2.0) * ft_D / std::pow(ft_L, 2);

	// Returns F(z)
	virtual double F(const double& z) {
		return a + b * std::pow(z, 2);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return 2 * b * z;
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) - 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) + 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) - 3.0 / 2.0 * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) + (9.0 / 8.0) * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) - 1.0 / 8.0 * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "D";
	}
};

// F(z) = a - b*z*z
class GeratrizE : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = a/(ft_L*ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return a - b * std::pow(z, 2);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return -2 * b * z;
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "E";
	}
};

// F(z) = a + b*z*z*z
class GeratrizF : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (1.0 / 2.0) * ft_D / std::pow(ft_L, 3);

	// Returns F(z)
	virtual double F(const double& z) {
		return a + b * std::pow(z, 3);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return 3 * b * std::pow(z, 2);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "F";
	}
};

// F(z) = a - b*z*z*z
class GeratrizG : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (1.0 / 2.0) * ft_D / std::pow(ft_L, 3);

	// Returns F(z)
	virtual double F(const double& z) {
		return a - b * std::pow(z, 3);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return -3 * b * std::pow(z, 2);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "G";
	}
};

// F(z) = a + b*sin(z)
class GeratrizH : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (1.0 / 2.0) * ft_D / std::sin(ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return a + b * std::sin(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return b * std::cos(z);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "H";
	}
};

// F(z) = a - b*sin(z)
class GeratrizI : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (1.0 / 2.0) * ft_D / std::sin(ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return a - b * std::sin(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return -b * std::cos(z);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "I";
	}
};

// F(z) = a + b*cosh(z)
class GeratrizJ : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = (1.0 / 2.0) * ft_D * (std::cosh(ft_L) - 2) / (std::cosh(ft_L) - 1);
	double b = (1.0 / 2.0) * ft_D / (std::cosh(ft_L) - 1);

	// Returns F(z)
	virtual double F(const double& z) {
		return a + b * std::cosh(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return b * std::sinh(z);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "J";
	}
};

// F(z) = a - b*cosh(z)
class GeratrizK : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = (1.0 / 2.0) * ft_D * std::cosh(ft_L) / (std::cosh(ft_L) - 1);
	double b = (1.0 / 2.0)* ft_D / (std::cosh(ft_L) - 1);

	// Returns F(z)
	virtual double F(const double& z) {
		return a - b * std::cosh(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return  - b * std::sinh(z);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "K";
	}
};

// F(z) = a + b*exp(z)
class GeratrizL : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = (1.0 / 2.0) * ft_D * (std::exp(ft_L) - 2) / (std::exp(ft_L) - 1);
	double b = (1.0 / 2.0) * ft_D / (std::exp(ft_L) - 1);

	// Returns F(z)
	virtual double F(const double& z) {
		return a + b * std::exp(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return   b * std::exp(z);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "L";
	}
};

// F(z) = a - b*exp(z)
class GeratrizM : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = (1.0 / 2.0) * ft_D * std::exp(ft_L) / (std::exp(ft_L) - 1);
	double b = (1.0 / 2.0) * ft_D / (std::exp(ft_L) - 1);

	// Returns F(z)
	virtual double F(const double& z) {
		return a - b * std::exp(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return  - b * std::exp(z);
	}

	// Return ID
	virtual std::string ID() {
		return "M";
	}
};

// ==== não altere mais! ======= !
#endif // GERATRIZES_TASK_HPP
