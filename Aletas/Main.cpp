#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include "Eigen/Dense"
#include "Aletas.hpp"
#include "AletaTask.hpp"
#include "FunDer.hpp"

void test_AletaTask();

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
	double b = (ft_D - a) / ft_L;

	// Returns F(z)
	virtual double F(const double& z) {
		return  a + b * z;
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return b;
	}

	// Returns dAs(z)/dz
	virtual double dAsdz(const double& z) {
		return 2.0 * PI * std::sqrt(1.0 + std::pow(b, 2)) * (a + b * z);
	}

	// Return ID
	virtual std::string ID() {
		return "B";
	}
};

// F(z) = a - b*z
class GeratrizC : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = a / ft_L;

	// Returns F(z)
	virtual double F(const double& z) {
		return  a - b * z;
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return -b;
	}

	// Returns dAs(z)/dz
	virtual double dAsdz(const double& z) {
		return 2.0 * PI * std::sqrt(1.0 + std::pow(b, 2)) * (a - b * z);
	}

	// Return ID
	virtual std::string ID() {
		return "C";
	}
};

// F(z) = a + b*z*z
class GeratrizD : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (ft_D -a)/(ft_L*ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return a + b * std::pow(z, 2);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return 2 * b * z;
	}

	// Returns dAs(z)/dz
	virtual double dAsdz(const double& z) {
		return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) - 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) + 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) - 3.0 / 2.0 * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) + (9.0 / 8.0) * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) - 1.0 / 8.0 * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	}

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
	virtual double dAsdz(const double& z) {
		return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	}

	// Return ID
	virtual std::string ID() {
		return "E";
	}
};

// F(z) = a + b*z*z*z
class GeratrizF : public AletaTaskInput {
	// Parâmetros!! (se houver)
	double a = ft_D / 2.0;
	double b = (ft_D - a)/(ft_L*ft_L*ft_L);

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
	double b = a/(ft_L*ft_L*ft_L);

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
	double b = (ft_D - a) / std::sin(ft_L);

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
	double b = a / std::sin(ft_L);

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
	//double a = (1.0 / 2.0) * ft_D * (std::exp(2 * ft_L) - 4 * std::exp(ft_L) + 1) / (std::exp(2 * ft_L) - 2 * std::exp(ft_L) + 1);
	double a = ft_D * (.5 - 1. / std::cosh(ft_L)) / (1. - 1. / std::cosh(ft_L));
	double b = (-a + ft_D)/ std::cosh(ft_L);

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
	double a = ft_D * .5 / (1. - 1. / std::cosh(ft_L));
	double b = a / std::cosh(ft_L);

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
	double a = ft_D * (.5 - 1. / std::exp(ft_L)) / (1. - 1. / std::exp(ft_L));
	double b = (ft_D - a) / std::exp(ft_L);

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
	double a = ft_D * .5 / (1. - 1. / std::exp(ft_L));
	double b = a / std::exp(ft_L);

	// Returns F(z)
	virtual double F(const double& z) {
		return a - b * std::exp(z);
	}

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) {
		return  - b * std::exp(z);
	}

	// Returns dAs(z)/dz
	//virtual double dAsdz(const double& z) {
	//	return 2 * M_PI * (2 * a * std::pow(b, 2) * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a * std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 2.0) * a / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + 4 * std::pow(b, 5) * std::pow(z, 6) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 5 * std::pow(b, 3) * std::pow(z, 4) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (3.0 / 2.0) * std::pow(b, 3) * std::pow(z, 4) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0) - 9.0 / 8.0 * b * std::pow(z, 2) / std::sqrt(4 * std::pow(b, 2) * std::pow(z, 2) + 1) + (1.0 / 8.0) * b * std::pow(z, 2) / std::pow(4 * std::pow(b, 2) * std::pow(z, 2) + 1, 3.0 / 2.0));
	//}

	// Return ID
	virtual std::string ID() {
		return "M";
	}
};

// ==== não altere mais! ======= !

int main()
{
	test_AletaTask();
	std::getchar();
}

void test_AletaTask()
{
	const int numberOfPoints = 9;

	std::unique_ptr<AletaTask> aleta;
	std::vector<std::shared_ptr<AletaTaskInput>> geratrizes;
	std::string filename;


	geratrizes.clear();

	// ========== Altere aqui! ==========
#pragma region  Geratrizes
	std::shared_ptr<AletaTaskInput> geratrizA = std::make_shared<GeratrizA>();
	geratrizes.push_back(geratrizA);

	std::shared_ptr<AletaTaskInput> geratrizB = std::make_shared<GeratrizB>();
	geratrizes.push_back(geratrizB);

	std::shared_ptr<AletaTaskInput> geratrizC = std::make_shared<GeratrizC>();
	geratrizes.push_back(geratrizC);

	std::shared_ptr<AletaTaskInput> geratrizD = std::make_shared<GeratrizD>();
	geratrizes.push_back(geratrizD);

	std::shared_ptr<AletaTaskInput> geratrizE = std::make_shared<GeratrizE>();
	geratrizes.push_back(geratrizE);

	std::shared_ptr<AletaTaskInput> geratrizF = std::make_shared<GeratrizF>();
	geratrizes.push_back(geratrizF);

	std::shared_ptr<AletaTaskInput> geratrizG = std::make_shared<GeratrizG>();
	geratrizes.push_back(geratrizG);

	std::shared_ptr<AletaTaskInput> geratrizH = std::make_shared<GeratrizH>();
	geratrizes.push_back(geratrizH);

	std::shared_ptr<AletaTaskInput> geratrizI = std::make_shared<GeratrizI>();
	geratrizes.push_back(geratrizI);

	std::shared_ptr<AletaTaskInput> geratrizJ = std::make_shared<GeratrizJ>();
	geratrizes.push_back(geratrizJ);

	std::shared_ptr<AletaTaskInput> geratrizK = std::make_shared<GeratrizK>();
	geratrizes.push_back(geratrizK);

	std::shared_ptr<AletaTaskInput> geratrizL = std::make_shared<GeratrizL>();
	geratrizes.push_back(geratrizL);

	std::shared_ptr<AletaTaskInput> geratrizM = std::make_shared<GeratrizM>();
	geratrizes.push_back(geratrizM);
#pragma endregion

	// ==== não altere mais! ======= !

	// Para cada geratriz
	for (int i = 0; i < geratrizes.size(); i++) {
		aleta = std::make_unique<AletaTask>(numberOfPoints, geratrizes[i]);

		std::cout << "Geratriz " << aleta->ID() << std::endl;
		Eigen::VectorXd T = aleta->getT();

		std::cout << T << std::endl;
		//std::cout << aleta.getA() << std::endl;

		std::cout << "\n " <<   std::endl;

		filename = "geratriz_" + aleta->ID() + ".txt";
		aleta->writeTtoFile(filename);
		aleta.release();
	}
}
