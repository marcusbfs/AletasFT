#ifndef ALETATASK_HPP
#define ALETATASK_HPP

#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "FunDer.hpp"
#include "Aletas.hpp"
#include "quad.hpp"

#define ft_D 0.005
#define ft_L 0.1
#define ft_k 14.0
#define ft_h 5.0
#define ft_T0 150.0
#define ft_Tinf 20.0
#define PI EIGEN_PI
#define M_PI EIGEN_PI


// Class for calcultation of Ac and dAc/dz based on F(z)
class AletaTaskInput {
private:
	const double h = 1e-8;
public:
	// Returns F(z)
	virtual double F(const double& z) = 0;

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) = 0;

	// Returns integrand
	double m_int(const double& z) {
		return F(z) * std::sqrt(1.0 + std::pow(dFdz(z), 2));
	}

	// Returns As(z)
	virtual double As(const double& z) {
		double sum = 0;
		//for (int i = 0; i < 7; i++) {
		//	sum += GAUSS_WEIGHTS[i] * m_int(z * .5 * (1. + GAUSS_NODES[i]));
		//}
		for (int i = 0; i < 15; i++) {
			sum += KRONROD_WEIGHTS[i] * m_int(z * .5 * (1. + KRONROD_NODES[i]));
		}
		return sum * M_PI * z;
	}

	// Returns dAs(z)/dz
	virtual double dAsdz(const double& z) {
		return (As(h + z) - As(z - h)) / (2. * h);
	}

	// Return ID
	virtual std::string ID() = 0;
};

class AcWrapper : public FunDer {
private:
	std::shared_ptr<AletaTaskInput> m_Fz;
public:
	AcWrapper(const std::shared_ptr<AletaTaskInput>& f)
		: m_Fz(f) {}

	// Returns the function value
	double fun(const double& z) const {
		return EIGEN_PI * std::pow(m_Fz->F(z), 2);
	}

	// Returns the derivative of function w.r.t. z value
	double der(const double& z) const {
		return 2.0 * EIGEN_PI * m_Fz->F(z) * m_Fz->dFdz(z);
	}
};

class AsWrapper : public FunDer {
private:
	std::shared_ptr<AletaTaskInput> m_Fz;
public:
	AsWrapper(const std::shared_ptr<AletaTaskInput>& f)
		: m_Fz(f) {}

	// Returns the function value
	double fun(const double& z) const {
		// It's not needed
		return 0.0;
	}

	// Returns the derivative of function w.r.t. z value
	double der(const double& z) const {
		return m_Fz->dAsdz(z);
	}
};

class AletaTask  {
public:
	Aleta m_aleta; 
	double Tinf = ft_Tinf;
	double T0 = ft_T0;
	std::shared_ptr<AletaTaskInput> m_geratriz;
	Eigen::VectorXd m_Theta;
	Eigen::VectorXd m_T;

public:
	//AletaTask(const std::shared_ptr<FunDer>& Fz, const std::shared_ptr<FunDer>& As) {
	AletaTask(const uint numberOfPoints,
		const std::shared_ptr<AletaTaskInput>& geratriz)
		: m_geratriz(geratriz)
	{
		// Initialize functions
		m_aleta.setAc(std::make_shared<AcWrapper>(m_geratriz));
		m_aleta.setAs(std::make_shared<AsWrapper>(m_geratriz));

		// Parameters
		m_aleta.setD(ft_D);
		m_aleta.setL(ft_L);
		m_aleta.seth(ft_h);
		m_aleta.setk(ft_k);
		m_aleta.setNumberOfPoints(numberOfPoints);
		m_aleta.setBoundaryConditions(T0-Tinf, 0.0);
	}

	void build() {
		m_aleta.build();
	}

	void solve() {
		m_aleta.solve();
		m_T = getTheta();
		for (int i = 0; i < m_T.size(); i++)
			m_T(i) += Tinf;
	}

	std::string ID() const {
		return m_geratriz->ID();
	}

	Eigen::VectorXd getTheta() const {
		return m_aleta.getx();
	}

	Eigen::MatrixXd getA() const {
		return m_aleta.getA();
	}

	Eigen::VectorXd getT() const {
		return m_T;
	}

	Eigen::VectorXd getFlux() const {
		return this->m_aleta.getFlux();
	}

	Eigen::VectorXd getRate() const {
		return this->m_aleta.getRate();
	}

	void writeDataToFile(const std::string& filename) const {

		std::ofstream file;
		Eigen::VectorXd T = this->getT();
		Eigen::VectorXd q = this->getFlux();
		Eigen::VectorXd qd = this->getRate();

		file.open(filename);

		for (size_t i = 0; i < T.size(); i++) {
			file << std::fixed <<
				std::setprecision(10) << 
				std::setw(20) << m_aleta.get_z(i) << 
				std::setw(20) << T(i) <<
				std::setw(20) << q(i) <<
				std::setw(20) << qd(i) <<
				std::endl;
		}
		file.close();
	}
};

#endif // ALETATASK_HPP