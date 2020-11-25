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
	const double m_h = 1e-8;
	std::function<double(const double&)> m_AsIntegrand =
		[this](const double& _z) {
		return this->F(_z) * std::sqrt(1.0 + std::pow(this->dFdz(_z), 2));
	};
	std::function<double(const double&)> m_VolumeIntegrand =
		[this](const double& _z) {
		return std::pow(this->F(_z), 2);
	};
public:
	// Returns F(z)
	virtual double F(const double& z) = 0;

	// Returns dF(z)/dz
	virtual double dFdz(const double& z) = 0;

	// Returns As(z)
	virtual double As(const double& z) {
		return quad15points(this->m_AsIntegrand, 0, static_cast<double>(z)) * 2.0 * M_PI;
	}

	// Returns dAs(z)/dz
	virtual double dAsdz(const double& z) {
		//return (As(m_h + z) - As(z - m_h)) / (2. * m_h);
		return 2.0 * M_PI * this->m_AsIntegrand(z);
	}

	// Returns PI* integral(ll, hl, F(z)^2)
	virtual double Volume(const double& ll, const double& hl) {
		return quad15points(this->m_VolumeIntegrand, ll, hl) * M_PI;
	}

	// Return ID
	virtual std::string ID() = 0;

	// destructor
	virtual ~AletaTaskInput() {};
};

class AcWrapper : public FunDer {
private:
	AletaTaskInput* m_Fz;
public:
	AcWrapper(AletaTaskInput* f)
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
	AletaTaskInput* m_Fz;
public:
	AsWrapper(AletaTaskInput* f)
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

class VolumeWrapper : public FunDer {
private:
	AletaTaskInput* m_Fz;
public:
	VolumeWrapper(AletaTaskInput* f)
		: m_Fz(f) {}

	// Returns the function value
	double fun(const double& z) const {
		// It's not needed
		return m_Fz->Volume(0.0, z);
	}

	// Returns the derivative of function w.r.t. z value
	double der(const double& z) const {
		return 0.0;
	}
};

class AletaTask  {
public:
	Aleta m_aleta; 
	double Tinf = ft_Tinf;
	double T0 = ft_T0;
	AletaTaskInput* m_geratriz;
	Eigen::VectorXd m_Theta;
	Eigen::VectorXd m_T;
private:
	std::unique_ptr<AcWrapper> m_Ac;
	std::unique_ptr<AsWrapper> m_As;
	std::unique_ptr<VolumeWrapper> m_Volume;

public:
	//AletaTask(const std::shared_ptr<FunDer>& Fz, const std::shared_ptr<FunDer>& As) {
	AletaTask(const uint numberOfPoints,
		AletaTaskInput* geratriz)
		: m_geratriz(geratriz)
	{
		// Initialize functions
		m_Ac = std::make_unique<AcWrapper>(m_geratriz);
		m_As = std::make_unique<AsWrapper>(m_geratriz);
		m_Volume = std::make_unique<VolumeWrapper>(m_geratriz);

		m_aleta.setAc(m_Ac.get());
		m_aleta.setAs(m_As.get());
		m_aleta.setVol(m_Volume.get());

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
		const double volume = this->m_aleta.getVolume();

		file.open(filename);

		file << 
			std::setw(20) << "Position" <<
			std::setw(20) << "Temperature" <<
			std::setw(20) << "Flux" <<
			std::setw(20) << "Flux derivative" <<
			std::setw(20) << "Volume" <<
			std::endl;
		for (size_t i = 0; i < T.size(); i++) {
			file << std::fixed <<
				std::setprecision(10) << 
				std::setw(20) << m_aleta.get_z(i) << 
				std::setw(20) << T(i) <<
				std::setw(20) << q(i) <<
				std::setw(20) << qd(i) <<
				std::setw(20) << volume <<
				std::endl;
		}
		file.close();
	}
};

#endif // ALETATASK_HPP