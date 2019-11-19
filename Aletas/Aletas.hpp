#ifndef ALETAS_HPP
#define ALETAS_HPP

#include <memory>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "FunDer.hpp"

typedef unsigned int uint;

class Aleta {
protected:
	// Parameters
	uint m_numberOfpoints;
	double m_theta_a, m_theta_b, m_k, m_h, m_D, m_L, m_delta;
	//Eigen::MatrixXd m_A;
	Eigen::SparseMatrix<double> m_A;
	Eigen::VectorXd m_b, m_x;
	std::shared_ptr<FunDer> m_Ac, m_As, m_vol;

	// Private functions
	inline double coeff_i_minus_1(const uint& i) const;
	inline double coeff_i(const uint& i) const;
	inline double coeff_i_plus_1(const uint& i) const;

public:
	// Return A matrix (Ax = b)
	Eigen::MatrixXd getA() const;
	// Return b vector (Ax = b)
	Eigen::VectorXd getb() const;
	// Return x vector (Ax = b)
	Eigen::VectorXd getx() const;
	// Return flux q = -kdT/dz
	Eigen::VectorXd getFlux() const;
	// Return reate q_dot = -k Ac dT/dz
	Eigen::VectorXd getRate() const;
	// Return the volume of aleta ayy
	double getVolume() const;

	// Get z position
	inline double get_z(const uint& i) const;

	// Set number of discretization points
	void setNumberOfPoints(const uint& n);

	// Set boundary conditions
	void setBoundaryConditions(const double& theta_a, const double& theta_b);

	// Set Ac
	void setAc(const std::shared_ptr<FunDer>& Ac);
	// Set As
	void setAs(const std::shared_ptr<FunDer>& As);
	// Set Volume
	void setVol(const std::shared_ptr<FunDer>& vol);
	// Set k
	void setk(const double& k);
	// Set h
	void seth(const double& h);
	// Set L
	void setL(const double& L);
	// Set D
	void setD(const double& D);

	// Build problem
	void build();

	// Solve problem
	void solve();

};

 void Aleta::build()
{
	// Resize arrays and zero'em
	m_A.resize(m_numberOfpoints, m_numberOfpoints);
	m_b.resize(m_numberOfpoints);
	m_x.resize(m_numberOfpoints);

	m_b.setZero();

	m_A.setZero();
	m_A.reserve(2 + 3 * (m_numberOfpoints - 2));

	// Define delta
	m_delta = m_L / (m_numberOfpoints - 1.0);

	// Boundary conditions

	// theta(0) = theta_a
	m_b(0) = m_theta_a;
	m_A.insert(0, 0) = 1.0;

	constexpr bool known_theta_b = false;

	if (known_theta_b) {
		// theta(L) = theta_b
		m_b(m_numberOfpoints - 1) = m_theta_b;
		m_A.insert(m_numberOfpoints - 1, m_numberOfpoints - 1) = 1.0; // theta(L) = theta_b
	}
	else {
		// h*theta(L) = - k * dtheta(L)/dz
		m_b(m_numberOfpoints - 1) = 0.0;

		constexpr bool version_1 = true;
		if (version_1) {
			// Version 1
			m_A.insert(m_numberOfpoints - 1, m_numberOfpoints - 2) = -m_k;
			m_A.insert(m_numberOfpoints - 1, m_numberOfpoints - 1) = m_k + m_delta * m_h;
		}
		else {
			// Version 2
			m_A.insert(m_numberOfpoints - 1, m_numberOfpoints - 3) = m_k * .5;
			m_A.insert(m_numberOfpoints - 1, m_numberOfpoints - 2) = -2.0 * m_k;
			m_A.insert(m_numberOfpoints - 1, m_numberOfpoints - 1) = -m_h*m_delta + 3. * m_k * .5;
		}
	}

	// Fill sparse matrix
	for (uint i = 1; i < m_numberOfpoints - 1; i++) {
		m_A.insert(i, i - 1) = this->coeff_i_minus_1(i);
		m_A.insert(i, i) = this->coeff_i(i);
		m_A.insert(i, i + 1) = this->coeff_i_plus_1(i);
	}
	//m_A.makeCompressed();

	// Fill dense matrix
	//for (uint i = 1; i < m_numberOfpoints - 1; i++) {
	//	m_A(i, i - 1) = this->coeff_i_minus_1(i);
	//	m_A(i, i) = this->coeff_i(i);
	//	m_A(i, i + 1) = this->coeff_i_plus_1(i);
	//}
}

 void Aleta::solve() {
	 // Solves system of algebraic equations
	 //m_x = m_A.fullPivHouseholderQr().solve(m_b);
	 //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	 Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
	 solver.compute(m_A);
	 m_x = solver.solve(m_b);
 }

 inline double Aleta::get_z(const uint& i) const
 {
	 //return m_L * static_cast<double>(i) / (m_numberOfpoints - 1.0);
	 return static_cast<double>(i) * m_delta;
 }

 inline double Aleta::coeff_i_minus_1(const uint& i) const
 {
	 double C1 = m_Ac->der(this->get_z(i));
	 return 1.0*m_Ac->fun(this->get_z(i)) - m_delta * .5 * C1;
 }
 inline double Aleta::coeff_i(const uint& i) const
 {
	 double C2 = m_h* m_As->der(this->get_z(i))/(  m_k);
	 return -2.0 * m_Ac->fun(this->get_z(i)) - m_delta * m_delta * C2;
 }
 inline double Aleta::coeff_i_plus_1(const uint& i) const
 {
	 double C1 = m_Ac->der(this->get_z(i));
	 return 1.0* m_Ac->fun(this->get_z(i)) + m_delta * .5 * C1;
 }

 Eigen::MatrixXd Aleta::getA() const
{
	//return m_A;
	return Eigen::MatrixXd(m_A);
}

inline Eigen::VectorXd Aleta::getb() const
{
	return m_b;
}

inline Eigen::VectorXd Aleta::getx() const
{
	return m_x;
}

inline Eigen::VectorXd Aleta::getFlux() const
{
	Eigen::VectorXd q(m_numberOfpoints);

	q(0) = -m_k*(m_x(1) - m_x(0))/ m_delta;
	q(m_numberOfpoints-1) = -m_k*(m_x(m_numberOfpoints-1) - m_x(m_numberOfpoints-2))/ m_delta;

	for (uint i = 1; i < m_numberOfpoints - 1; i++) {
		q(i) = - m_k* (m_x(i + 1) - m_x(i - 1)) * .5 / m_delta;
	}

	return q;
}

inline Eigen::VectorXd Aleta::getRate() const
{
	Eigen::VectorXd qd = this->getFlux();

	for (uint i = 0; i < m_numberOfpoints; i++) {
		qd(i) *= this->m_Ac->fun(this->get_z(i));
	}

	return qd;
}

inline double Aleta::getVolume() const
{
	return m_vol->fun(m_L);
}


inline void Aleta::setNumberOfPoints(const uint& n)
{
	m_numberOfpoints = n;
}

inline void Aleta::setBoundaryConditions(const double& theta_a, const double& theta_b)
{
	m_theta_a = theta_a;
	m_theta_b = theta_b;
}

inline void Aleta::setAc(const std::shared_ptr<FunDer>& Ac)
{
	m_Ac = Ac;
}

inline void Aleta::setAs(const std::shared_ptr<FunDer>& As)
{
	m_As = As;
}

inline void Aleta::setVol(const std::shared_ptr<FunDer>& vol)
{
	m_vol = vol;
}

inline void Aleta::seth(const double& h)
{
	m_h = h;
}

inline void Aleta::setL(const double& L)
{
	m_L = L;
}

inline void Aleta::setD(const double& D)
{
	m_D = D;
}

inline void Aleta::setk(const double& k)
{
	m_k = k;
}

#endif // ALETAS_HPP