#ifndef ALETAS_HPP
#define ALETAS_HPP

#include <memory>
#include "Eigen/Dense"
#include "FunDer.hpp"

typedef unsigned int uint;

class Aleta {
protected:
	// Parameters
	uint m_numberOfpoints;
	double m_theta_a, m_theta_b, m_k, m_h, m_D, m_L, m_delta;
	Eigen::MatrixXd m_A;
	Eigen::VectorXd m_b, m_x;
	std::shared_ptr<FunDer> m_Ac, m_As;

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

	// Define delta
	m_delta = m_L / (m_numberOfpoints - 1.0);

	// Initialize b vector
	m_b(0) = m_theta_a;
	m_b(m_numberOfpoints - 1) = m_theta_b;

	// Initialize A matrix
	m_A(0, 0) = m_A(m_numberOfpoints - 1, m_numberOfpoints - 1) = 1.0;
	for (uint i = 1; i < m_numberOfpoints - 1; i++) {
		m_A(i, i - 1) = this->coeff_i_minus_1(i);
		m_A(i, i) = this->coeff_i(i);
		m_A(i, i + 1) = this->coeff_i_plus_1(i);
	}
}

 void Aleta::solve() {
	 // Solves system of algebraic equations
	 m_x = m_A.colPivHouseholderQr().solve(m_b);
 }

 inline double Aleta::get_z(const uint& i) const
 {
	 //return m_L * static_cast<double>(i) / (m_numberOfpoints - 1.0);
	 return static_cast<double>(i) * m_delta;
 }

 inline double Aleta::coeff_i_minus_1(const uint& i) const
 {
	 //double C1 = m_Ac->der(this->get_z(i)) / m_Ac->fun(this->get_z(i));
	 //return 1.0 / (m_delta * m_delta) - C1 * .5 / m_delta;

	 double C1 = m_Ac->der(this->get_z(i)) / m_Ac->fun(this->get_z(i));
	 return 1.0 - m_delta * .5 * C1;
 }
 inline double Aleta::coeff_i(const uint& i) const
 {
	 //double C2 = (m_h/m_k)* m_As->der(this->get_z(i)) / m_Ac->fun(this->get_z(i));
	 //return -2.0 / (m_delta * m_delta) - C2;

	 double C2 = m_h* m_As->der(this->get_z(i))/( m_Ac->fun(this->get_z(i)) * m_k);
	 return -2.0 - m_delta * m_delta * C2;
 }
 inline double Aleta::coeff_i_plus_1(const uint& i) const
 {
	 //double C1 = m_Ac->der(this->get_z(i)) / m_Ac->fun(this->get_z(i));
	 //return 1.0 / (m_delta * m_delta) + C1 * .5 / m_delta;

	 double C1 = m_Ac->der(this->get_z(i)) / m_Ac->fun(this->get_z(i));
	 return 1.0 + m_delta * .5 * C1;
 }

 Eigen::MatrixXd Aleta::getA() const
{
	return m_A;
}

inline Eigen::VectorXd Aleta::getb() const
{
	return m_b;
}

inline Eigen::VectorXd Aleta::getx() const
{
	return m_x;
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