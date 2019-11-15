#ifndef FUNDER_HPP
#define FUNDER_HPP

class FunDer {
public:
	// Returns the function value
	virtual double fun(const double& z) const = 0;

	// Returns the derivative of function w.r.t. z value
	virtual double der(const double& z) const = 0;
};

#endif // FUNDER_HPP
