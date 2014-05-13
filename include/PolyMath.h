#ifndef POLYMATH_H
#define POLYMATH_H

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <vector>
#include <string>
//#include <numeric> // inner_product
//#include <sstream>
//#include "float.h"

namespace CoolProp{

/// The base class for Polynomials
class BasePolynomial{

protected:
	bool DEBUG;

public:
	// Constructor
	BasePolynomial();
	// Destructor.  No implementation
	virtual ~BasePolynomial(){};

protected:
	/// Basic checks for coefficient vectors.
	/** Starts with only the first coefficient dimension
	 *  and checks the vector length against parameter n. */
	bool checkCoefficients(const std::vector<long double> &coefficients, const unsigned int n);
	bool checkCoefficients(const std::vector< std::vector<long double> > &coefficients, const unsigned int rows, const unsigned int columns);

	/** Integrating coefficients for polynomials is done by dividing the
	 *  original coefficients by (i+1) and elevating the order by 1
	 *  through adding a zero as first coefficient.
	 *  Some reslicing needs to be applied to integrate along the x-axis.
	 *  In the brine/solution equations, reordering of the parameters
	 *  avoids this expensive operation. However, it is included for the
	 *  sake of completeness.
	 */
	std::vector<long double> integrateCoeffs(const std::vector<long double> &coefficients);
	std::vector< std::vector<long double> > integrateCoeffs(const std::vector< std::vector<long double> > &coefficients, bool axis);

	/** Deriving coefficients for polynomials is done by multiplying the
	 *  original coefficients with i and lowering the order by 1.
	 *
	 *  It is not really deprecated, but untested and therefore a warning
	 *  is issued. Please check this method before you use it.
	 */
	DEPRECATED(std::vector<long double> deriveCoeffs(const std::vector<long double> &coefficients));
	DEPRECATED(std::vector< std::vector<long double> > deriveCoeffs(const std::vector< std::vector<long double> > &coefficients, unsigned int axis));

private:
	/** The core of the polynomial wrappers are the different
	 *  implementations that follow below. In case there are
	 *  new calculation schemes available, please do not delete
	 *  the implementations, but mark them as deprecated.
	 *  The old functions are good for debugging since the
	 *  structure is easier to read than the backward Horner-scheme
	 *  or the recursive Horner-scheme.
	 */

	/// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficient vector.
	 *  Starts with only the first coefficient at T^0. */
	DEPRECATED(long double simplePolynomial(const std::vector<long double> &coefficients, long double T));
	DEPRECATED(long double simplePolynomial(const std::vector<std::vector<long double> > &coefficients, long double x, long double T));

	/// Simple integrated polynomial function generator.
	/** Base function to produce integrals of n-th order polynomials based on
	 *  the length of the coefficient vector.
	 *  Starts with only the first coefficient at T^0 */
	///Indefinite integral in T-direction
	long double simplePolynomialInt(const std::vector<long double> &coefficients, long double T);
	///Indefinite integral in T-direction only
	long double simplePolynomialInt(const std::vector<std::vector<long double> > &coefficients, long double x, long double T);

	/// Simple integrated polynomial function generator divided by independent variable.
	/** Base function to produce integrals of n-th order
	 *  polynomials based on the length of the coefficient
	 *  vector. Starts with only the first coefficient at T^0 */
	///Indefinite integral of a polynomial divided by its independent variable
	long double simpleFracInt(const std::vector<long double> &coefficients, long double T);
	///Indefinite integral of a polynomial divided by its 2nd independent variable
	long double simpleFracInt(const std::vector<std::vector<long double> > &coefficients, long double x, long double T);

	/** Simple integrated centred(!) polynomial function generator divided by independent variable.
	 *  We need to rewrite some of the functions in order to
	 *  use central fit. Having a central temperature Tbase
	 *  allows for a better fit, but requires a different
	 *  formulation of the fracInt function group. Other
	 *  functions are not affected.
	 *  Starts with only the first coefficient at T^0 */
	///Helper function to calculate the D vector:
	long double factorial(long double nValue);
	long double binom(long double nValue, long double nValue2);
	std::vector<long double> fracIntCentralDvector(int m, long double T, long double Tbase);
	///Indefinite integral of a centred polynomial divided by its independent variable
	long double fracIntCentral(const std::vector<long double> &coefficients, long double T, long double Tbase);

	/// Horner function generator implementations
	/** Represent polynomials according to Horner's scheme.
	 *  This avoids unnecessary multiplication and thus
	 *  speeds up calculation.
	 */
	long double baseHorner(const std::vector<long double> &coefficients, long double T);
	long double baseHorner(const std::vector< std::vector<long double> > &coefficients, long double x, long double T);
	///Indefinite integral in T-direction
	long double baseHornerInt(const std::vector<long double> &coefficients, long double T);
	///Indefinite integral in T-direction only
	long double baseHornerInt(const std::vector<std::vector<long double> > &coefficients, long double x, long double T);
	///Indefinite integral of a polynomial divided by its independent variable
	long double baseHornerFra(const std::vector<long double> &coefficients, long double T);
	///Indefinite integral of a polynomial divided by its 2nd independent variable
	long double baseHornerFra(const std::vector<std::vector<long double> > &coefficients, long double x, long double T);

	/** Alternatives
	 *  Simple functions that heavily rely on other parts of this file.
	 *  We still need to check which combinations yield the best
	 *  performance.
	 */
	///Derivative in T-direction
	long double deriveIn2Steps(const std::vector<long double> &coefficients, long double T);
	///Derivative in terms of x(axis=true) or T(axis=false).
	long double deriveIn2Steps(const std::vector< std::vector<long double> > &coefficients, long double x, long double T, bool axis);
	///Indefinite integral in T-direction
	long double integrateIn2Steps(const std::vector<long double> &coefficients, long double T);
	///Indefinite integral in terms of x(axis=true) or T(axis=false).
	long double integrateIn2Steps(const std::vector< std::vector<long double> > &coefficients, long double x, long double T, bool axis);
	///Indefinite integral in T-direction of a polynomial divided by its independent variable
	long double fracIntIn2Steps(const std::vector<long double> &coefficients, long double T);
	///Indefinite integral in T-direction of a polynomial divided by its 2nd independent variable
	long double fracIntIn2Steps(const std::vector<std::vector<long double> > &coefficients, long double x, long double T);
	///Indefinite integral of a centred polynomial divided by its 2nd independent variable
	long double fracIntCentral2Steps(const std::vector<std::vector<long double> > &coefficients, long double x, long double T, long double Tbase);

public:
	/** Here we define the functions that should be used by the
	 *  respective implementations. Please do no use any other
	 *  method since this would break the purpose of this interface.
	 *  Note that the functions below are supposed to be aliases
	 *  to implementations declared elsewhere in this file.
	 */

	/// Evaluates a one-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input
	inline long double polyval(const std::vector<long double> &coefficients, long double x){
		return baseHorner(coefficients,x);
	}

	/// Evaluates a two-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input in the 1st dimension
	/// @param y long double value that represents the current input in the 2nd dimension
	inline long double polyval(const std::vector< std::vector<long double> > &coefficients, long double x, long double y){
		//return simplePolynomial(coefficients,x,y);
		return baseHorner(coefficients,x,y);
	}

	/// Evaluates the indefinite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param T long double value that represents the current input
	inline long double polyint(const std::vector<long double> &coefficients, long double T){
		//return simplePolynomialInt(coefficients,T);
		return baseHornerInt(coefficients,T);
		//return integrateIn2Steps(coefficients,T);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial along the 2nd axis (T)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input in the 1st dimension
	/// @param T long double value that represents the current input in the 2nd dimension
	inline long double polyint(const std::vector< std::vector<long double> > &coefficients, long double x, long double T){
		//return simplePolynomialInt(coefficients,x,T);
		return baseHornerInt(coefficients,x,T);
		//return integrateIn2Steps(coefficients,x,T,false);
	}

	/// Evaluates the derivative of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param T long double value that represents the current input
	inline long double polyder(const std::vector<long double> &coefficients, long double T){
		return deriveIn2Steps(coefficients,T);
	}

	/// Evaluates the derivative of a two-dimensional polynomial along the 2nd axis (T)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input in the 1st dimension
	/// @param T long double value that represents the current input in the 2nd dimension
	inline long double polyder(const std::vector< std::vector<long double> > &coefficients, long double x, long double T){
		return deriveIn2Steps(coefficients,x,T,false);
	}

	/// Evaluates the indefinite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param T long double value that represents the current position
	inline long double polyfracint(const std::vector<long double> &coefficients, long double T){
		//return simpleFracInt(coefficients,T);
		return baseHornerFra(coefficients,T);
		//return fracIntIn2Steps(coefficients,T);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input in the 1st dimension
	/// @param T long double value that represents the current input in the 2nd dimension
	inline long double polyfracint(const std::vector< std::vector<long double> > &coefficients, long double x, long double T){
		//return simpleFracInt(coefficients,x,T);
		return baseHornerFra(coefficients,x,T);
		//return fracIntIn2Steps(coefficients,x,T);
	}

	/// Evaluates the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param T long double value that represents the current position
	/// @param Tbase central temperature for fitted function
	inline long double polyfracintcentral(const std::vector<long double> &coefficients, long double T, long double Tbase){
		return fracIntCentral(coefficients,T,Tbase);
	}

	/// Evaluates the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input in the 1st dimension
	/// @param T long double value that represents the current input in the 2nd dimension
	/// @param Tbase central temperature for fitted function
	inline long double polyfracintcentral(const std::vector< std::vector<long double> > &coefficients, long double x, long double T, long double Tbase){
		return fracIntCentral2Steps(coefficients,x,T,Tbase);
	}

	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param T long double value that represents the current input
	/// @param n int value that determines the kind of exponential function
	long double expval(const std::vector<long double> &coefficients, long double T, int n);

	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x long double value that represents the current input in the 1st dimension
	/// @param T long double value that represents the current input in the 2nd dimension
	/// @param n int value that determines the kind of exponential function
	long double expval(const std::vector< std::vector<long double> > &coefficients, long double x, long double T, int n);
};

/// The classes for Polynomials
class PolynomialImpl1D : public BasePolynomial{
protected:
	std::vector<long double> coefficients;
public:
	PolynomialImpl1D();
	PolynomialImpl1D(const std::vector<long double> &coefficients);
	virtual ~PolynomialImpl1D(){};
	long double eval(long double x);
	long double integ(long double x);
	long double deriv(long double x);
	long double solve(long double y, long double x0);
	//long double solveInt(long double y);
};

class PolynomialImpl2D : public BasePolynomial{
protected:
	std::vector< std::vector<long double> > coefficients;
public:
	PolynomialImpl2D();
	PolynomialImpl2D(const std::vector< std::vector<long double> > &coefficients);
	virtual ~PolynomialImpl2D(){};
	long double eval(long double x, long double z);
	long double solve(long double y, long double z);
};

class PolynomialInt1D : public BasePolynomial{
protected:
	std::vector<long double> coefficients;
public:
	PolynomialInt1D();
	PolynomialInt1D(const std::vector<long double> &coefficients);
	virtual ~PolynomialInt1D(){};
	long double eval(long double x);
	//long double integ(long double x);
	//long double deriv(long double x);
	long double solve(long double y);
	//long double solveInt(long double y);
};

class PolynomialInt2D : public BasePolynomial{
protected:
	std::vector< std::vector<long double> > coefficients;
public:
	PolynomialInt2D();
	PolynomialInt2D(const std::vector< std::vector<long double> > &coefficients);
	virtual ~PolynomialInt2D(){};
	long double eval(long double x, long double z);
	long double solve(long double y, long double z);
};

class PolynomialFracInt1D : public BasePolynomial{
protected:
	std::vector<long double> coefficients;
public:
	PolynomialFracInt1D();
	PolynomialFracInt1D(const std::vector<long double> &coefficients);
	virtual ~PolynomialFracInt1D(){};
	long double eval(long double x);
	long double solve(long double y);
};

class PolynomialFracInt2D : public BasePolynomial{
protected:
	std::vector< std::vector<long double> > coefficients;
public:
	PolynomialFracInt2D();
	PolynomialFracInt2D(const std::vector< std::vector<long double> > &coefficients);
	virtual ~PolynomialFracInt2D(){};
	long double eval(long double x, long double z);
	long double solve(long double y, long double z);
};



///// The classes for Polynomials
//class PolynomialImpl1D : public BasePolynomial{
//protected:
//	std::vector<long double> coefficients;
//public:
//	PolynomialImpl1D(const std::vector<long double> &coefficients);
//	virtual ~PolynomialImpl1D(){};
//	long double eval(long double x);
//	long double integ(long double x);
//	long double deriv(long double x);
//	long double solve(long double y);
//	long double solveInt(long double y);
//};
//
//class PolynomialImpl2D : public BasePolynomial{
//private:
//	std::vector< std::vector<long double> > coefficients;
//public:
//	PolynomialImpl2D(const std::vector< std::vector<long double> > &coefficients);
//	virtual ~PolynomialImpl2D(){};
//	long double eval(long double x, long double z);
//	long double integ(long double x, long double z);
//	long double deriv(long double x, long double z);
//	long double solve(long double y, long double z);
//	long double solveInt(long double y, long double z);
//};


}; /* namespace CoolProp */
#endif
