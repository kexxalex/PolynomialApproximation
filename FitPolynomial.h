#pragma once
#include "Polynomial.h"
#include "Matrix.h"

template<class T> struct point {
	T x, y;
};

template <class T> inline T pow(T val, unsigned int exp) {
	T p(1);
	for (unsigned int k = 0; k < exp; k++)
		p *= val;

	return p;
}

template<class T> class PolyFitter
{
public:
	PolyFitter(unsigned int degree) : fit_poly(degree) {
		for (unsigned d = 0; d < degree; d++)
			coefficients.push_back(d);
	}
	PolyFitter(const std::vector<unsigned int>& coefficients) : fit_poly(), coefficients(coefficients.begin(), coefficients.end()) {}

	const Polynomial<T>& do_fit(const std::vector<point<T>>& measurements) {
		unsigned cc = coefficients.size();

		Matrix<T> M(cc, cc);
		Matrix<T> Y(cc, 1);

		for (unsigned k = 0; k < cc; k++) {
			for (unsigned j = k; j < cc; j++) {
				T& val = M[k][j];
				for (const auto& p : measurements)
					val += pow<T>(p.x, coefficients[k] + coefficients[j]);
			}

			T& val = Y[k][0];
			for (const auto& p : measurements)
				val += p.y * pow<T>(p.x, coefficients[k]);
		}
		Matrix<T> I = M.makeSymmetric().invertSymmetric();
		Matrix<T> solution = I * Y;

		for (unsigned i = 0; i < cc; i++)
			fit_poly.get(coefficients[i]) = *solution[i];

		variance = T(0);
		for (const auto& p : measurements) {
			T val(fit_poly(p.x) - p.y);
			variance += val*val;
		}
		variance /= measurements.size() - 1;

		return fit_poly.shrink();
	}

	const Polynomial<T>& getPolynomial() const { return fit_poly; }
	const T& getVariance() const { return variance; }

private:
	Polynomial<T> fit_poly;
	std::vector<unsigned int> coefficients;

	std::vector<T> approx_coeff;
	T variance{ 0 };
};

typedef PolyFitter<float> PolyFitterF;
typedef PolyFitter<double> PolyFitterD;
typedef PolyFitter<long double> PolyFitterLD;
