#include <iostream>
#include "FitPolynomial.h"
#include "Matrix.h"

using namespace std;

typedef long double APPROX;

int main(int argc, char** args) {
	if (argc <= 2)
		return 1;
	
	unsigned int coeff_count = stoi(args[1]);
	unsigned int data_length = stoi(args[2]);

	if (coeff_count == 0)
		return 2;

	if (data_length == 0 || 3 + coeff_count > argc)
		return 3;

	std::vector<unsigned int> coefficients;
	for (unsigned i = 0; i < coeff_count; i++)
		coefficients.emplace_back(stoi(args[3+i]));

	PolyFitter<APPROX> fitter(coefficients);

	std::vector<point<APPROX>> measurements;
	for (unsigned i = 0; i < data_length; i++) {
		point<APPROX> p;
		cin >> p.x;
		cin >> p.y;
		measurements.emplace_back(p);
	}

	auto p = fitter.do_fit(measurements);

	for (unsigned int i = 0; i < p.degree; i++)
		cout << p[i] << ";";
	cout << p[p.degree] << ";" << fitter.getVariance();

	coefficients.clear();
	measurements.clear();

	return 0;
}
