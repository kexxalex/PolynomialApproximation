#pragma once
#include <vector>
#include <iostream>
#include <string>

inline unsigned int max(unsigned int a, unsigned int b) {
	return (a > b) ? a : b;
}


template<class T> class Polynomial
{
public:
	unsigned int degree{ 0 };

	Polynomial(const std::vector<T>& coeff = { 0 }) : degree(coeff.size() - 1), coefficients(coeff.begin(), coeff.end()) { shrink(); }
	Polynomial(const Polynomial<T>& poly) : Polynomial<T>(poly.coefficients) { shrink(); }
	Polynomial(unsigned int deg) : degree(deg), coefficients(deg + 1, T(0)) {}

	static Polynomial<T> basis(unsigned int deg) {
		Polynomial p(deg);
		p.get(deg) = T(1);
		return p;
	}

	static Polynomial<T> ones(unsigned int deg) {
		Polynomial p(std::vector<T>(deg+1, T(1)));
		return p;
	}

	Polynomial<T>& operator=(const Polynomial<T>& other) {
		coefficients.assign(other.coefficients.begin(), other.coefficients.end());
		return shrink();
	}

	~Polynomial() { coefficients.clear(); }

	void print(const std::string& fname = "p", const std::string& variable = "x", bool newline = true) const {
		std::cout << fname << "(" << variable << ") = ";
		for (unsigned int i = degree; i >= 2; i--) {
			const T& val = (*this)[i];
			if (val == T(0))
				continue;
			if (i == degree && val < 0)
				std::cout << "-";
			else if (i < degree && val < 0)
				std::cout << " - ";
			else if (i < degree && val >= 0)
				std::cout << " + ";
			std::cout << ((val < 0) ? -val : val) << variable << "^" << i;
		}

		T val = (*this)[1];
		if (val != T(0)) {
			if (degree > 1) {
				std::cout << ((val < 0) ? " - " : " + ") << ((val < 0) ? -val : val) << variable;
			}
			else if (degree == 1) {
				std::cout << ((val < 0) ? "-" : "") << ((val < 0) ? -val : val) << variable;
			}
		}

		val = (*this)[0];
		if (degree == 0) {
			if (val == T(0))
				std::cout << 0;
			else
				std::cout << ((val < 0) ? "-" : "") << ((val < 0) ? -val : val);
		}
		else {
			if (val != T(0))
				std::cout << ((val < 0) ? " - " : " + ") << ((val < 0) ? -val : val);
		}
		if (newline)
			std::cout << std::endl;
	}

	Polynomial<T> differentiate() const {
		Polynomial p(degree);

		for (unsigned int i = 0; i < degree; i++)
			p.get(i) = (i + 1) * (*this)[i + 1];

		return p.shrink();
	}

	Polynomial<T>& operator+=(const Polynomial<T>& other) {
		unsigned int deg = max(degree, other.degree);

		for (unsigned int i = 0; i <= deg; i++)
			get(i) += other[i];

		return *this;
	}

	Polynomial<T>& operator-=(const Polynomial<T>& other) {
		return (*this) += -other;
	}

	Polynomial<T>& operator*=(const Polynomial<T>& other) {
		*this = (*this) * other;
		return *this;
	}

	Polynomial<T>& operator*=(const T& scalar) {
		for (unsigned int i = 0; i <= degree; i++)
			get(i) *= scalar;
		return shrink();
	}

	Polynomial<T> operator+(const Polynomial<T>& other) const {
		Polynomial p(*this);
		p += other;
		return p;
	}

	Polynomial<T> operator-(const Polynomial<T>& other) const {
		Polynomial p(*this);
		p -= other;
		return p;
	}

	Polynomial<T> operator+(const T& constant) {
		Polynomial p(*this);
		p.get(0) += constant;
		return p;
	}

	Polynomial<T> operator-(const T& constant) {
		Polynomial p(*this);
		p.get(0) -= constant;
		return p;
	}

	Polynomial<T> operator-() const {
		Polynomial p(*this);
		for (unsigned int i = 0; i <= degree; i++)
			p.get(i) *= T(-1);

		return p;
	}

	Polynomial<T> operator*(const Polynomial<T>& other) const {
		std::vector<T> new_coeff(degree + other.degree + 1, T(0));

		for (unsigned int i = 0; i <= degree; i++) {
			for (unsigned int j = 0; j <= other.degree; j++) {
				new_coeff[i + j] += (*this)[i] * other[j];
			}
		}
		return Polynomial(new_coeff);
	}

	Polynomial<T> operator*(T scalar) const {
		Polynomial p(*this);
		return p *= scalar;
	}

	T operator[](unsigned int index) const {
		if (index > degree)
			return T(0);

		return coefficients[index];
	}

	T operator()(const T& value) const {
		T result(0);
		for (unsigned int i = 0; i <= degree; i++) {
			T product((*this)[i]);
			for (unsigned int j = 0; j < i; j++)
				product *= value;

			result += product;
		}
		return result;
	}

	T& get(unsigned int index) {
		if (index > degree) {
			coefficients.insert(coefficients.end(), index - degree, T(0));
			degree = index;
		}

		return coefficients[index];
	}

	Polynomial<T>& shrink() {
		if (degree == coefficients.size() - 1 && get(degree) != T(0))
			return *this;


		unsigned int deg = coefficients.size() - 1;
		unsigned int i = deg;
		for (; i <= deg; i--) {
			if (get(i) != T(0))
				break;
		}

		if (i == 0xFFFFFFFF) {
			coefficients.clear();
			coefficients.emplace_back(T(0));
		}

		else if (i < deg) {
			coefficients.erase(coefficients.begin() + i + 1, coefficients.end());
		}
		degree = coefficients.size() - 1;

		return *this;
	}

private:
	std::vector<T> coefficients; // Basis (1, x, x^2, x^3, ...)
};
