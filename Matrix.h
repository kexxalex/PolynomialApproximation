#pragma once


template<class T> class Matrix {
public:
	Matrix(unsigned int r, unsigned int c, const std::vector<T>& values) : rows(r), columns(c) {
		entries = new T*[rows];
		for (unsigned _r = 0; _r < rows; _r++) {
			entries[_r] = new T[columns];
			for (unsigned _c = 0; _c < columns; _c++) {
				entries[_r][_c] = values[_r * columns + _c];
			}
		}
	}

	Matrix(unsigned int r, unsigned int c) : rows(r), columns(c) {
		entries = new T*[rows];
		for (unsigned _r = 0; _r < rows; _r++) {
			entries[_r] = new T[columns]{ 0 };
		}
	}

	Matrix(const Matrix<T>& m) : rows(m.rows), columns(m.columns) {
		entries = new T*[rows];
		for (unsigned r = 0; r < rows; r++) {
			entries[r] = new T[columns];
			for (unsigned c = 0; c < columns; c++)
				entries[r][c] = m.entries[r][c];
		}
	}

	void operator=(const Matrix<T>& m) {
		if (m.columns != columns || m.rows != rows)
			return;

		for (unsigned r = 0; r < rows; r++) {
			for (unsigned c = 0; c < columns; c++)
				entries[r][c] = m.entries[r][c];
		}
	}

	Matrix<T>& makeSymmetric(bool fromUpperTri = true) {
		for (unsigned i = 0; i < rows; i++) {
			for (unsigned j = i+1; j < columns; j++) {
				if (fromUpperTri)
					entries[j][i] = entries[i][j];
				else
					entries[i][j] = entries[j][i];
			}
		}
		return *this;
	}
	Matrix<T> transpose() const {
		Matrix m(columns, rows);
		for (unsigned r = 0; r < rows; r++) {
			for (unsigned c = 0; c < columns; c++) {
				m.entries[r][c] = entries[c][r];
			}
		}
		return m;
	}

	Matrix<T> invertSymmetric() const {
		Matrix<T> A(rows, columns * 2); // augmented matrix (this | 1)
		Matrix<T> L(rows, rows); // Lower triangular
		Matrix<T> DI(rows, rows); // Diagonal Inverse
		for (unsigned r = 0; r < rows; r++) {
			for (unsigned c = 0; c < columns; c++)
				A.entries[r][c] = entries[r][c];
			A.entries[r][r + rows] = 1;
		}

		A = A.reduce();
		for (unsigned r = 0; r < rows; r++) {
			for (unsigned c = 0; c < columns; c++)
				L.entries[r][c] = A.entries[r][c + rows];

			if (A.entries[r][r] == T(0))
				continue;
			DI.entries[r][r] = 1 / A.entries[r][r];
		}

		return L.transpose() * DI * L;
	}

	Matrix<T> operator*(const Matrix<T> A) const {
		Matrix<T> M(rows, A.columns);
		for (unsigned r = 0; r < rows; r++) {
			for (unsigned c = 0; c < A.columns; c++) {
				T val(0);
				for (unsigned i = 0; i < columns; i++) {
					val += entries[r][i] * A.entries[i][c];
				}
				M.entries[r][c] = val;
			}
		}
		return M;
	}

	Matrix<T>& reduce() {
		unsigned int col = 0;
		for (unsigned int r = 0; r < rows - 1; r++) {
			col = r;
			while (entries[r][col] == T(0) && col < columns)
				col++;
			if (col == columns)
				continue;

			for (unsigned int c_row = r + 1; c_row < rows; c_row++) {
				if (entries[c_row][col] == T(0))
					continue;

				T factor(-entries[c_row][col] / entries[r][col]);
				for (unsigned int c = col + 1; c < columns; c++)
					entries[c_row][c] += factor * entries[r][c];

				entries[c_row][col] = T(0);
			}
		}
		return *this;
	}

	T* const & operator[](unsigned index) { return entries[index]; }

	~Matrix() {
		for (unsigned r = 0; r < rows; r++)
			delete[] entries[r];
		delete[] entries;
	}

	const unsigned rows;
	const unsigned columns;

private:
	T** entries;
};