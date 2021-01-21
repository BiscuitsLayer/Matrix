#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <exception>

#define DEBUG(var) std::cout << "DEBUG: " << #var << " = " << var << std::endl;

std::ostream& operator << (std::ostream& stream, const std::vector <double>& vec) {
	for (size_t i = 0; i < vec.size (); ++i) {
		stream << vec [i] << ' ';
	}
	return stream;
}

//	Генератор тестов (+ добавить ответов для E2E тестов)
//	Exception-safety в операторах копирования и присваивания

namespace Linear {

	const double EPS = 1e-3;
	using pair_size_t = std::pair <size_t, size_t>;

	template <typename T>
	class Matrix;

}

namespace Determinant {

	enum class Type {
		ERROR = 0,
		FULL = 1,
		GAUSS = 2
	};

	template <typename T>
	T Full (const Linear::Matrix <T>& matrix);
    double Gauss (const Linear::Matrix <double>& matrix);
}

namespace Linear {

	template <typename T>
	class Matrix final {
		private:
			size_t nRows_ = 0, nCols_ = 0;
			T** data_ = nullptr;
			int gaussFactor_ = 1;

		public:
			//	CTORS AND DTORS
						Matrix 	(size_t rows, size_t cols, T value = T{});
			explicit 	Matrix 	(size_t rows = 0);
						~Matrix ();

			//	FROM VECTOR
			Matrix (size_t rows, size_t cols, std::vector <T> &vec);

			//	FROM ANOTHER MATRIX
			Matrix (const Matrix& rhs);
			Matrix (Matrix&& rhs) noexcept;

			//	OVERLOADED OPERATORS AND METHODS
			Matrix& operator = 	(const Matrix& rhs);
			Matrix& operator = 	(Matrix&& rhs);
			bool 	operator == (const Matrix& rhs) const;
			bool 	operator != (const Matrix& rhs) const;
			Matrix& operator *= (const T number) 	&;
			Matrix& operator *= (const Matrix& rhs) &;
			Matrix& operator += (const Matrix& rhs) &;
			Matrix& operator -= (const Matrix& rhs) &;
			operator std::vector <T>  ();
			Matrix& Transpose 	() &;
			Matrix& Negate 		() &;
			Matrix& Clear 		() &;
			Matrix& Diagonalize () &;

			//	ROW AND COLUMN OPERATIONS
			Matrix& SwapRows 	(size_t lhs, size_t rhs);
			Matrix& AddRows		(size_t source, size_t destination, T factor);
			Matrix& SwapCols 	(size_t lhs, size_t rhs);
			Matrix& AddCols 	(size_t source, size_t destination, T factor);

			//	BASIC TYPES
			static Matrix Zeros (size_t n);
			static Matrix Eye 	(size_t n);

			//	GETTERS
			pair_size_t Shape		() const;
			size_t 		Size 		() const;
			T 			Trace 		() const;
			const T& 	At 			(size_t i, size_t j) const;
			void 		Dump 		(std::ostream& stream) const;
			int			GaussFactor () const;

			//	SETTERS
			T& At (size_t i, size_t j);

			//	ALGEBRA
			T Determinant (Determinant::Type type = Determinant::Type::ERROR) const;
			size_t Rank () const;
	};

	//	INPUT AND OUTPUT
	template <typename T>
	std::istream& operator >> (std::istream& stream, Matrix <T>& rhs);
	template <typename T>
	std::ostream& operator << (std::ostream& stream, const Matrix <T>& rhs);

	//	OVERLOADED OPERATORS
	template <typename T>
	Matrix <T> operator + (const Matrix <T>& lhs, const Matrix <T>& rhs);
	template <typename T>
	Matrix <T> operator - (const Matrix <T>& lhs, const Matrix <T>& rhs);
	template <typename T>
	Matrix <T> operator * (const T number, const Matrix <T>& rhs);
	template <typename T>
	Matrix <T> operator * (const Matrix <T>& lhs, const T number);
	template <typename T>
	Matrix <T> operator * (const Matrix <T>& lhs, const Matrix <T>& rhs);
}


// MAIN CODE //


template <typename T>
T Determinant::Full (const Linear::Matrix <T>& matrix) {
	auto shape = matrix.Shape ();
	size_t nRows = shape.first, nCols = shape.second;
	if (nRows != nCols) {
		throw (std::invalid_argument ("Trying to calcute non-square matrix determinant."));
	}
	if (nRows == 1) {
		return matrix.At (0, 0);
	}
	else {
		T ans {};
		for (size_t i = 0; i < nCols; ++i) {
			Linear::Matrix <T> temp { nRows - 1, nCols - 1 };
			//	Creating minor
			for (size_t j = 0; j < nRows - 1; ++j) {
				for (size_t k = 0; k < nCols; ++k) {
					if (k < i) {
						temp.At (j, k) = matrix.At (j + 1, k);
					}
					else if (k == i) {
						continue;
					}
					else {
						temp.At (j, k - 1) = matrix.At (j + 1, k);
					}
				}
			}
			//	Main sum
			ans += ( i % 2 == 0 ? + 1 : - 1 ) * matrix.At (0, i) * Determinant::Full (temp);
		}
		ans = (ans < Linear::EPS ? 0 : ans);
		return std::round (ans);
	}
}

double Determinant::Gauss (const Linear::Matrix <double>& matrix) {
	auto shape = matrix.Shape ();
	size_t nRows = shape.first, nCols = shape.second;
	if (nRows != nCols) {
		throw (std::invalid_argument ("Trying to calcute non-square matrix determinant."));
	}
	if (nRows == 1) {
		return matrix.At (0, 0);
	}
	else {
		double ans = 1.0;
		Linear::Matrix <double> temp = matrix;
		temp.Diagonalize ();
		ans *= temp.GaussFactor ();

		for (size_t i = 0; i < nRows; ++i) {
			ans *= temp.At (i, i);
		}
		ans = (ans < Linear::EPS ? 0 : ans);
		return std::round (ans);
	}
}

template <typename T>
Linear::Matrix <T>::Matrix (size_t rows, size_t cols, T value):
	nRows_ (rows),
	nCols_ (cols),
	data_ (new T* [nRows_] {}),
	gaussFactor_ (1)
	{
		if (nRows_ * nCols_ == 0) {
			nRows_ = nCols_ = 0;
		}
		else {
			for (size_t i = 0; i < nRows_; ++i) {
				data_[i] = new T [nCols_] {};
			}
		}
	}

template <typename T>
Linear::Matrix <T>::Matrix (size_t rows):
	Matrix (rows, rows)
	{}

template <typename T>
Linear::Matrix <T>::~Matrix () {
	for (size_t i = 0; i < nRows_; ++i) {
		delete [] data_[i];
	}
	delete [] data_;
}

template <typename T>
Linear::Matrix <T>::Matrix (size_t rows, size_t cols, std::vector <T> &vec):
	Matrix (rows, cols)
	{
		if (vec.size () != nRows_ * nCols_) {
			throw (std::invalid_argument ("Vector and Matrix sizes do not match."));
		}
		for (size_t pos = 0; pos < std::min (nRows_ * nCols_, vec.size ()); ++pos) {
			data_[pos / nCols_][pos % nCols_] = vec [pos];
		}
	}

template <typename T>
Linear::Matrix <T>::Matrix (const Matrix& rhs) {
	nRows_ = rhs.nRows_;
	nCols_ = rhs.nCols_;
	data_ = new T* [nRows_] {};
	for (size_t i = 0; i < nRows_; ++i) {
		data_[i] = new T [nCols_] {};
		for (size_t j = 0; j < nCols_; ++j) {
			data_[i][j] = rhs.data_[i][j];
		}
	}
	gaussFactor_ = rhs.gaussFactor_;
	//	No return value in ctor
}

template <typename T>
Linear::Matrix <T>::Matrix (Matrix&& rhs) noexcept {
	std::swap (nRows_, rhs.nRows_);	
	std::swap (nCols_, rhs.nCols_);
	std::swap (data_, rhs.data_);
	std::swap (gaussFactor_, rhs.gaussFactor_);
	//	No return value in ctor
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::operator = (const Matrix& rhs) {
	if (this != &rhs) {
		delete [] data_;
		nRows_ = rhs.nRows_;
		nCols_ = rhs.nCols_;
		data_ = new T* [nRows_] {};
		for (size_t i = 0; i < nRows_; ++i) {
			data_[i] = new T [nCols_] {};
			for (size_t j = 0; j < nCols_; ++j) {
				data_[i][j] = rhs.data_[i][j];
			}
		}
		gaussFactor_ = rhs.gaussFactor_;
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::operator = (Matrix&& rhs) {
	if (this != &rhs) {
		std::swap (nRows_, rhs.nRows_);	
		std::swap (nCols_, rhs.nCols_);
		std::swap (data_, rhs.data_);
		std::swap (gaussFactor_, rhs.gaussFactor_);
	}
	return *this;
}

template <typename T>
bool Linear::Matrix <T>::operator == (const Matrix& rhs) const {
	if (Shape () != rhs.Shape()) {
		return false;
	}
	for (size_t i = 0; i < nRows_; ++i) {
		for (size_t j = 0; j < nCols_; ++j) {
			if (data_[i][j] != rhs.data_[i][j]) {
				return false;
			}
		}
	}
	return true;
}

template <typename T>
bool Linear::Matrix <T>::operator != (const Matrix& rhs) const {
	return !(*this == rhs);
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::operator *= (const T number) & { 
	//	MULTIPLY BY THE NUMBER (OF THE SAME TYPE)
	for (size_t i = 0; i < nRows_; ++i) {
		for (size_t j = 0; j < nCols_; ++j) {
			data_[i][j] *= number;
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::operator *= (const Matrix& rhs) & {	
	// MULTIPLY BY ANOTHER MATRIX (OF THE SAME TYPE)
	auto lhsShape = Shape (), rhsShape = rhs.Shape ();
	if (lhsShape.second != rhsShape.first) {
		throw (std::invalid_argument ("Trying to multiply by a matrix with inappropriate size."));
	}
	else {
		Matrix <T> ans { lhsShape.first, rhsShape.second };
		for (size_t i = 0; i < lhsShape.first; ++i) {
			for (size_t j = 0; j < rhsShape.second; ++j) {
				for (size_t k = 0; k < lhsShape.second; ++k) {
					ans.data_[i][j] += data_[i][k] * rhs.data_[k][j];
				}
			}
		}
		*this = std::move (ans);
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::operator += (const Matrix& rhs) & {
	if (this->Shape () != rhs.Shape ()) {
		throw (std::invalid_argument ("Matrix sizes do not match."));
	}
	else {
		for (size_t i = 0; i < nRows_; ++i) {
			for (size_t j = 0; j < nCols_; ++j) {
				data_[i][j] += rhs.data_[i][j];
			}
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::operator -= (const Matrix& rhs) & {
	if (this->Shape () != rhs.Shape ()) {
		throw (std::invalid_argument ("Matrix sizes do not match."));
	}
	else {
		for (size_t i = 0; i < nRows_; ++i) {
			for (size_t j = 0; j < nCols_; ++j) {
				data_[i][j] -= rhs.data_[i][j];
			}
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>::operator std::vector <T> () {
	std::vector <T> ans {};
	for (size_t i = 0; i < nRows_; ++i) {
		for (size_t j = 0; j < nCols_; ++j) {
			ans.push_back (At (i, j));
		}
	}
	return ans;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::Transpose () & {
	Matrix temp {nCols_, nRows_};
	for (size_t i = 0; i < nRows_; ++i) {
		for (size_t j = 0; j < nCols_; ++j) {
			temp.data_[j][i] = data_[i][j];
		}
	}
	*this = std::move (temp);
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::Negate () & {
	for (size_t i = 0; i < nRows_; ++i) {
		for (size_t j = 0; j < nCols_; ++j) {
			data_[i][j] = (- 1) * data_[i][j];
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::Clear () & {
	*this = std::move (Matrix <T> {});
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::Diagonalize () & {
	for (size_t i = 0; i < std::min (nRows_, nCols_); ++i) {
		T maxElement = {};
		size_t maxIdx = i + 1;
		for (size_t j = i + 1; j < nRows_; ++j) {
			if (At (j, i) > maxElement) {
				maxElement = At (j, i);
				maxIdx = j;
			}
		}
		if (maxElement == T {}) {
			//	Matrix has a zero-column
			continue;
		}
		else {
			SwapRows (i, maxIdx);
		}
		for (size_t j = i + 1; j < nRows_; ++j) {
			AddRows (i, j, (- 1) * (At (j, i) / At (i, i)));
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::SwapRows (size_t lhs, size_t rhs) {
	if (lhs >= nRows_ || rhs >= nRows_ || lhs == rhs) {
		throw (std::invalid_argument ("Wrong lhs / rhs value."));
	}
	else {
		std::swap (data_[lhs], data_[rhs]);
		gaussFactor_ *= -1;
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::AddRows (size_t source, size_t destination, T factor) {
	if (source >= nRows_ || destination >= nRows_ || source == destination) {
		throw (std::invalid_argument ("Wrong lhs / rhs value."));
	}
	else {
		for (size_t i = 0; i < nCols_; ++i) {
			data_[destination][i] += data_[source][i] * factor;
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::SwapCols (size_t lhs, size_t rhs) {
	if (lhs >= nCols_ || rhs >= nCols_ || lhs == rhs) {
		throw (std::invalid_argument ("Wrong lhs / rhs value."));
	}
	else {
		for (size_t i = 0;  i < nRows_; ++i) {
			std::swap (data_[i][lhs], data_[i][rhs]);
			gaussFactor_ *= -1;
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix <T>::AddCols (size_t source, size_t destination, T factor) {
	if (source >= nCols_ || destination >= nCols_ || source == destination) {
		throw (std::invalid_argument ("Wrong source / destination value."));
	}
	else {
		for (size_t i = 0; i < nRows_; ++i) {
			data_[i][destination] += data_[i][source] * factor;
		}
	}
	return *this;
}

template <typename T>
Linear::Matrix <T> Linear::Matrix <T>::Zeros (size_t n) {
	Matrix <T> temp {n};
	//	No std::move here because of RVO
	return temp;
}

template <typename T>
Linear::Matrix <T> Linear::Matrix <T>::Eye (size_t n) {
	Matrix <T> temp {n};
	for (size_t i = 0; i < n; ++i) {
		temp.data_[i][i] = static_cast <T> (1);
	}
	//	No std::move here because of RVO
	return temp;
}

template <typename T>
Linear::pair_size_t Linear::Matrix <T>::Shape () const {
	return std::make_pair (nRows_, nCols_);
}

template <typename T>
size_t Linear::Matrix <T>::Size () const {
	return nRows_ * nCols_;
}

template <typename T>
T Linear::Matrix <T>::Trace () const {
	T ans {};
	for (size_t i = 0; i < std::min (nRows_, nCols_); ++i) {
		ans += data_[i][i];
	}
	return ans;
}

template <typename T>
const T& Linear::Matrix <T>::At (size_t i, size_t j) const {
	if (i >= nRows_ || j >= nCols_) {
		throw (std::invalid_argument ("Wrong i / j value."));
		return data_[0][0];
	}
	else {
		return data_[i][j];
	}
}

template <typename T>
void Linear::Matrix <T>::Dump (std::ostream& stream) const {
	stream << "Matrix Dump: rows = " << nRows_ << ", columns = " << nCols_;
	if (nRows_ * nCols_ != 0) {
		stream << std::endl;
	}
	for (size_t i = 0; i < nRows_; ++i) {
		for (size_t j = 0; j < nCols_; ++j) {
			stream << (j == 0 ? "" : " ") << data_[i][j];
		}
		if (i != nRows_ - 1)
			stream << std::endl;
	}
}

template <typename T>
int Linear::Matrix <T>::GaussFactor () const {
	return gaussFactor_;
}

template <typename T>
T& Linear::Matrix <T>::At (size_t i, size_t j) {
	return const_cast <T&> (static_cast <const Matrix<T>*> (this)->At (i, j));
}

template <typename T>
T Linear::Matrix <T>::Determinant (Determinant::Type type) const {
	switch (type) {
		case Determinant::Type::ERROR: {
			throw (std::runtime_error ("Invalid determinant type."));
			return 0;
		}
		case Determinant::Type::FULL: {
			return Determinant::Full (*this);
		}
		case Determinant::Type::GAUSS: {
			return Determinant::Gauss (*this);
		}
	}
	return 0;
}

template <typename T>
size_t Linear::Matrix <T>::Rank () const {
	auto shape = Shape ();
	Linear::Matrix <T> temp = *this;
	temp.Diagonalize ();
	DEBUG (temp);

	size_t ans = 0;
	for (size_t i = 0; i < temp.nRows_; ++i) {
		if (std::fabs (temp.At (i, i)) > EPS) {
			++ans;
			continue;
		}
	}
	return ans;
}

template <typename T>
std::istream& Linear::operator >> (std::istream& stream, Matrix <T>& rhs) {
	size_t rows = 0, cols = 0;
	stream >> rows >> cols;

	std::vector <T> temp (rows * cols);
	for (T& element : temp) {
		stream >> element;
	}
	rhs = {rows, cols, temp};
	return stream;
}

template <typename T>
std::ostream& Linear::operator << (std::ostream& stream, const Matrix <T>& rhs) {
	rhs.Dump (stream);
	return stream;
}

template <typename T>
Linear::Matrix <T> Linear::operator + (const Matrix <T>& lhs, const Matrix <T>& rhs) {
	Matrix <T> temp { lhs };
	temp += rhs;
	//	No std::move here because of RVO
	return temp;
}

template <typename T>
Linear::Matrix <T> Linear::operator - (const Matrix <T>& lhs, const Matrix <T>& rhs) {
	Matrix <T> temp { lhs };
	temp -= rhs;
	//	No std::move here because of RVO
	return temp;
}

template <typename T>
Linear::Matrix <T> Linear::operator * (const T number, const Matrix <T>& rhs) {
	Matrix <T> temp { rhs };
	temp *= number;
	//	No std::move here because of RVO
	return temp;
}

template <typename T>
Linear::Matrix <T> Linear::operator * (const Matrix <T>& lhs, const T number) {
	return number * lhs;
}

template <typename T>
Linear::Matrix <T> Linear::operator * (const Matrix <T>& lhs, const Matrix <T>& rhs) {
	Matrix <T> temp { lhs };
	temp *= rhs;
	//	No std::move here because of RVO
	return temp;
}