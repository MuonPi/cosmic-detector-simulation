/*
http://www.chesterproductions.net.nz/blogs/it/code/c-valarray-matrix-madness/11/

C++ Valarray Matrix Madness
by admin under c++, code

Numeric STL container valarray began life partially as an attempt to make C++ appealing to the supercomputing community. At the time the big thing in those big machines was, the ironically named, vector processing. However the valarray fell by the wayside as the people driving its development left the STL development group. Perhaps they realised that it didn’t really fit 100% with the STL, or maybe they just got sidetracked; who knows. But it is still useful, and here are a few reasons why:

    Can be used to write faster code for numeric spaces than possible with other STL types like the ubiquitous vector template.
    Offers the coder the possibility of staying within the comfort zone of the STL and familiar object oriented concepts with a small speed sacrifice over hand carved C.
    Allows library developers a way of writing optimized libraries for different environments so that coder can concentrate on the development at hand and not loose track in the complexity of the target environment.

So I’ve been playing around with valarray. It seemed that the best example to play with is that old classic the matrix. So here it is. Yes I know that there are some particularly hairy things wrong with it, but its not meant as a copy and paste solution. Its here as the results of a learning exercise and an example of what’s possible with valarray. There are one or two places which are not implemented or have not been tested, but you should be able to complete or test these by just looking at the rest of the examples. It should compile and run as is.

Have a look and have fun.
*/

#pragma once

#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>
#include <string>
#include <valarray>
#include <vector>

#include "algebra_types.h"

/*------------------------------------------------------------------------------
  matrix2d interface

  This is the template for the matrix class.  Its not that fancy and there could
  be a lot more added to it but that was not the point of the exercise.  It is
  also limited by the number of types supported by the underlying valarray data
  container.  From an OO point of view its really not that good either given
  that a lot of the manipulation or generator methods really don't need to have
  access directly to the data part of the construct.  However on a practical
  basis including those methods in this class provide small benefits in speed
  and allow things to be a little more obvious.  There are also problems in the
  way that the generators return new matrix2d objects.  But I've no intention to
  fix them as this was just a learning exercise.  I've tossed the implementation
  of the methods into seperate compilation objects to make the interface cleaner
  for inspection purposes.  In keeping with the valarray perspective there is no
  bounds checking anywhere - you have been warned.
  ----------------------------------------------------------------------------*/

template<typename T>

class matrix2d {
public:

  // creates based on the rows and data size
  matrix2d(std::size_t rows, std::valarray<T> data);
  // creates an empty rows x size matrix
  matrix2d(std::size_t rows, std::size_t cols);
  // direct initialisation - beware that rows x cols must equal data.size()
  matrix2d(std::size_t rows, std::size_t cols, std::valarray<T> data);

  // get the number of rows in the matrix2d
  std::size_t rows() const;
  // get the number of columns in the matrix2d
  std::size_t cols() const;
  // get a copy of the data in the matrix2d
  std::valarray<T> array() const;

  // retrieve the data from row r of the matrix
  std::valarray<T> row(std::size_t r) const;
  // retrieve the data from col c of the matrix
  std::valarray<T> col(std::size_t c) const;

  // retrieve reference to the data from row r of the matrix
  std::slice_array<T> row(std::size_t r);
  // retrieve reference to the data from col c of the matrix
  std::slice_array<T> col(std::size_t c);

  // basic item reference
  T & operator()(std::size_t r, std::size_t c);
  // basic item retrieval
  T operator()(std::size_t r, std::size_t c) const;

  friend matrix2d operator*(const matrix2d& A, const matrix2d& B)
  {
        assert(A.cols() == B.rows());
        matrix2d result { A.rows(), B.cols() };
        for(size_t i = 0; i < A.rows(); ++i) {
            for(size_t j = 0; j < B.cols(); ++j) {
                // Take dot product of row a[i] and col b[j]
                result(i,j) = (A.row(i)*B.col(j)).sum();
            }
        }
        return result;
  }

  friend matrix2d operator*(T c, const matrix2d& A)
  {
      matrix2d result { A };
      result.data_ *= c;
      return result;
  }

  friend matrix2d operator+(const matrix2d& A, const matrix2d& B)
  {
        assert(A.cols() == B.cols() && A.rows() == B.rows());
        matrix2d result { A.rows(), A.cols() };
        for(size_t i = 0; i < A.rows(); ++i) {
            for(size_t j = 0; j < B.cols(); ++j) {
                result(i,j) = A(i,j) + B(i,j);
            }
        }
        return result;
  }
  
  friend matrix2d operator-(const matrix2d& A, const matrix2d& B)
  {
        assert(A.cols() == B.cols() && A.rows() == B.rows());
        matrix2d result { A.rows(), A.cols() };
        for(size_t i = 0; i < A.rows(); ++i) {
            for(size_t j = 0; j < B.cols(); ++j) {
                result(i,j) = A(i,j) - B(i,j);
            }
        }
        return result;
  }

  // generate a new matrix that is the transposition of this one
  matrix2d<T> transpose();
  // generate a new matrix with this matrix's data and sort each row
  matrix2d<T> rowItmSort();
  // generate a new matrix with this matrix's data and sort each row in reverse
  matrix2d<T> rowItmSortRev();
  // generate a new matrix with this matrix's data and sort each col
  matrix2d<T> colItmSort();
  // generate a new matrix with this matrix's data and sort each col in reverse
  matrix2d<T> colItmSortRev();

  // generate a new matrix of this one with m appended below
  matrix2d<T> appendRows(matrix2d<T> &m);
  // generate a new matrix of this one with m appended to the right
  matrix2d<T> appendCols(matrix2d<T> &m);
  // generate a matrix of this one, upper left corner at row t col l - UNTESTED
  matrix2d<T> extractMatrix2d(size_t t, size_t l, size_t w, size_t h);

protected:

  std::size_t rows_;
  std::size_t cols_;
  std::valarray<T> data_;

};

/*------------------------------------------------------------------------------
  matrix2d implementation
  ----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  matrix2d constructors
  ----------------------------------------------------------------------------*/

template<class T>
matrix2d<T>::matrix2d(std::size_t rows, std::valarray<T> data) :
rows_(rows), cols_(data.size() / rows), data_(data) {}

template<class T>
matrix2d<T>::matrix2d(std::size_t rows, std::size_t columns) :
rows_(rows), cols_(columns), data_(rows * columns) {}

template<class T>
matrix2d<T>::matrix2d(std::size_t rows, std::size_t columns,
std::valarray<T> data) :
rows_(rows), cols_(columns), data_(data) {}

/*------------------------------------------------------------------------------
  matrix2d operations
  ----------------------------------------------------------------------------*/

template<class T>
std::size_t matrix2d<T>::rows() const {
  return rows_;
}

template<class T>
std::size_t matrix2d<T>::cols() const {
  return cols_;
}

template<class T>
std::valarray<T> matrix2d<T>::array() const {
  return data_;
}

template<class T>
std::valarray<T> matrix2d<T>::row(std::size_t r) const {
  return data_[std::slice(r * cols(), cols(), 1)];
}

template<class T>
std::valarray<T> matrix2d<T>::col(std::size_t c) const {
  return data_[std::slice(c, rows(), cols())];
}

template<class T>
std::slice_array<T> matrix2d<T>::row(std::size_t r) {
  return data_[std::slice(r * cols(), cols(), 1)];
}

template<class T>
std::slice_array<T> matrix2d<T>::col(std::size_t c) {
  return data_[std::slice(c, rows(), cols())];
}

template<class T>
T& matrix2d<T>::operator()(std::size_t r, std::size_t c) {
  return data_[r * cols() + c];
}

template<class T>
T matrix2d<T>::operator()(std::size_t r, std::size_t c) const {
  return row(r)[c];
}

/*------------------------------------------------------------------------------
  matrix2d generators
  ----------------------------------------------------------------------------*/

template<class T>
matrix2d<T> matrix2d<T>::transpose() {
    matrix2d<T> result(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        result.col(i) = static_cast<std::valarray<T> > (row(i));
    }
    return result;
}

template<class T>
matrix2d<T> matrix2d<T>::rowItmSort() {
    matrix2d<T> result(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        std::valarray<T> x = static_cast<std::valarray<T> > (row(i));
        std::sort(&x[0], &x[cols_]);
        result.row(i) = x;
    }
    return result;
}

template<class T> bool rev (const T & a, const T & b) { return a > b; }

template<class T>
matrix2d<T> matrix2d<T>::rowItmSortRev() {
    matrix2d<T> result(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        std::valarray<T> x = static_cast<std::valarray<T> > (row(i));
        std::sort(&x[0], &x[cols_], rev<T>);
        result.row(i) = x;
    }
    return result;
}

template<class T>
matrix2d<T> matrix2d<T>::colItmSort() {
    matrix2d<T> result(rows_, cols_);
    for (std::size_t i = 0; i < cols_; ++i) {
        std::valarray<T> x = static_cast<std::valarray<T> > (col(i));
        std::sort(&x[0], &x[rows_]);
        result.col(i) = x;
    }
    return result;
}

template<class T>
matrix2d<T> matrix2d<T>::colItmSortRev() {
    matrix2d<T> result(rows_, cols_);
    for (std::size_t i = 0; i < cols_; ++i) {
        std::valarray<T> x = static_cast<std::valarray<T> > (col(i));
        std::sort(&x[0], &x[rows_], rev<T>);
        result.col(i) = x;
    }
    return result;
}

template<class T>
matrix2d<T> matrix2d<T>::appendRows(matrix2d<T> &m) {
    matrix2d<T> result(rows_ + m.rows_, cols_);
    result.data_[std::slice(0, rows_ * cols_, 1)] = data_;
    result.data_[std::slice(rows_ * cols_, m.rows_ * m.cols_, 1)] = m.data_;
    return result;
}

template<class T>
matrix2d<T> matrix2d<T>::appendCols(matrix2d<T> &m) {
    matrix2d<T> result(rows_, cols_ + m.cols_);
    std::size_t s1[] = {rows_,cols_}; // shape of left matrix
    std::size_t p1[] = {result.cols_,1}; // position of left matrix in result
    std::size_t s2[] = {m.rows_,m.cols_}; // shape of right matrix
    std::size_t p2[] = {result.cols_,1}; // position or right matrix in result
    std::valarray<std::size_t> sv1(s1, 2);
    std::valarray<std::size_t> pv1(p1, 2);
    std::valarray<std::size_t> sv2(s2, 2);
    std::valarray<std::size_t> pv2(p2, 2);
    result.data_[std::gslice(0, sv1, pv1)] = data_; // copy left matrix into place
    result.data_[std::gslice(cols_, sv2, pv2)] = m.data_; // repeat for m
    return result;
}

template<class T>
matrix2d<T> matrix2d<T>::extractMatrix2d(size_t x, size_t y, size_t w,
size_t h) {

  /* TEST ME TEST ME TEST ME TEST ME TEST ME TEST ME TEST ME TEST ME TEST ME */

  matrix2d<T> result(h, w);

  size_t x2[] = {h, w}, s[] = {rows_, 1};
  std::valarray<size_t> xa(x2, 2), sa(s, 2);

  result.data_ = data_[(const std::gslice)std::gslice(y * rows_ + x, xa, sa)];

  return result;
}
