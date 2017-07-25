/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#ifndef _EXAHYPE_KERNELS_KERNEL_UTILS_H_
#define _EXAHYPE_KERNELS_KERNEL_UTILS_H_

#include "../../Peano/tarch/Assertions.h"

#include <sstream>

namespace kernels {

// If you have a C++11-enabled compiler, use this:
//#define C11CONSTEXPR constexpr
// else this:
#define C11CONSTEXPR

/**
 * This is a single successor class for the idx2, idx3, idx4, idx5, idx6 classes.
 * It works basically like idx6. If you work with less than 6 dimensions, nothing
 * will change and the compiler optimizes away the unused variables.
 * 
 * In contrast to the index classes below, this new index class can also compute
 * the inverse index.
 **/
struct index {
	// add further indices if you want to support more than maximal 6 indices.
	const int i0, i1, i2, i3, i4, i5; // Lenghts
	const int b0, b1, b2, b3, b4, b5; // Basis
	const int size;
	
	index(int j0=1, int j1=1, int j2=1, int j3=1, int j4=1, int j5=1) :
		i0(j0), i1(j1), i2(j2), i3(j3), i4(j4), i5(j5),
		b0(i1 * i2 * i3 * i4 * i5),
		b1(i2 * i3 * i4 * i5),
		b2(i3 * i4 * i5),
		b3(i4 * i5),
		b4(i5),
		b5(1),
		size(j0*j1*j2*j3*j4*j5) {}
	
	/**
	 * Compute a single index ("superindex", global index, ...) from the tuples.
	 **/
	C11CONSTEXPR int get(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
		assertion2(j0 < i0, j1, i1);
		assertion2(j1 < i1, j2, i2);
		assertion2(j2 < i2, j3, i3);
		assertion2(j3 < i3, j4, i4);
		assertion2(j4 < i4, j5, i5);
		
		return j0 * b0 + j1 * b1 + j2 * b2 + j3 * b3 + j4 * b4 + j5 * b5;
	}
	
	/**
	 * Inverse index: Get the index tuple for a given index. This is the inverse
	 * function of get(), so such a call will not do anything:
	 *    myidx.rev(myidx.get(a,b,c...), a, b, c, ...)
	 **/
	void rev(int pos, int& j0, int& j1, int& j2, int& j3, int& j4, int& j5) const {
		// The algorithm exploits the short notation for int a,b that
		// a/b == std::floor(a/b). The c's are the contributions at position
		// i while the j's are the digits for position i. If you want to
		// increase the overall amount of positions, extend logically.
		j0 = (pos) / b0;
		j1 = (pos-j0*b0) / b1;
		j2 = (pos-j0*b0-j1*b1) / b2;
		j3 = (pos-j0*b0-j1*b1-j2*b2) / b3;
		j4 = (pos-j0*b0-j1*b1-j2*b2-j3*b3) / b4;
		j5 = (pos-j0*b0-j1*b1-j2*b2-j3*b3-j4*b4) / b5;
	}
	
	/* As there are no reference default values, ie. we cannot write
	 *   void rev(int pos, int& j0=NULL) or
	 *   void rev(int pos, int& j1=int(), ...)
	 * we do it with pointers. Now you can write instead
	 *   int i, j, k;
	 *   myidx.rev(position, &i, &j, &k);
	 * which is much more convenient than tracking the non-used indices
	 * for your own.
	 */
	void rev(int pos, int* const j0=NULL, int* const j1=NULL, int* const j2=NULL, int* const j3=NULL, int* const j4=NULL, int* const j5=NULL) const {
		int a0=0, a1=0, a2=0, a3=0, a4=0, a5=0; // default storages
		rev(pos, j0?*j0:a0, j1?*j1:a1, j2?*j2:a2, j3?*j3:a3, j4?*j4:a4, j5?*j5:a5);
	}
	
	// syntactic sugar:
	C11CONSTEXPR int operator()(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
		return get(j0, j1, j2, j3, j4, j5);
	}
	
	/**
	 * Checks if the given indices are in the valid range. Works also if assertions are
	 * not enabled. Can be handy to check access to variables.
	 **/
	C11CONSTEXPR bool check(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
		// assertion way
		assertion2(j0 < i0, j1, i1);
		assertion2(j1 < i1, j2, i2);
		assertion2(j2 < i2, j3, i3);
		assertion2(j3 < i3, j4, i4);
		assertion2(j4 < i4, j5, i5);
		// non-assertion way:
		return (j0 < i0) && (j1 < i1) && (j2 < i2) && (j3 < i3) && (j4 < i4) && (j5 < i5);
	}
	
	/**
	 * Allows to print index tuples.
	 **/
	static std::string strIndex(int min=0, int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) {
		std::stringstream s;
		s << "(";
		s << j0;
		if(j1 > min) s << ", " << j1;
		if(j2 > min) s << ", " << j2;
		if(j3 > min) s << ", " << j3;
		if(j4 > min) s << ", " << j4;
		if(j5 > min) s << ", " << j5;
		s << ")";
		return s.str();
	}
	
	/// Some index to string
	C11CONSTEXPR std::string getStr(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
		return strIndex(/*min*/ -1, j0,j1,j2,j3,j4,j5);
	}
	C11CONSTEXPR std::string revStr(int pos) const {
		int j0,j1,j2,j3,j4,j5;
		rev(pos, j0,j1,j2,j3,j4,j5);
		return strIndex(/*min*/ -1, j0,j1,j2,j3,j4,j5);
	}

	/// Object to string
	C11CONSTEXPR std::string toString() const {
		return strIndex(/*min*/ +1, i0,i1,i2,i3,i4,i5);
	}
};

/**
 * A tiny class for storing a continous n-dimensional array together with its sizes.
 * You can call instances as getter and setter for elements, ie.
 * 
 * 	array<double> L(2,3); // a 2x3 matrix of doubles
 *	// use as setter
 *	L(1,2) = 7;
 *	// or as getter:
 *	double val = L(1,2);
 * 
 **/
template<typename T>
struct array {
	index idx;
	T *data;
	
	array(int j0=1, int j1=1, int j2=1, int j3=1, int j4=1, int j5=1) : idx(j0,j1,j2,j3,j4,j5) {
		data = new T[idx.size];
	}
	
	T& operator()(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
		return data[idx(j0, j1, j2, j3, j4, j5)];
	}

	~array() { delete data; }
};

/**
 * The same as array, just without local storage but instead shadowing data which is not owned.
 * This can be handy for accessing big arrays stored somewhere else. A tiny example is
 * 
 *	double storage[5*7*9];
 *	dshadow convenient(storage, 5, 7, 9);
 *	convenient(2,3,8) = 10;
 *      double val = convenient(2,3,8);
 *
 **/
template<typename T>
struct shadow {
	index idx;
	T *data;
	
	shadow(T* _data, int j0=1, int j1=1, int j2=1, int j3=1, int j4=1, int j5=1) : data(_data), idx(j0,j1,j2,j3,j4,j5) {}
	
	T& operator()(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
		return data[idx(j0, j1, j2, j3, j4, j5)];
	}
};

typedef array<double> darray;
typedef shadow<double> dshadow;


struct idx2 {
  idx2(int I, int J, int line = -1) : I_(I), J_(J), size(I*J), line_(line) {}

  int operator()(int i, int j) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    return i * J_ + j;
  }

  void rev(int pos, int& i, int& j) const {
    assertion(pos < size);
    i = pos % J_;
    j = pos - i * J_;
  }

  const int I_, J_, size, line_;
};

struct idx3 {
  idx3(int I, int J, int K, int line = -1) : I_(I), J_(J), K_(K), size(I*J*K), line_(line) {}

  int operator()(int i, int j, int k) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    assertion3(k < K_, k, K_, line_);
    return i * (J_ * K_) + j * K_ + k;
  }

  const int I_, J_, K_, size, line_;
};

struct idx4 {
  idx4(int I, int J, int K, int L, int line = -1)
      : I_(I), J_(J), K_(K), L_(L), size(I*J*K*L), line_(line) {}

  int operator()(int i, int j, int k, int l) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    assertion3(k < K_, k, K_, line_);
    assertion3(l < L_, l, L_, line_);
    return i * (J_ * K_ * L_) + j * (K_ * L_) + k * L_ + l;
  }

  const int I_, J_, K_, L_, size, line_;
};

struct idx5 {
  idx5(int I, int J, int K, int L, int M, int line = -1)
      : I_(I), J_(J), K_(K), L_(L), M_(M), size(I*J*K*L*M), line_(line) {}

  int operator()(int i, int j, int k, int l, int m) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    assertion3(k < K_, k, K_, line_);
    assertion3(l < L_, l, L_, line_);
    assertion3(m < M_, m, M_, line_);
    return i * (J_ * K_ * L_ * M_) + j * (K_ * L_ * M_) + k * (L_ * M_) +
           l * M_ + m;
  }

  const int I_, J_, K_, L_, M_, size, line_;
};

struct idx6 {
  idx6(int I, int J, int K, int L, int M, int N, int line = -1)
      : I_(I), J_(J), K_(K), L_(L), M_(M), N_(N), size(I*J*K*L*M*N), line_(line) {}

  int operator()(int i, int j, int k, int l, int m, int n) {
    assertion3(i < I_, i, I_, line_);
    assertion3(j < J_, j, J_, line_);
    assertion3(k < K_, k, K_, line_);
    assertion3(l < L_, l, L_, line_);
    assertion3(m < M_, m, M_, line_);
    assertion3(n < N_, n, N_, line_);
    return i * (J_ * K_ * L_ * M_ * N_) + j * (K_ * L_ * M_ * N_) +
           k * (L_ * M_ * N_) + l * (M_ * N_) + m * N_ + n;
  }

  const int I_, J_, K_, L_, M_, N_, size, line_;
};

}  // namespace kernels

#endif  // _EXAHYPE_KERNELS_KERNEL_UTILS_H_
