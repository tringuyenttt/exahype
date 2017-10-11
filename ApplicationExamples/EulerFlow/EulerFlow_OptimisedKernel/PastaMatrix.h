#ifndef SMATRIX_H
#define SMATRIX_H

/**
 * Pasta is a small helper library which resembles parts of the closed-source
 * Pizza Library which is used in the Whisky code by means of an Open Source
 * clone.
 * 
 * Currently, we only implement basic header-only vector classes in order to
 * overcome obstacles from Peano's tarch. We also offer bindings to the Vector
 * classes there!
 * 
 ********************************* ATTENTION ********************************
 * Currently, we copy Pizza source code here. This is should be changed soon.
 * We could switch over to alternative OSS code doing something similar, eg.
 * https://bitbucket.org/whiskydevs/tensortemplates/
 ****************************************************************************
 * 
 * I also would like to have classes here which allow an in memory vector
 * representation only for compile time of stuff like an existing double*.
 * Thus we would avoid unnessessarily copying stuff.
 *
 * This library is (CC0) / Public Domain.
 **/

#include <complex>
#include <initializer_list>
#include <assert.h>

namespace Pasta {

// PizzaNumUtils/src/mathe.h
typedef double pz_real;
typedef std::complex<pz_real> pz_cmplx;

template <class T> inline T sqr(T t) {return t*t;}

// PizzaNumUtils/src/smatrix.h

typedef pz_real sm_real;
template<int N,bool UP> class sm_tensor2_sym;

//---------------------------------------------------------------------------------------
//  Helper
//---------------------------------------------------------------------------------------
enum zero_literal {ZERO=0}; //we want to write matrix=ZERO but not matrix= 14.0
enum one_literal {ONE=1};

template<int N> class sm_array {
  typedef sm_array<N> me;
  public:
  sm_real v[N];
  sm_array(){}
  sm_real &operator[](int j) {return v[j];}
  const sm_real &operator[](int j) const {return v[j];}

  void assign_prod(const me &a,const sm_real z)
    {for(int i=0;i<N;i++) v[i]=a.v[i]*z;}
  void assign_div(const me &a,sm_real z)
    {assign_prod(a,1.0/z);}
  void assign_sum(const me &a,const me &b)
    {for(int i=0;i<N;i++) v[i]=a.v[i]+b.v[i];}
  void assign_diff(const me &a,const me &b)
    {for(int i=0;i<N;i++) v[i]=a.v[i]-b.v[i];}
  void assign_minus(const me &a)
    {for(int i=0;i<N;i++) v[i]=-a.v[i];}

  void operator+=(const me &a) {assign_sum(*this,a);}
  void operator-=(const me &a) {assign_diff(*this,a);}
  void operator*=(sm_real z) {assign_prod(*this,z);}
  void operator/=(sm_real z) {assign_div(*this,z);}

  void zero() {for(int i=0;i<N;i++) v[i]=0.0;}
};

//---------------------------------------------------------------------------------------
//  N-vector, co/contra variant
//---------------------------------------------------------------------------------------

template<int N,bool UP> class sm_tensor1 {
  typedef sm_tensor1<N,UP> me;
  public:
  sm_array<N> c;
  enum {SIZE=N};
  sm_tensor1(){}
  sm_tensor1(zero_literal z) {c.zero();}
  sm_real &operator()(int j) {return c[j];}
  const sm_real &operator()(int j) const {return c[j];}

  // added by sven. This is of course stupid.
  sm_tensor1(const sm_real* const prefill) {
	  for(int i=0; i<N; i++) c[i]=prefill[i];
  }
  // this requires C++ and allows writing sm_tensor1{1,2,3} or sm_tensor({1,2,3}) but not sm_tensor(1,2,3).
  sm_tensor1(std::initializer_list<sm_real> prefill) {
	assert(prefill.size() == N);
	int i=0;
	for (auto it : prefill) //  for (auto it = prefill.begin() ; it != prefill.end() ; ++it)
		c[i++] = it;
  }

  void assign_sum(const me &a,const me &b) {c.assign_sum(a.c,b.c);}
  void assign_diff(const me &a,const me &b) {c.assign_diff(a.c,b.c);}
  void assign_minus(const me &a) {c.assign_minus(a.c);}
  void assign_prod(const me &a,sm_real z) {c.assign_prod(a.c,z);}
  void assign_div(const me &a,sm_real z) {c.assign_div(a.c,z);}

  void operator+=(const me &a) {assign_sum(*this,a);}
  void operator-=(const me &a) {assign_diff(*this,a);}
  void operator*=(sm_real z) {assign_prod(*this,z);}
  void operator/=(sm_real z) {assign_div(*this,z);}

  me operator+(const me &a) const {me e; e.assign_sum(*this,a); return e;}
  me operator-(const me &a) const {me e; e.assign_diff(*this,a); return e;}
  me operator*(sm_real z) const {me e; e.assign_prod(*this,z); return e;}
  me operator/(sm_real z) const {me e; e.assign_div(*this,z); return e;}

  void assign_prod(const sm_tensor2_sym<N,UP> &m,const sm_tensor1<N,!UP> &w);
  sm_real norm() const;
};

template<int N,bool UP>
sm_real sm_tensor1<N,UP>::norm() const
{
  sm_real e=c[0]*c[0];
  for (int i=1;i<N;i++) e+=c[i]*c[i];
  return sqrt(e);
}

template<int N,bool UP>
sm_tensor1<N,UP> operator-(sm_tensor1<N,UP> &v)
{
  sm_tensor1<N,UP> erg;
  erg.assign_minus(v);
  return erg;
}
//---------------------------------------------------------------------------------------
// Contraction v^i w_i  resp. v_i w^i
//---------------------------------------------------------------------------------------

template<int N,bool UP>
inline sm_real operator*(const sm_tensor1<N,UP> &v,const sm_tensor1<N,!UP> &w)
{
  sm_real erg=v(0)*w(0);
  for(int i=1;i<N;i++) erg+=v(i)*w(i);
  return erg;
}

//---------------------------------------------------------------------------------------
// Scalar * vector
//---------------------------------------------------------------------------------------

template<int N,bool UP>
sm_tensor1<N,UP> inline operator*(sm_real z,const sm_tensor1<N,UP> &a) {
  sm_tensor1<N,UP> erg;
  erg.assign_prod(a,z);
  return erg;
}


// Standard 2D/3D variables for doing non-relativistic math
typedef sm_tensor1<3,true> vec3;
typedef sm_tensor1<2,true> vec2;



} // namespace Pasta
#endif /* SMATRIX_H */