#ifndef _FUNC_DEF_H_H_H_H_
#define _FUNC_DEF_H_H_H_H_

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <time.h>
#include <cmath>

using namespace std;
typedef double real;

#define LENTH 5
#define NMESH 503
#define NSTP 101
#define NWARMING 301

extern const real pi;
extern const int N;
extern real ergg[], rdfdata[], *rdf;
extern real deltar, rlim, rlimsqr;


class VEC3
{
public:
	VEC3();
	VEC3(real xx, real yy, real zz);
	~VEC3();
	real sx, sy, sz;
	void RandVec(real);
	void operator+=(const VEC3 &);
	void operator-=(const VEC3 &);
	void Set(real, real, real);
};

unsigned int Rnd();

real Sqr(real);
real Cub(real);
unsigned int Rnd();
real Rnd01();
void ShowGraph(real *dat, int npts, int pls = 0);
void FuncPlot(real (**func)(real), int pls, real xstt, real xedd, int npts = 50);
real PtlErg(real rsqr);
real DEdr(real rsqr);
real ItgrtRDF();
#endif
