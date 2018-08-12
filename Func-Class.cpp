#include "Func-Def.h"

unsigned int seed0 = (int)time(NULL), seed1 = 0, seed2 = seed0;

real Sqr(real z)
{
	return z * z;
}

real Cub(real z)
{
	return z * z*z;
}

unsigned int Rnd()
{
	return seed0 = (seed0 * 40014u) % 2147483563u;
}

real Rnd01()
{
	seed1 = (seed1 * 40014u) % 2147483563u;
	seed2 = (seed2 * 40692u) % 2147483399u;
	return seed1 < seed2 ?
		((real)(2147483562u + seed1 - seed2)) / ((real)2147483563u)
		: ((real)(seed1 - seed2)) / ((real)2147483563u);
}

void VEC3::Set(real xx, real yy, real zz)
{
	sx = xx; sy = yy; sz = zz;
}

void VEC3::RandVec(real len)
{
	real zs = 2.0 * Rnd01() - 1.0,
		xy = sqrt(1 - zs * zs),
		theta = 2 * pi * Rnd01(),
		lth = len * Rnd01();
	sz = lth * zs;
	sx = lth * xy * cos(theta);
	sy = lth * xy * sin(theta);
}

void VEC3::operator+=(const VEC3 &xx)
{
	sx += xx.sx; sy += xx.sy; sz += xx.sz;
}

void VEC3::operator-=(const VEC3 &xx)
{
	sx -= xx.sx; sy -= xx.sy; sz -= xx.sz;
}

VEC3::VEC3()
{
}

VEC3::VEC3(real xx, real yy, real zz)
{
	Set(xx, yy, zz);
}

VEC3::~VEC3()
{
}

void ShowGraph(real *dat, int npts, int pls)
{
	char cmd[npts * 20 + 80];
	int pos;
	if(pls<1)
	{
#ifdef _WIN32
        pos = sprintf(cmd, "python arg-plot.py %f", dat[0]);
#else
        pos = sprintf(cmd, "python3 arg-plot.py %f", dat[0]);
#endif
        for (int k = 1; k < npts; k++)
            pos += sprintf(cmd + pos, " %f", dat[k]);
    }
    else
    {
#ifdef _WIN32
        pos = sprintf(cmd, "python arg-plot.py x %d %f", pls, dat[0]);
#else
        pos = sprintf(cmd, "python3 arg-plot.py x %d %f", pls, dat[0]);
#endif
        for (int k = 1; k < npts; k++)
            pos += sprintf(cmd + pos, " %f", dat[k]);
    }
	system(cmd);
}

void FuncPlot(real (**func)(real), int pls, real xstt, real xedd, int npts)
{
    real *das=new real[(pls+1)*npts];
    for(int ii=0;ii<npts;ii++)
        das[ii]=xstt+ii * (xedd-xstt)/(npts-1);
    int pos=npts;
    for(int jj=0;jj<pls;jj++)
        for(int ii=0;ii<npts;ii++)
            das[pos++]=(*func[jj])(das[ii]);
    ShowGraph(das,npts*(pls+1),pls);
    delete[] das;
}

real ItgrtRDF()
{
    real rst=0,rr;
    for(int ks=0;ks<NMESH;ks++)
    {
        rr=(ks+.5)*deltar;
        rst+=(rdf[ks]*rr*DEdr(rr));
    }
    rst*=deltar;
    return rst;
}
