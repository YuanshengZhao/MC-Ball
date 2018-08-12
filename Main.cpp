#include "Func-Def.h"

const int N = 4 * LENTH * LENTH * LENTH;
unsigned int cnt[NMESH];
const real pi = 3.141592653589793238462643383279502884197;
real deltar = 1.0, rlim = .000005, rlimsqr = rlim * rlim;
real ergg[N], rdfdata[NMESH*2], *rdf=&rdfdata[NMESH];
real bx, ul, lvg;

VEC3 ball[N];

real distsqr(int ii, int jj)
{
	real dx = ball[ii].sx - ball[jj].sx,
		dy = ball[ii].sy - ball[jj].sy,
		dz = ball[ii].sz - ball[jj].sz;

	if (dx > ul) dx = bx - dx;
	else if (dx < -ul) dx = bx + dx;

	if (dy > ul) dy = bx - dy;
	else if (dy < -ul) dy = bx + dy;

	if (dz > ul) dz = bx - dz;
	else if (dz < -ul) dz = bx + dz;

	return dx * dx + dy * dy + dz * dz;
}

//L-J Potential
real PtlErg(real rsqr)
{
	real r_s6 = Cub(rlimsqr / (rsqr + 1.e-10));
	return r_s6 * r_s6 - 2 * r_s6;
}

real DEdr(real rs)
{
    real r_s = rlimsqr / (rs*rs + 1.e-10),
    r_s2=r_s * r_s,
    r_s4=r_s2*r_s2,
    r_s7=r_s4*r_s2*r_s;
    return 12*rs/rlimsqr*(r_s4-r_s7);
}


real Erg(int ii)
{
	real erg = 0;
	for (int i = 0; i < N; i++)
	{
		if (i == ii) continue;
		erg += PtlErg(distsqr(i, ii));
	}
	return erg;
}

int bmove(int step, real beta)
{
	int nn, nmv = 0, nsp = 0;
	real temp;
	VEC3 dr, tp;
	step *= N;
	while (nsp < step)
	{
		nsp++;
		nn = Rnd() % N;
		dr.RandVec(lvg);
		tp = ball[nn];
		ball[nn] += dr;
		if (ball[nn].sx < 0) ball[nn].sx += bx;
		if (ball[nn].sx > bx) ball[nn].sx -= bx;
		if (ball[nn].sy < 0) ball[nn].sy += bx;
		if (ball[nn].sy > bx) ball[nn].sy -= bx;
		if (ball[nn].sz < 0) ball[nn].sz += bx;
		if (ball[nn].sz > bx) ball[nn].sz -= bx;

		temp = Erg(nn);

		if (Rnd01() > exp(beta*(ergg[nn] - temp)))//discard change
			ball[nn] = tp;
		else//accept change
		{
			nmv++;
			ergg[nn] = temp;
		}
	}

	return nmv;
}

void Init()
{
	bx = LENTH;
	ul = bx / 2.0;
	deltar = .6 * bx / NMESH;
	lvg = .005;

	int cct = 0;
	//Set FCC Lattice
	for (int ii = 0; ii < LENTH; ii++)
		for (int jj = 0; jj < LENTH; jj++)
			for (int kk = 0; kk < LENTH; kk++)
			{
				ball[cct++].Set(ii + .001, jj + .001, kk + .001);
				ball[cct++].Set(ii + .5 + .001, jj + .5 + .001, kk + .001);
				ball[cct++].Set(ii + .5 + .001, jj + .001, kk + .5 + .001);
				ball[cct++].Set(ii + .001, jj + .5 + .001, kk + .5 + .001);
			}

	rlim = .499*sqrt(2.0);
	rlimsqr = Sqr(rlim);

	for (int i = 0; i < N; i++)
		ergg[i] = Erg(i);

    for(int i=0;i<NMESH;i++)
        rdfdata[i]=(i+0.5)*deltar;
}

void smallerball(real ddx)
{
	rlim -= ddx;
	rlimsqr = rlim * rlim;
}

void stt()
{
	int nx;
	for (int i = 0; i < N; i++)
		for (int j = i + 1; j < N; j++)
		{
			nx = (int)(sqrt(distsqr(i, j)) / deltar);
			if (nx < NMESH) cnt[nx]++;
		}
}

real test(bool b_out, real beta, int step = 50)
{
	FILE *fp;
	real s = pi * N *pow(rlim / bx, 3) / 6;
	for (int i = 0; i < NMESH; i++) cnt[i] = 0;
    real toterg;

	for (int i = 0; i < NWARMING; i++)
	{
        toterg=0;
        for(int i=0;i<N;i++)
            toterg+=ergg[i];
		printf("Warming up... %d/%d, erg=%lf, nmv=%d\n", i + 1, NWARMING, toterg, bmove(step, beta));
	}

	for (int i = 0; i < NSTP; i++)
	{
        toterg=0;
        for(int i=0;i<N;i++)
            toterg+=ergg[i];
		printf("Counting... %d/%d, erg=%lf, nmv=%d\n", i + 1, NSTP, toterg, bmove(step, beta));
		stt();
	}
	s = 2 * bx * bx * bx * 3 / (N * N * 4 * pi * NSTP);

	for (int i = 1; i < NMESH - 1; i++)
		rdf[i] = cnt[i] / (Cub(deltar * (i + 1)) - Cub(deltar * i)) * s;

	if (b_out)
	{
		fp = fopen("rst.csv", "w");
		for (int i = 1; i < NMESH - 1; i++)
			fprintf(fp, "%f,", rdf[i]);
		fprintf(fp, "%f\n", cnt[NMESH - 1] / (pow(deltar * (NMESH), 3) - pow(deltar*(NMESH - 1), 3))*s);
		fclose(fp);
	}

	return 0;
}

int main()
{
	Init();
	//for (real tp = .001; tp < .011; tp += .003)
	//{
	//	test(false, 1.0 / tp);
	//	cout << tp << endl;
	//	ShowGraph();
	//}
    real tp=.03;
    test(false, 1.0 / tp,10);
    cout << tp << endl;
    ShowGraph(rdfdata,NMESH*2,1);
//
//    real tp=.1;
//    test(false, 1.0 / tp,0);
//    cout << tp << endl;
//    ShowGraph(rdfdata,NMESH*2,1);

//    tp=.001;
//    test(false, 1.0 / tp,200);
//    cout << tp << endl;
//    ShowGraph(rdf,NMESH);
//    real (*func[2])(real);
//    func[0]=PtlErg;
//    func[1]=DEdr;
//    FuncPlot(func,1,rlim/1.2,rlim*2);
	return 0;
}
