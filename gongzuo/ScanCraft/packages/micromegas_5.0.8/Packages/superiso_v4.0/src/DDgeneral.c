#include "include.h"
#include "include_dm.h"

void DDbraket(int n, double (*func)(double[], double[]), double xtra[], double* a0, double*x0, double*b0, double xinit[], double xi[], double xlimmin[], double xlimmax[], double factor)
{
	int BRAKETMAX=50;
	int nmax=5;

	int i;
	double a0tmp,b0tmp;
	
	double xr,xrand[nmax],frand;
	double xa[nmax], xb[nmax];
	double fa,fb,fmin;

/*	*a0=-1.0e+99;*/
/*	*b0=1.0e+99;*/
	int afirst=1;
	int bfirst=1;
	for(i=0;i<n;i++)
	{
		if(xi[i]!=0)
		{
			a0tmp=(xlimmin[i]-xinit[i])/xi[i];
			b0tmp=(xlimmax[i]-xinit[i])/xi[i];
			
			if(a0tmp>*a0||afirst) {*a0=a0tmp; afirst=0;}
			if(b0tmp<*b0||bfirst) {*b0=b0tmp; bfirst=0;}
		}	
	}
	
	if(*a0>*b0)
	{
		*a0=0;
		*b0=0;
		*x0=0;
	}
	else
	{ 
		for(i=0;i<n;i++)
		{
			xa[i]=xinit[i]+*a0*xi[i];
			xb[i]=xinit[i]+*b0*xi[i];
		}
		fa=factor*func(xa,xtra);
		fb=factor*func(xb,xtra);
		if(fa<fb)
		{ 
			*x0=*a0;
			fmin=fa;
		}
		else
		{
			*x0=*b0;
			fmin=fb;
		}
		int irand;
		
		int icount=0;
		for(irand=0;irand<BRAKETMAX;irand++)
		{
			xr=*a0+irand*(*b0-*a0)/BRAKETMAX;
			for(i=0;i<n;i++)xrand[i]=xinit[i]+xr*xi[i];
			frand=factor*func(xrand,xtra);
			if(frand<fmin)
			{
				*x0=xr;
				fmin=frand;
				icount++;
			}
			if(icount==5) break;
		}
	}
	
	return;
}

/*--------------------------------------------------------------------*/

int DDbrentmethod(int n, double (*func)(double[], double[]), double xtra[], double xmin[], double xi[], double xlimmin[], double xlimmax[], double *fmin, double factor)
//minimize function func along vector xi
{
	double a0,b0,x0, xmin1d;
	double xitemp[n];
	int i;

	DDbraket(n,func,xtra,&a0,&x0,&b0,xmin,xi,xlimmin,xlimmax,factor);
	
	int test=DDbrentmethod1D(func,xtra,a0,  x0, b0, fmin, &xmin1d, factor, xmin, xi, n);
	for(i=0;i<n; i++) xmin[i]=xmin[i]+xmin1d*xi[i];
	
	return test;
}

/*--------------------------------------------------------------------*/

double DDbrentmethod1Dfunc(double (*func)(double[], double[]), double xtra[], double x, double xmini[], double xi[], int n)
{
	int i;
	double xitemp[n];
	
	for(i=0;i<n;i++) xitemp[i]=xmini[i]+x*xi[i];
	
	return func(xitemp, xtra);
}

int DDbrentmethod1D(double (*func)(double[], double[]), double xtra[], double a0, double x0,double b0, double *fmin, double *xmin, double factor, double xmini[], double xi[], int n)
{	
	int BrentIterMax=1000;
	double ZEPS=1.0e-70;
	double CGOLD=0.1;
	double tol=1.0e-3;

	double tol1, tol2;
	double a,b;
	
	if(a0<b0)
	{
		a=a0;
		b=b0;
	}
	else
	{
		a=b0;
		b=a0;
	}
	
	double x,w,v; // points used in the parabolic fit
	double u;
	double  pa1, pa2; // parabola parameters pa2*X^2+pa1*X+pa0
	double pamin; // parabola minimum 
	double xm;
	double fx,fw,fv, fu;
	double e=0, etmp,d;
	int iter=0;
	
	x=w=v=x0;
	fx=fw=fv=factor*DDbrentmethod1Dfunc(func,xtra,x0,xmini,xi,n);
	xm=0.5*(a+b);
	
	do
	{
		tol1=tol*fabs(x)+ZEPS ;
		tol2=e*tol1;

		if (fabs(e)>tol1)
		{
			pa2=((fw-fx)/(w-x)-(fv-fx)/(v-x))/(w-v); // parabolic fit
			pa1=(fv-fx)/(v-x) -pa2*(v+x);
			pamin=-0.5*pa1/pa2;
			etmp=e;
			e=d;
			if(fabs(x-pamin)>0.5*etmp||pamin<a||pamin>b)
			{// cannot use the parabolic fit
				if(x>=xm) e=a-x;
				else e=b-x;
				d=CGOLD*e;
			}
			else
			{
				d=pamin-x;
				u=pamin;
				if((u-a)<tol2||(b-u)<tol2)
				{// pamin too close to borders
					d=fabs(tol1);
					if(xm<x) d=-d;
				}
			}
		}
		else
		{
			if(x>=xm) e=a-x;
			else e=b-x;
			d=CGOLD*e;
		}	
		
		if(fabs(d)>=tol1) u=x+d;
		else
		{
			if(d<0) u=x-fabs(tol1);
			else u=x+fabs(tol1);
		}
		fu=factor*DDbrentmethod1Dfunc(func,xtra,u,xmini,xi,n);

		if(fu<=fx)
		{
			if(u>=x) a=x;
			else b=x;
			v=w;w=x;x=u;
			fv=fw;fw=fx;fx=fu;
		}
		else
		{
			if(u<x) a=u;
			else b=u;
			if(fu<=fw||w==x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if(fu<=fv||v==x||v==w)
			{
				v=u;
				fv=fu;
			}
		}

		xm=0.5*(a+b);

		iter++;
	}
	while(iter<BrentIterMax&&fabs(x-xm)+0.5*(b-a)> 2*tol1);
	
	*xmin=x;
	*fmin=fx;

	if(iter==BrentIterMax)
	{
		printf("WARNING: too many iterations in Brent's method, precision decreased\n");
		return 0;
	}
	else return 1;
}

/*----------------------------------------------------------------*/
#define debug 0

int DDpowellaux(int n, double (*func)(double[], double[]), double xtra[], double x0[], double xlimmin[], double xlimmax[], double* fmin,  double xmin[], double ftol, double factor)
{	
	double TINY=1.0e-70;
	int ITMAX=500;
	int nmax=5;

	int i,j,iter;
	double vect[nmax][nmax];
	//initiate set of vectors
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j) vect[i][j]=1;
			else vect[i][j]=0;
		}
	}	

	for(j=0;j<n;j++) xmin[j]=x0[j]; // initiate minimum position

	double vecti[nmax];// current vector in the iteration
	double fmintemp;
	int imax;
	double fdeltamax;
	*fmin=factor*func(x0,xtra);
	double fiter, xiter[nmax];
	int testbrent; 
	
	double meanvect[nmax], xextrapol[nmax];

	double fextrapol;
	double t;
	
	for(iter=0;iter<=ITMAX;iter++)
	{
			
		fiter=*fmin;
		for(j=0;j<n;j++) xiter[j]=xmin[j];
		fdeltamax=0;
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++) vecti[j]=vect[j][i];
			fmintemp=*fmin;
			testbrent=DDbrentmethod(n,func,xtra,xmin,vecti,xlimmin,xlimmax,fmin,factor);
			if(testbrent==0)
			{
				if(debug) printf("Brent method unsuccessful\n");
				return 0;
			}

			if((fmintemp-*fmin)>fdeltamax)
			{
				imax=i;
				fdeltamax=fmintemp-*fmin;	
			}
		}


		if(2.*(fiter-(*fmin)) <= ftol*(fabs(fiter)+fabs(*fmin))+TINY) return 0;
		
		if (iter == ITMAX)
		{
			if(debug) printf("Powell exceeding maximum iterations.\n");
			return 0;
		}
		
		for (j=0;j<n;j++)
		{ //check if it is worth to keep the mean deplacement vector
			meanvect[j]=xmin[j]-xiter[j]; 
			xextrapol[j]=xmin[j]+meanvect[j];
		}

		fextrapol=factor*func(xextrapol,xtra);
		if (fextrapol < fiter)
		{
			t=2.*(fiter-2.*(*fmin)+fextrapol)*pow(fiter-(*fmin)-fdeltamax,2.)-fdeltamax*pow(fiter-fextrapol,2.); 

			if (t < 0.)
			{
				DDbrentmethod(n,func,xtra,xmin,meanvect, xlimmin, xlimmax, fmin, factor);
				vect[j][imax]=vect[j][n]; 
				vect[j][n]=meanvect[j];
			}
		}

	}
	
	return 1;
}


/*--------------------------------------------------------------------*/

int DDpowell(int n, double (*func)(double[], double[]), double xtra[], double xlimmin[], double xlimmax[], double* fmin, double xmin[], double ftol, char* option)
{
	int nmax=5;
	
	srand((unsigned int)time(NULL));
	double x0[nmax];
	int i;

	double factor=0.;
	if(strcmp(option,"min")==0) factor=1.;
	else if (strcmp(option,"max")==0) factor=-1.;

	for(i=0;i<n;i++) x0[i]=(xlimmin[i]+xlimmax[i])/2.;
	int test=DDpowellaux(n,func,xtra,x0,xlimmin, xlimmax, fmin,  xmin, ftol, factor);
	if(!strcmp(option, "max")) *fmin=-*fmin; 
	
	return test;
}

/*--------------------------------------------------------------------*/

void DDsimpson_rec(double (*func)(double, double[], struct DDparameters*), double xtra[], struct DDparameters* DDparam, double *f, double *x, int *n, int n1, int n2, double *restot, double err)
{
	int Nmax=500;
	int i;

	if((*n)<Nmax&&(n1<n2)&&n2<=*n)
	{
		double res4=simpson_noerr(f,x,n1,n2);
		double xn2=x[n2];
		double xn1=x[n1];

		for(i=(*n);i>n2; i--)
		{
			x[i+4]=x[i];
			f[i+4]=f[i];
		}
		*n+=4;
		n2+=4;

		double delta=(xn2-xn1)/((double)(n2-n1));

		for(i=n2;i>=n1; i--)
		{
			if(!(i&1))
			{
				x[i]=x[(i+n1)/2];
				f[i]=f[(i+n1)/2];
			}
			else
			{
				x[i]=x[n1]+(i-n1)*delta;
				f[i]=func(x[i],xtra,DDparam);
			}
		}

		double res8=simpson_noerr(f,x,n1,n2);
		double err0;
		if(*restot==0) err0=fabs(res8-res4);
		else err0=fabs((res8-res4)/(*restot));
		*restot+=(res8-res4);

		if(err0>err) 
		{
			int n0=(*n);
			DDsimpson_rec(func,xtra,DDparam, f,x, n,n1,n1+4,restot,err);
			int dn=*n-n0;
			for(i=dn+n1+4;i<=*n; i++)
			{
				x[i-dn-4]=x[i];
				f[i-dn-4]=f[i];
			}
			*n-=dn+4;
			DDsimpson_rec(func,xtra,DDparam, f,x,n, n1,n1+4,restot,err);
		}
	}
	
	return;
}

/*--------------------------------------------------------------------*/

double DDsimpson(struct DDparameters* DDparam, double (*func)(double, double[], struct DDparameters*), double xtra[], double a, double b, double err)
{
	int n=4;
	double * f=malloc(500*sizeof(double));
	double * x=malloc(500*sizeof(double));
	int i;
	double delta=(b-a)/((double)n);
	for(i=0;i<=n; i++)
	{
		x[i]=a+i*delta;
		f[i]=(*func)(a+i*delta,xtra,DDparam);
	}
	double restot=simpson_noerr(f,x,0,n);
	DDsimpson_rec(func,xtra,DDparam, f,x, &n,0,n,&restot,err);
	free(f);
	free(x);

	return restot;
}

/*--------------------------------------------------------------------*/

void IDbraket(int n, double (*func)(double[], double[], struct array, struct propagation_parameters*, struct fermi*), double xtra[], struct array spect, struct propagation_parameters* pparam, struct fermi* fe, double* a0, double*x0, double*b0, double xinit[], double xi[], double xlimmin[], double xlimmax[], double factor)
{
	int BRAKETMAX=50;
	int nmax=5;

	int i;
	double a0tmp,b0tmp;
	
	double xr,xrand[nmax],frand;
	double xa[nmax], xb[nmax];
	double fa,fb,fmin;

/*	*a0=-1.0e+99;*/
/*	*b0=1.0e+99;*/
	int afirst=1;
	int bfirst=1;
	for(i=0;i<n;i++)
	{
		if(xi[i]!=0)
		{
			a0tmp=(xlimmin[i]-xinit[i])/xi[i];
			b0tmp=(xlimmax[i]-xinit[i])/xi[i];
			
			if(a0tmp>*a0||afirst) {*a0=a0tmp; afirst=0;}
			if(b0tmp<*b0||bfirst) {*b0=b0tmp; bfirst=0;}
		}	
	}
	
	if(*a0>*b0)
	{
		*a0=0;
		*b0=0;
		*x0=0;
	}
	else
	{ 
		for(i=0;i<n;i++)
		{
			xa[i]=xinit[i]+*a0*xi[i];
			xb[i]=xinit[i]+*b0*xi[i];
		}
		fa=factor*func(xa,xtra,spect,pparam,fe);
		fb=factor*func(xb,xtra,spect,pparam,fe);
		if(fa<fb)
		{ 
			*x0=*a0;
			fmin=fa;
		}
		else
		{
			*x0=*b0;
			fmin=fb;
		}
		int irand;
		
		int icount=0;
		for(irand=0;irand<BRAKETMAX;irand++)
		{
			xr=*a0+irand*(*b0-*a0)/BRAKETMAX;
			for(i=0;i<n;i++)xrand[i]=xinit[i]+xr*xi[i];
			frand=factor*func(xrand,xtra,spect,pparam,fe);
			if(frand<fmin)
			{
				*x0=xr;
				fmin=frand;
				icount++;
			}
			if(icount==5) break;
		}
	}
	
	return;
}

/*--------------------------------------------------------------------*/

int IDbrentmethod(int n, double (*func)(double[], double[], struct array, struct propagation_parameters*, struct fermi*), double xtra[], struct array spect, struct propagation_parameters* pparam, struct fermi* fe, double xmin[], double xi[], double xlimmin[], double xlimmax[], double *fmin, double factor)
//minimize function func along vector xi
{
	double a0,b0,x0, xmin1d;
	double xitemp[n];
	int i;

	IDbraket(n,func,xtra,spect,pparam,fe,&a0,&x0,&b0,xmin,xi,xlimmin,xlimmax,factor);
	
	int test=IDbrentmethod1D(func,xtra,spect,pparam,fe,a0,  x0, b0, fmin, &xmin1d, factor, xmin, xi, n);
	for(i=0;i<n; i++) xmin[i]=xmin[i]+xmin1d*xi[i];
	
	return test;
}

/*--------------------------------------------------------------------*/

double IDbrentmethod1Dfunc(double (*func)(double[], double[], struct array, struct propagation_parameters*, struct fermi*), double xtra[], struct array spect, struct propagation_parameters* pparam, struct fermi* fe, double x, double xmini[], double xi[], int n)
{
	int i;
	double xitemp[n];
	
	for(i=0;i<n;i++) xitemp[i]=xmini[i]+x*xi[i];
	
	return func(xitemp, xtra,spect,pparam,fe);
}

int IDbrentmethod1D(double (*func)(double[], double[], struct array, struct propagation_parameters*, struct fermi*), double xtra[], struct array spect, struct propagation_parameters* pparam, struct fermi* fe, double a0, double x0,double b0, double *fmin, double *xmin, double factor, double xmini[], double xi[], int n)
{	
	int BrentIterMax=1000;
	double ZEPS=1.0e-70;
	double CGOLD=0.1;
	double tol=1.0e-3;

	double tol1, tol2;
	double a,b;
	
	if(a0<b0)
	{
		a=a0;
		b=b0;
	}
	else
	{
		a=b0;
		b=a0;
	}
	
	double x,w,v; // points used in the parabolic fit
	double u;
	double  pa1, pa2; // parabola parameters pa2*X^2+pa1*X+pa0
	double pamin; // parabola minimum 
	double xm;
	double fx,fw,fv, fu;
	double e=0, etmp,d;
	int iter=0;
	
	x=w=v=x0;
	fx=fw=fv=factor*IDbrentmethod1Dfunc(func,xtra,spect,pparam,fe,x0,xmini,xi,n);
	xm=0.5*(a+b);
	
	do
	{
		tol1=tol*fabs(x)+ZEPS ;
		tol2=e*tol1;

		if (fabs(e)>tol1)
		{
			pa2=((fw-fx)/(w-x)-(fv-fx)/(v-x))/(w-v); // parabolic fit
			pa1=(fv-fx)/(v-x) -pa2*(v+x);
			pamin=-0.5*pa1/pa2;
			etmp=e;
			e=d;
			if(fabs(x-pamin)>0.5*etmp||pamin<a||pamin>b)
			{// cannot use the parabolic fit
				if(x>=xm) e=a-x;
				else e=b-x;
				d=CGOLD*e;
			}
			else
			{
				d=pamin-x;
				u=pamin;
				if((u-a)<tol2||(b-u)<tol2)
				{// pamin too close to borders
					d=fabs(tol1);
					if(xm<x) d=-d;
				}
			}
		}
		else
		{
			if(x>=xm) e=a-x;
			else e=b-x;
			d=CGOLD*e;
		}	
		
		if(fabs(d)>=tol1) u=x+d;
		else
		{
			if(d<0) u=x-fabs(tol1);
			else u=x+fabs(tol1);
		}
		fu=factor*IDbrentmethod1Dfunc(func,xtra,spect,pparam,fe,u,xmini,xi,n);

		if(fu<=fx)
		{
			if(u>=x) a=x;
			else b=x;
			v=w;w=x;x=u;
			fv=fw;fw=fx;fx=fu;
		}
		else
		{
			if(u<x) a=u;
			else b=u;
			if(fu<=fw||w==x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if(fu<=fv||v==x||v==w)
			{
				v=u;
				fv=fu;
			}
		}

		xm=0.5*(a+b);

		iter++;
	}
	while(iter<BrentIterMax&&fabs(x-xm)+0.5*(b-a)> 2*tol1);
	
	*xmin=x;
	*fmin=fx;

	if(iter==BrentIterMax)
	{
		printf("WARNING: too many iterations in Brent's method, precision decreased\n");
		return 0;
	}
	else return 1;
}

/*----------------------------------------------------------------*/

int IDpowellaux(int n, double (*func)(double[], double[], struct array, struct propagation_parameters*, struct fermi*), double xtra[], struct array spect, struct propagation_parameters* pparam, struct fermi* fe, double x0[], double xlimmin[], double xlimmax[], double* fmin,  double xmin[], double ftol, double factor)
{	
	double TINY=1.0e-70;
	int ITMAX=500;
	int nmax=5;

	int i,j,iter;
	double vect[nmax][nmax];
	//initiate set of vectors
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j) vect[i][j]=1;
			else vect[i][j]=0;
		}
	}	

	for(j=0;j<n;j++) xmin[j]=x0[j]; // initiate minimum position

	double vecti[nmax];// current vector in the iteration
	double fmintemp;
	int imax;
	double fdeltamax;
	*fmin=factor*func(x0,xtra,spect,pparam,fe);
	double fiter, xiter[nmax];
	int testbrent; 
	
	double meanvect[nmax], xextrapol[nmax];

	double fextrapol;
	double t;
	
	for(iter=0;iter<=ITMAX;iter++)
	{
			
		fiter=*fmin;
		for(j=0;j<n;j++) xiter[j]=xmin[j];
		fdeltamax=0;
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++) vecti[j]=vect[j][i];
			fmintemp=*fmin;
			testbrent=IDbrentmethod(n,func,xtra,spect,pparam,fe,xmin,vecti,xlimmin,xlimmax,fmin,factor);
			if(testbrent==0)
			{
				if(debug) printf("Brent method unsuccessful\n");
				return 0;
			}

			if((fmintemp-*fmin)>fdeltamax)
			{
				imax=i;
				fdeltamax=fmintemp-*fmin;	
			}
		}


		if(2.*(fiter-(*fmin)) <= ftol*(fabs(fiter)+fabs(*fmin))+TINY) return 0;
		
		if (iter == ITMAX)
		{
			if(debug) printf("Powell exceeding maximum iterations.\n");
			return 0;
		}
		
		for (j=0;j<n;j++)
		{ //check if it is worth to keep the mean deplacement vector
			meanvect[j]=xmin[j]-xiter[j]; 
			xextrapol[j]=xmin[j]+meanvect[j];
		}

		fextrapol=factor*func(xextrapol,xtra,spect,pparam,fe);
		if (fextrapol < fiter)
		{
			t=2.*(fiter-2.*(*fmin)+fextrapol)*pow(fiter-(*fmin)-fdeltamax,2.)-fdeltamax*pow(fiter-fextrapol,2.); 

			if (t < 0.)
			{
				IDbrentmethod(n,func,xtra,spect,pparam,fe,xmin,meanvect, xlimmin, xlimmax, fmin, factor);
				vect[j][imax]=vect[j][n]; 
				vect[j][n]=meanvect[j];
			}
		}

	}
	
	return 1;
}


/*--------------------------------------------------------------------*/

int IDpowell(int n, double (*func)(double[], double[], struct array, struct propagation_parameters*, struct fermi*), double xtra[], struct array spect, struct propagation_parameters* pparam, struct fermi* fe, double xlimmin[], double xlimmax[], double* fmin, double xmin[], double ftol, char* option)
{
	int nmax=5;
	
	srand((unsigned int)time(NULL));
	double x0[nmax];
	int i;

	double factor=0.;
	if(strcmp(option,"min")==0) factor=1.;
	else if (strcmp(option,"max")==0) factor=-1.;

	for(i=0;i<n;i++) x0[i]=(xlimmin[i]+xlimmax[i])/2.;
	int test=IDpowellaux(n,func,xtra,spect,pparam,fe,x0,xlimmin, xlimmax, fmin,  xmin, ftol, factor);
	if(!strcmp(option, "max")) *fmin=-*fmin; 
	
	return test;
}
