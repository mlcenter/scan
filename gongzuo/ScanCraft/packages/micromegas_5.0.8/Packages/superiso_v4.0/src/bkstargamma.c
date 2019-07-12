#include "include.h"


/*----------------------------------------------------------------------*/

double BR_BKstargamma(int charge, double complex C0b[], double complex C1b[], double complex C2b[], double complex Cpb[], struct parameters* param, double mu_b)
{	
	double mB,mKstar,eq;
	if(charge==0) 
	{
		mB=param->m_Bd;
		mKstar=param->m_Kstar0;
		eq=-1./3.;
	}
	else
	{
		mB=param->m_B;
		mKstar=param->m_Kstar;
		eq=2./3.;
	}

	double mc=mc_pole_1loop(param);
	double mbpole=mb_pole_1loop(param);
	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);

	int ie,je;
	
	double alpha_em=1./133.;
// 	double alpha_em=1./137.; /*SN: second line after Eq. 59 in 0106067*/

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double complex Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=mB/2.;

	int nf=5;
	double f_Kstar_perp=param->f_Kstar_perp;
	f_Kstar_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_Kstar_par=param->f_Kstar_par;
	
	double T1_0=param->T1_BKstar;
	double xi_perp=T1_0;

	double complex C7eff=Cmub[7];
	double complex C8eff=Cmub[8];
				
	double complex C7effp=Cpb[7];

	double complex C1bar=Cmub[1]/2.;
	double complex C2bar=Cmub[2]-Cmub[1]/6.;
	double complex C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
	double complex C4bar=Cmub[4]/2.+8.*Cmub[6];
	double complex C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
	double complex C6bar=Cmub[4]/2.+2.*Cmub[6];

	double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
	double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */

	double complex Cperpp0=0.;
	double complex Cperpm0=0.;

	double complex Cperpp0u=0.;
	double complex Cperpm0u=0.;

	double logb=log(mu_b/mb);

	double epsilon=1.e-10;
	double shat=epsilon; 

	double mchat=mc/mb; 
	double z=mchat*mchat;	

	double complex Cperppf=0.;
	double complex Cperpmf=0.;

	double complex F27=F27_bsll(shat,z,logb);
	double complex F87=F87_bsll(shat,logb);
	double complex F27_u=F27u(shat,logb);
	
	double complex Cperpnf=(-C2bar*F27-C8eff*F87)/4.*3.;
	
	double complex Cperpnfu=(-C2bar*(F27+F27_u))/4.*3.;
	
	double complex Cperpp1=Cperppf+Cperpnf;
	double complex Cperpm1=Cperpmf+Cperpnf;

	double complex Cperpp=Cperpp0+alphas_mub*4./3./4./pi*Cperpp1;
	double complex Cperpm=Cperpm0+alphas_mub*4./3./4./pi*Cperpm1;
			
	double complex Cperppu=Cperpp0u+alphas_mub*4./3./4./pi*Cperpnfu; 
	double complex Cperpmu=Cperpm0u+alphas_mub*4./3./4./pi*Cperpnfu; 
			
	double Xi_perp=1.;

	double eu=2./3.;
	double ed=-1./3.;

	double a1perp=param->a1perp;
	double a2perp=param->a2perp;
	double a1par=param->a1par;
	double a2par=param->a2par;

	a1perp*=pow(eta,4./3.*(4.*1./2.)/(11.-2./3.*nf));
	a2perp*=pow(eta,4./3.*(4.*(1./2.+1./3.))/(11.-2./3.*nf));
	a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
	a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

	double u;

	double complex int_perppp,int_perppm,int_perpmp,int_perpmm;
	double complex int_perppu;
	double complex Tperppp0,Tperpppf,Tperpppnf,Tperppp;
	double complex Tperppm0,Tperppmf,Tperppmnf,Tperppm;
	double complex Tperpmp0,Tperpmpf,Tperpmpnf,Tperpmp;
	double complex Tperpmm0,Tperpmmf,Tperpmmnf,Tperpmm;
	double complex Tperppnfu,Tperppu;
	
	int_perppp=int_perppm=int_perpmp=int_perpmm=0.;
	int_perppu=0.;

	double lambda_Bp=param->lambda_Bp;
	lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

	double phiKstar_perp;
	double complex tperp_mb,tperp_mc,tperp_0;

	double complex integ3=0;

	double complex Fperp=0.;
	double complex Xperp=0.;
	double x;
	double complex integ4=0.;
	double complex FV;
	double complex FVu;
	double complex integ5=0.;

	double zeta3A=param->zeta3A;
	double zeta3V=param->zeta3V;
	double wA10=param->wA10;
	double deltatp=param->deltatp;
	double deltatm=param->deltatm;
		
	int n1=10;
	int n1sav=n1;
	for(ie=0;ie<=n1;ie++)
	{
		u=(double)ie/n1;
		if(ie==0) n1*=2;
		if(ie==n1){u=(double)(ie-1)/n1;n1*=2;}

	/* Tperp */		
		Tperppp0=Tperpmp0=0.;
		Tperpppf=0.;
		Tperpmpf=0.;

		Tperppm0=Tperpmm0=0.;
		Tperppmf=Tperpmmf=0.;
		
		phiKstar_perp=phi_Kstar(u,a1perp,a2perp);
		tperp_mc=tperp_bkll(u,mc,0,E_Kstar,param);
		tperp_mb=tperp_bkll(u,mb,0,E_Kstar,param);
		tperp_0=tperp_bkll(u,0.,0,E_Kstar,param);
		
		Tperpppnf=Tperpmpnf=
		// -4.*ed*C8eff/u /*SN: I have removed this term which would result in an infinite value ---> see first bulletpoint in page 6 of 1608.02556v3*/
		+mB/2./mb*(eu*tperp_mc*(C2bar+C4bar-C6bar)
		+ed*tperp_mb*(C3bar+C4bar-C6bar-4.*mb/mB*C5bar)
		+ed*tperp_0*C3bar);
			
		Tperppmnf=Tperpmmnf=0.;

		Tperppnfu=mB/2./mb*eu*(tperp_mc-tperp_0)*(Cmub[2]-Cmub[1]/6.);
	
		Tperppp=Tperppp0+alphas_muf*4./3./4./pi*(Tperpppf+Tperpppnf); 
		Tperppm=Tperppm0+alphas_muf*4./3./4./pi*(Tperppmf+Tperppmnf);
		Tperpmp=Tperpmp0+alphas_muf*4./3./4./pi*(Tperpmpf+Tperpmpnf);
		Tperpmm=Tperpmm0+alphas_muf*4./3./4./pi*(Tperpmmf+Tperpmmnf);

		Tperppu=alphas_muf*4./3./4./pi*Tperppnfu;

		int_perppp+=phiKstar_perp*Tperppp/n1/lambda_Bp; 
		int_perpmp+=phiKstar_perp*Tperpmp/n1/lambda_Bp;

		int_perppu+=phiKstar_perp*Tperppu/n1/lambda_Bp;
		
		
		integ3+=phiKstar_perp/(1.-u)/n1;

		x=(1.-u)*mB*mB;

		double h_mc=h_bkll(x,mc,mu_b);
		double h_mb=h_bkll(x,mbpole,mu_b);
		double h_0=h_bkll(x,0.,mu_b);

		FV=3./4.*(h_mc*(C2bar+C4bar+C6bar)+h_mb*(C3bar+C4bar+C6bar)+h_0*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));

		integ4+=phiKstar_perp/(1.-u)*FV/n1;

		Fperp+=phiKstar_perp/(1.-u)/3./n1;
		Xperp+=(u<=1.-0.5/mB)*phiKstar_perp/pow((1.-u),2.)/3.*(0.5/mB)/n1;

		integ5+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;

		if(ie==0||ie==n1sav) n1=n1sav;
	}
	
	/* Tau_perp */		

	double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perppp+int_perppm); 
	
	double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_B*f_Kstar_perp/mB*Xi_perp*(int_perpmp+int_perpmm); 
						
	double complex DeltaTauperpWA=-eq*4.*pi*pi/3.*param->f_B*f_Kstar_perp/mb/mB*(Cmub[3]+4./3.*Cmub[4]+4.*Cmub[5]+16./3.*Cmub[6])*integ3
	+eq*2.*pi*pi/3.*param->f_B*f_Kstar_par/mb/mB*mKstar/lambda_Bp*(Cmub[3]+4./3.*Cmub[4]+16.*Cmub[5]+64./3.*Cmub[6]);

	double rho=0.;
	double phi=0.;
	Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;

	double complex DeltaTauperpHSA=eq*4./3.*alphas_muf/4./pi*pi*pi*param->f_B/3./mb/mB*(12.*C8eff*mb/mB*f_Kstar_perp*Xperp
	+8.*f_Kstar_perp*integ4-4.*mKstar*f_Kstar_par/lambda_Bp*integ5);

	Tauperpp+= DeltaTauperpWA+DeltaTauperpHSA;
	Tauperpm+= DeltaTauperpWA+DeltaTauperpHSA;			
				
	double complex hplus0,hminus0;
	hplus0=hminus0=0.;

	double Tauperpmu,Tauperppu,Tauperpmu_bar,Tauperppu_bar;
	Tauperpmu=Tauperppu=Tauperpmu_bar=Tauperppu_bar=0.;
	
	double complex Nprimecoeff=(-4.*param->Gfermi)/sqrt(2.)*mB*alpha_em/4./pi*param->Vtb*conj(param->Vts);
	double complex Nprimecoeff_bar=(-4.*param->Gfermi)/sqrt(2.)*mB*alpha_em/4./pi*conj(param->Vtb)*param->Vts;
	
	double complex HVplus,HVminus,HVplus_bar,HVminus_bar;
	HVplus=HVminus=HVplus_bar=HVminus_bar=0.;
	
	/*SN: 1) When considering only the non-factorisable piece Tauperpm = Tauperpp, so there is not need to have them separately (I have written it this way to be similar to when using the soft FF method) **/
	/*SN: 2)currently the Cabbibo suppressed contributions Tauperpmu and Tauperppu have been neglected but I have kept them in the general formula in case we want to add it later  **/
	HVplus=I*Nprimecoeff*mB*mB/(2.*sqrt(alpha_em*pi))*( 2.*param->mass_b/mB*(mB*mB-mKstar*mKstar)/mB/mB*(-C7effp*T1_0 + 1./2.*((Tauperpm+Tauperpmu)-(Tauperpp+Tauperppu)) ) + (-16.*pi*pi*hplus0) );
	HVminus=I*Nprimecoeff*mB*mB/(2.*sqrt(alpha_em*pi))*( 2.*param->mass_b/mB*(mB*mB-mKstar*mKstar)/mB/mB*(C7eff*T1_0 + 1./2.*((Tauperpm+Tauperpmu)+(Tauperpp+Tauperppu)) )+(-16.*pi*pi*hminus0));

	HVplus_bar=I*Nprimecoeff_bar*mB*mB/(2.*sqrt(alpha_em*pi))*( 2.*param->mass_b/mB*(mB*mB-mKstar*mKstar)/mB/mB*(-C7effp*T1_0 + 1./2.*((Tauperpm+Tauperpmu_bar)-(Tauperpp+Tauperppu_bar)) ) + (-16.*pi*pi*hplus0) );
	HVminus_bar=I*Nprimecoeff_bar*mB*mB/(2.*sqrt(alpha_em*pi))*( 2.*param->mass_b/mB*(mB*mB-mKstar*mKstar)/mB/mB*(C7eff*T1_0 + 1./2.*((Tauperpm+Tauperpmu_bar)+(Tauperpp+Tauperppu_bar)) )+(-16.*pi*pi*hminus0));
	
	double Gamma_BKstargamma=(mB*mB-mKstar*mKstar)/16./pi/mB/mB/mB*(cabs(HVplus*conj(HVplus))+cabs(HVminus*conj(HVminus)));
	double Gamma_BKstargamma_bar=(mB*mB-mKstar*mKstar)/16./pi/mB/mB/mB*(cabs(HVplus_bar*conj(HVplus_bar))+cabs(HVminus_bar*conj(HVminus_bar)));
				
	if(charge==0) return param->life_Bd/hbar*(Gamma_BKstargamma+Gamma_BKstargamma_bar)/2.;
	else return param->life_B/hbar*(Gamma_BKstargamma+Gamma_BKstargamma_bar)/2.;
}

/*----------------------------------------------------------------------*/

double BR_BKstargamma_calculator(int charge, char name[])
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQpb[3];
	
	struct parameters param;
				
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BR_BKstargamma(charge,C0b,C1b,C2b,Cpb,&param,mu_b);
}
