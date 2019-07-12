#include "include.h"

/*----------------------------------------------------------------------*/
/*------------------------------ SOFT ----------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_BKll_dq2_soft(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double mB,mK,eq;
	if(charge==0) 
	{
		mB=param->m_Bd;
		mK=param->m_K0;
		eq=-1./3.;
	}
	else
	{
		mB=param->m_B;
		mK=param->m_K;
		eq=2./3.;
	}
	
	double mc=mc_pole_1loop(param);
	double mbpole=mb_pole_1loop(param);
	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
	
	double beta_l=sqrt(1.-4.*ml*ml/q2);

	double alpha_em=1./133.;

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	int ie,je;
	double complex Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_K=(mB*mB+mK*mK-q2)/2./mB;

	int nf=5;
	
	
	/********LCSR+Lattice fit parameters from Altmannshofer and Straub 1411.3161***************/
	double tau_plus=pow(mB+mK,2.);
	double tau_0=(mB+mK)*pow(sqrt(mB)-sqrt(mK),2.);
	double z_q2=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));

	double P_p=1.-q2/pow(mB+param->DmBp_BK,2.);

	double P_T=1.-q2/pow(mB+param->DmBT_BK,2.);

	double f0=param->a00_BK+param->a10_BK*z_q2+param->a20_BK*pow(z_q2,2.)+param->a30_BK*pow(z_q2,3.);
	double fp=(1./P_p)*(param->a0p_BK+param->a1p_BK*z_q2+param->a2p_BK*pow(z_q2,2.)+(-param->a1p_BK/3.+2.*param->a2p_BK/3.)*pow(z_q2,3.));
	double fT=(1./P_T)*(param->a0T_BK+param->a1T_BK*z_q2+param->a2T_BK*pow(z_q2,2.)+(-param->a1T_BK/3.+2.*param->a2T_BK/3.)*pow(z_q2,3.));

	
	double CT=0.;
	//double CT5=0.;
	
	double complex FV=0.;
	double complex FA=0.;
	double complex FS=0.;
	double complex FP=0.;
	double complex FT=0.;
	double complex FT5=0.;

	double a_l=0.;
	double b_l=0.;
	double c_l=0.;

	if(q2<14.)
	{	
		double complex C7eff=Cmub[7];
		double complex C8eff=Cmub[8];
		double complex C9=Cmub[9];
		double complex C10=Cmub[10];
	
		double complex C7effp=Cpb[7];
		double complex C9p=Cpb[9];
		double complex C10p=Cpb[10];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double xi_P=fp;
	
		double fT_fp=0.;
		double f0_fp=0.;
	
	
		double complex C1bar=Cmub[1]/2.;
		double complex C2bar=Cmub[2]-Cmub[1]/6.;
		double complex C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
		double complex C4bar=Cmub[4]/2.+8.*Cmub[6];
		double complex C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
		double complex C6bar=Cmub[4]/2.+2.*Cmub[6];
	
		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */

		double complex h_mc=h_bkll(q2,mc,mu_b);
		double complex h_mb=h_bkll(q2,mbpole,mu_b);
		double complex h_0=h_bkll(q2,0.,mu_b);

		double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
		+h_mc*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
		+h_mb*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
		+h_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);
		


		double complex Cparp0=-(C7eff+C7effp)-mB/2./mb*Y;

		double logb=log(mu_b/mb);
		double DeltaM=-6.*logb-4.*(1.-mu_f/mb);
		double L=-(mb*mb-q2)/q2*log(1.-q2/mb/mb);
	
		double shat=q2/mb/mb;
		double mchat=mc/mb;
		double z=mchat*mchat;	

		double complex Cparpf=-(C7eff+C7effp)*(-2.*logb+2.*L+DeltaM);

		double complex F27=F27_bsll(shat,z,logb);
		double complex F87=F87_bsll(shat,logb);
		double complex F29=F29_bkll(shat,z,logb);
		double complex F19=F19_bkll(shat,z,logb);
		double complex F89=F89_bsll(shat);

		double complex Cparnf=(C2bar*F27+C8eff*F87+mB/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		
		double complex Cparp1=Cparpf+Cparnf;
	
		double complex Cparp=Cparp0+alphas_mub*4./3./4./pi*Cparp1;
		
		double eu=2./3.;
		double ed=-1./3.;

		double a1K=param->a1K;
		double a2K=param->a2K;

		a1K*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2K*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_parpp,int_parpm;
		double complex Tparpp0,Tparppf,Tparppnf,Tparpp;
		double complex Tparpm0,Tparpmf,Tparpmnf,Tparpm;
		int_parpp=int_parpm=0.;

		double int_DeltaFK=0.;

		double lambda_Bp=param->lambda_Bp;
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(mB-mb)/3.;
		double complex lambda_Bm=1./(exp(-q2/mB/omega0)/omega0*(-Ei(q2/mB/omega0)+I*pi));

		double phiK;
		double complex tpar_mb,tpar_mc,tpar_0;


		int n1=10;
		int n1sav=n1; 
		for(ie=0;ie<=n1;ie++) 
		{
			u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=0.99;n1*=2;}
			

		/* Tpar */		

			phiK=phi_Kstar(u,a1K,a2K);
			tpar_mc=tpar_bkll(u,mc,q2,E_K,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_K,param);
			tpar_0=tpar_bkll(u,0.,q2,E_K,param);
		
			Tparpp0=0.; 
		
			Tparppf=(C7eff+C7effp)*4.*mB/E_K/(1.-u);

			Tparppnf=mB/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
					
			Tparpm0=-eq*4.*mB/mb*(C3bar+3.*C4bar);
					
			Tparpmf=0.;

			h_mc=h_bkll((1.-u)*mB*mB+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*mB*mB+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*mB*mB+u*q2,0.,mu_b);

			Tparpmnf=eq*(8.*C8eff/((1.-u)+u*q2/mB/mB)
			+6.*mB/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf); 
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);

			int_parpp+=(phiK*Tparpp/lambda_Bp)/n1;
			int_parpm+=(phiK*Tparpm/lambda_Bm)/n1;
			
			int_DeltaFK+=(phiK/lambda_Bp/(1.-u))/n1;

			if(ie==0||ie==n1sav) n1=n1sav; 
		}
		
		
		/* Tau_P */		
		double complex TauP=xi_P*(-Cparp)+pi*pi/3.*param->f_B*param->f_K/mB*(-int_parpp-int_parpm);
		
		double DeltaFK=8.*pi*pi/3.*param->f_B*param->f_K/mB*int_DeltaFK;
		
		f0_fp=2.*E_K/mB*(1.+alphas_muf*4./3./2./pi*(1.-L)+alphas_muf*4./3./16./pi*mB*(mB-2.*E_K)/E_K/E_K*DeltaFK/xi_P);
		
		fT_fp=(mK+mB)/mB*(1.+alphas_muf*4./3./2./pi*(-logb+L)-alphas_muf*4./3./8./pi*mB/E_K*DeltaFK/xi_P);

		
		double lambda=pow(mB,4.)+pow(mK,4.)+q2*q2-2.*(mB*mB*mK*mK+mB*mB*q2+mK*mK*q2);
	
		FV=xi_P*((C9+C9p)+2./(mB+mK)*(mb*TauP/xi_P+4.*CT*ml*fT_fp));
		FA=xi_P*(C10+C10p);
		FS=xi_P*(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*f0_fp*(CQ1+CQ1p);
		FP=xi_P*((mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*f0_fp*(CQ2+CQ2p)-ml*(C10+C10p)*(1.-(mB*mB-mK*mK)/q2*(f0_fp-1.)));

		//FT=xi_P*2.*sqrt(lambda)*beta_l/(mB+mK)*CT*fT_fp;
		//FT5=xi_P*2.*sqrt(lambda)*beta_l/(mB+mK)*CT5*fT_fp;

		/* hadronic uncertainties */		
		FV*=1.+param->BtoKlow_FV_err_noq2+q2/6.*param->BtoKlow_FV_err_q2;
		FA*=1.+param->BtoKlow_FA_err_noq2+q2/6.*param->BtoKlow_FA_err_q2;
		FS*=1.+param->BtoKlow_FS_err_noq2+q2/6.*param->BtoKlow_FS_err_q2;
		FP*=1.+param->BtoKlow_FP_err_noq2+q2/6.*param->BtoKlow_FP_err_q2;

		double Gamma_0=pow(cabs(param->Vtb*conj(param->Vts)),2.)*param->Gfermi*param->Gfermi*alpha_em*alpha_em/512./pow(pi,5.)/pow(mB,3.);

		double C_K=Gamma_0*sqrt(lambda)*beta_l;
		
		a_l=C_K*(q2*(beta_l*beta_l*cabs(FS*conj(FS))+cabs(FP*conj(FP)))+lambda/4.*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*(mB*mB-mK*mK+q2)*creal(FP*conj(FA))+4.*ml*ml*mB*mB*cabs(FA*conj(FA))); 
		
		b_l=2.*C_K*(q2*(beta_l*beta_l*creal(FS*conj(FT))+creal(FP*conj(FT5)))+ml*(sqrt(lambda)*beta_l*creal(FS*conj(FV))+(mB*mB-mK*mK+q2)*creal(FT5*conj(FA)))); 
		
		c_l=C_K*(q2*(beta_l*beta_l*cabs(FT*conj(FT))+cabs(FT5*conj(FT5)))-lambda/4.*beta_l*beta_l*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*sqrt(lambda)*beta_l*creal(FT*conj(FV))); 
	}
	

	if(q2>6.)
	{
		double complex C10=Cmub[10];
	
		double complex C7effp=Cpb[7];
		double complex C9p=Cpb[9];
		double complex C10p=Cpb[10];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];
	
		
		double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
		double shat=q2/mB/mB;
		double z=4.*mb*mb/q2;
		
		double complex x1=0.5+0.5*I*csqrt(z-1.);
		double complex x2=0.5-0.5*I*csqrt(z-1.);
		double complex x3=0.5+0.5*I/csqrt(z-1.);
		double complex x4=0.5-0.5*I/csqrt(z-1.);
		
		double complex A=
		-104./243.*log(mb*mb/mu_b/mu_b)+4.*shat/27./(1.-shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*log(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(csqrt(z-1.)))+9.*shat*shat*shat*log(shat)*log(shat)+18.*pi*I*shat*(1.-2.*shat)*log(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
		
		double complex B=
		8./243./shat*((4.-34.*shat-17.*pi*I*shat)*log(mb*mb/mu_b/mu_b)+8.*shat*pow(log(mb*mb/mu_b/mu_b),2.)+17.*shat*log(shat)*log(mb*mb/mu_b/mu_b))
		+(2.+shat)*csqrt(z-1.)/729./shat*(-48.*log(mb*mb/mu_b/mu_b)*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
		-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
		-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
		-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
		+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
		-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
		+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
		double complex C=-16./81.*log(q2/mu_b/mu_b)+428./243.-64./27.*zeta3+16./81.*pi*I;
		
		double complex C9eff=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(shat))
		+8.*mc*mc/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts))+2.*Cmub[3]+20.*Cmub[5]);
		
		double complex C7eff=Cmub[7]
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(shat,log(mu_b/mb)));
		
		//double kappa=1.-2.*alphas_mub/3./pi*log(mu_b/mb);
		
		
		FA=(C10+C10p)*fp;
		FV=(C9eff+C9p)*fp+2.*mb/(mB+mK)*(C7eff+C7effp+4.*ml/mb*CT)*fT;
		FS=(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*(CQ1+CQ1p)*f0;
		FP=(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*(CQ2+CQ2p)*f0-ml*(C10+C10p)*(fp-(mB*mB-mK*mK)/q2*(f0-fp));

		//FT=2.*sqrt(lambda)*beta_l/(mB+mK)*CT*fT;
		//FT5=2.*sqrt(lambda)*beta_l/(mB+mK)*CT5*fT;

		/* hadronic uncertainties */		
		FV*=1.+param->BtoKhigh_FV_err;
		FA*=1.+param->BtoKhigh_FA_err;
		FS*=1.+param->BtoKhigh_FS_err;
		FP*=1.+param->BtoKhigh_FP_err;
	
		double lambda=pow(mB,4.)+pow(mK,4.)+q2*q2-2.*(mB*mB*mK*mK+mB*mB*q2+mK*mK*q2);
	
		double Gamma_0=pow(cabs(param->Vtb*conj(param->Vts)),2.)*param->Gfermi*param->Gfermi*alpha_em*alpha_em/512./pow(pi,5.)/pow(mB,3.);
	
		double C_K=Gamma_0*sqrt(lambda)*beta_l;
	
	
		double a_l_high=C_K*(q2*(beta_l*beta_l*cabs(FS*conj(FS))+cabs(FP*conj(FP)))+lambda/4.*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*(mB*mB-mK*mK+q2)*creal(FP*conj(FA))+4.*ml*ml*mB*mB*cabs(FA*conj(FA))); 
		
		double b_l_high=2.*C_K*(q2*(beta_l*beta_l*creal(FS*conj(FT))+creal(FP*conj(FT5)))+ml*(sqrt(lambda)*beta_l*creal(FS*conj(FV))+(mB*mB-mK*mK+q2)*creal(FT5*conj(FA)))); 
		
		double c_l_high=C_K*(q2*(beta_l*beta_l*cabs(FT*conj(FT))+cabs(FT5*conj(FT5)))-lambda/4.*beta_l*beta_l*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*sqrt(lambda)*beta_l*creal(FT*conj(FV))); 
	
		if(q2>14.)
		{
			a_l=a_l_high;
			b_l=b_l_high;
			c_l=c_l_high;
		}
		else
		{
			a_l=a_l*(14.-q2)/8.+a_l_high*(q2-6.)/8.;
			b_l=b_l*(14.-q2)/8.+b_l_high*(q2-6.)/8.;
			c_l=c_l*(14.-q2)/8.+c_l_high*(q2-6.)/8.;
		}

	}


	double dGamma_BKll_dq2=2.*(a_l+c_l/3.);

	double AFB[3],FH[3];

	AFB[0]=b_l/dGamma_BKll_dq2;
	AFB[1]=b_l;
	AFB[2]=dGamma_BKll_dq2;

	FH[0]=2.*(a_l+c_l)/dGamma_BKll_dq2;
	FH[1]=2.*(a_l+c_l);
	FH[2]=dGamma_BKll_dq2;

	for(je=0;je<=Nobs_BKll;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=FH[ie];
	}

	return dGamma_BKll_dq2;
}

/*----------------------------------------------------------------------*/

double dGamma_BKll_dq2_soft_calculator(int gen, int charge, double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(B->K mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(gen,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(gen,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(gen,Cpb,CQpb,mu_W,mu_b,&param);

	return dGamma_BKll_dq2_soft(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_soft(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_BKll+1],obs_den[Nobs_BKll+1];
	for(je=0;je<=Nobs_BKll;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated FH */
	
	double dobs[Nobs_BKll+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=s0m=0.;
	s0p=1.;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGamma_BKll_dq2_soft(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	
	for(ie=1;ie<nmax;ie++)	
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s+=h;
		s0p=s;
		
		Gamma+=4.*dGamma_BKll_dq2_soft(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKll;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		Gamma+=2.*dGamma_BKll_dq2_soft(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKll;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}

		if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	
	s=smax;
	Gamma+=4.*dGamma_BKll_dq2_soft(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGamma_BKll_dq2_soft(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	
	obs[0]=s0;
	for(je=1;je<=Nobs_BKll;je++) obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_BKll;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	if(charge==0) return param->life_Bd/hbar*Gamma;
	else return param->life_B/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double RK_BKll_soft(int charge, double smin, double smax, double complex C0be[], double complex C1be[], double complex C2be[], double complex CQ0be[], double complex CQ1be[], double complex Cpbe[], double complex CQpbe[], double complex C0bmu[], double complex C1bmu[], double complex C2bmu[], double complex CQ0bmu[], double complex CQ1bmu[], double complex Cpbmu[], double complex CQpbmu[], struct parameters* param, double mu_b)
{
	double obs[Nobs_BKll+1];
	
	return BRBKll_soft(2,charge,smin,smax,obs,C0bmu,C1bmu,C2bmu,CQ0bmu,CQ1bmu,Cpbmu,CQpbmu,param,mu_b)/BRBKll_soft(1,charge,smin,smax,obs,C0be,C1be,C2be,CQ0be,CQ1be,Cpbe,CQpbe,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_lowq2_soft(int gen,int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKll_soft(gen,charge,1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_highq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKll_soft(gen,charge,14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKmumu_lowq2_soft_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->K mu+ mu-) and all the other observables */

	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11], CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKll_lowq2_soft(2,1,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKmumu_highq2_soft_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->K mu+ mu-) and all the other observables */
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKll_highq2_soft(2,1,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double RK_BKll_soft_calculator(int charge, double smin, double smax, char name[])
{
	double complex C0bmu[11],C1bmu[11],C2bmu[11],C0wmu[11],C1wmu[11],C2wmu[11],Cpbmu[11],CQ0bmu[3],CQ1bmu[3],CQpbmu[3];
	double complex C0be[11],C1be[11],C2be[11],C0we[11],C1we[11],C2we[11],Cpbe[11],CQ0be[3],CQ1be[3],CQpbe[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(2,C0wmu,C1wmu,C2wmu,mu_W,&param);
	C_calculator_base1(C0wmu,C1wmu,C2wmu,mu_W,C0bmu,C1bmu,C2bmu,mu_b,&param);
	CQ_calculator(2,CQ0bmu,CQ1bmu,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpbmu,CQpbmu,mu_W,mu_b,&param);

	CW_calculator(1,C0we,C1we,C2we,mu_W,&param);
	C_calculator_base1(C0we,C1we,C2we,mu_W,C0be,C1be,C2be,mu_b,&param);
	CQ_calculator(1,CQ0be,CQ1be,mu_W,mu_b,&param);
	Cprime_calculator(1,Cpbe,CQpbe,mu_W,mu_b,&param);
	
	return RK_BKll_soft(charge,smin,smax,C0be,C1be,C2be,CQ0be,CQ1be,Cpbe,CQpbe,C0bmu,C1bmu,C2bmu,CQ0bmu,CQ1bmu,Cpbmu,CQpbmu,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*------------------------------ FULL ----------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_BKll_dq2_full(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double mB,mK,eq;
	if(charge==0) 
	{
		mB=param->m_Bd;
		mK=param->m_K0;
		eq=-1./3.;
	}
	else
	{
		mB=param->m_B;
		mK=param->m_K;
		eq=2./3.;
	}
	
	double mc=mc_pole_1loop(param);
	double mbpole=mb_pole_1loop(param);
	double mb_mub=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
	
	double beta_l=sqrt(1.-4.*ml*ml/q2);

	double alpha_em=1./133.;

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	int ie,je;
	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_K=(mB*mB+mK*mK-q2)/2./mB;

	int nf=5;
	
	
	/********LCSR+Lattice fit parameters from Altmannshofer and Straub 1411.3161***************/
	double tau_plus=pow(mB+mK,2.);
	//double tau_minus=pow(mB-mK,2.);
	double tau_0=(mB+mK)*pow(sqrt(mB)-sqrt(mK),2.);
	double z_q2=(sqrt(tau_plus-q2)-sqrt(tau_plus-tau_0))/(sqrt(tau_plus-q2)+sqrt(tau_plus-tau_0));

	double P_p=1.-q2/pow(mB+param->DmBp_BK,2.);

	double P_T=1.-q2/pow(mB+param->DmBT_BK,2.);

	double f0=param->a00_BK+param->a10_BK*z_q2+param->a20_BK*pow(z_q2,2.)+param->a30_BK*pow(z_q2,3.);
	double fp=(1./P_p)*(param->a0p_BK+param->a1p_BK*z_q2+param->a2p_BK*pow(z_q2,2.)+(-param->a1p_BK/3.+2.*param->a2p_BK/3.)*pow(z_q2,3.));
	double fT=(1./P_T)*(param->a0T_BK+param->a1T_BK*z_q2+param->a2T_BK*pow(z_q2,2.)+(-param->a1T_BK/3.+2.*param->a2T_BK/3.)*pow(z_q2,3.));

	
	double CT=0.;
	//double CT5=0.;
	
	double complex FV=0.;
	double complex FA=0.;
	double complex FS=0.;
	double complex FP=0.;
	double complex FT=0.;
	double complex FT5=0.;

	double a_l=0.;
	double b_l=0.;
	double c_l=0.;

	if(q2<14.)
	{	
		double complex C7eff=Cmub[7];
		double complex C8eff=Cmub[8];
		double complex C9=Cmub[9];
		double complex C10=Cmub[10];
	
		double complex C7effp=Cpb[7];
		double complex C9p=Cpb[9];
		double complex C10p=Cpb[10];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];

		double xi_P=fp;
	
		//double fT_fp=0.;
		//double f0_fp=0.;
	
		double complex C1bar=Cmub[1]/2.;
		double complex C2bar=Cmub[2]-Cmub[1]/6.;
		double complex C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
		double complex C4bar=Cmub[4]/2.+8.*Cmub[6];
		double complex C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
		double complex C6bar=Cmub[4]/2.+2.*Cmub[6];
	
		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */

		double complex h_mc=h_bkll(q2,mc,mu_b);
		double complex h_mb=h_bkll(q2,mbpole,mu_b);
		double complex h_0=h_bkll(q2,0.,mu_b);

		double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
		+h_mc*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
		+h_mb*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
		+h_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);
		


		double complex Cparp0=0.;

		double logb=log(mu_b/mb);
		double DeltaM=-6.*logb-4.*(1.-mu_f/mb);
		//double L=-(mb*mb-q2)/q2*log(1.-q2/mb/mb);
	
		double shat=q2/mb/mb;
		double mchat=mc/mb;
		double z=mchat*mchat;	

		double complex Cparpf=0.;

		double complex F27=F27_bsll(shat,z,logb);
		double complex F87=F87_bsll(shat,logb);
		double complex F29=F29_bkll(shat,z,logb);
		double complex F19=F19_bkll(shat,z,logb);
		double complex F89=F89_bsll(shat);

		double complex Cparnf=(C2bar*F27+C8eff*F87+mB/2./mb*(C2bar*F29+2.*C1bar*(F19+1./6.*F29)+C8eff*F89))/4.*3.;
		
		double complex Cparp1=Cparpf+Cparnf;
	
		double complex Cparp=Cparp0+alphas_mub*4./3./4./pi*Cparp1;
		
		double eu=2./3.;
		double ed=-1./3.;

		double a1K=param->a1K;
		double a2K=param->a2K;

		a1K*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2K*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_parpp,int_parpm;
		double complex Tparpp0,Tparppf,Tparppnf,Tparpp;
		double complex Tparpm0,Tparpmf,Tparpmnf,Tparpm;
		int_parpp=int_parpm=0.;

		double int_DeltaFK=0.;

		double lambda_Bp=param->lambda_Bp;
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_f,2.))*(1.-2.*1.4);

		double omega0=2.*(mB-mb)/3.;
		double complex lambda_Bm=1./(exp(-q2/mB/omega0)/omega0*(-Ei(q2/mB/omega0)+I*pi));

		double phiK;
		double complex tpar_mb,tpar_mc,tpar_0;


		int n1=10;
		int n1sav=n1; 
		for(ie=0;ie<=n1;ie++) 
		{
			u=(double)ie/n1;
			if(ie==0) n1*=2;
			if(ie==n1){u=0.99;n1*=2;}
			

		/* Tpar */		

			phiK=phi_Kstar(u,a1K,a2K);
			tpar_mc=tpar_bkll(u,mc,q2,E_K,param);
			tpar_mb=tpar_bkll(u,mb,q2,E_K,param);
			tpar_0=tpar_bkll(u,0.,q2,E_K,param);
		
			Tparpp0=0.; 
		
			Tparppf=0.;

			Tparppnf=mB/mb*(eu*tpar_mc*(C2bar+C4bar-C6bar)
			+ed*tpar_mb*(C3bar+C4bar-C6bar)
			+ed*tpar_0*C3bar);
					
			Tparpm0=-eq*4.*mB/mb*(C3bar+3.*C4bar);
					
			Tparpmf=0.;

			h_mc=h_bkll((1.-u)*mB*mB+u*q2,mc,mu_b);
			h_mb=h_bkll((1.-u)*mB*mB+u*q2,mbpole,mu_b);
			h_0=h_bkll((1.-u)*mB*mB+u*q2,0.,mu_b);

			Tparpmnf=eq*(8.*C8eff/((1.-u)+u*q2/mB/mB)
			+6.*mB/mb*(h_mc*(C2bar+C4bar+C6bar)
			+h_mb*(C3bar+C4bar+C6bar)
			+h_0*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			Tparpp=Tparpp0+alphas_muf*4./3./4./pi*(Tparppf+Tparppnf); 
			Tparpm=Tparpm0+alphas_muf*4./3./4./pi*(Tparpmf+Tparpmnf);

			int_parpp+=(phiK*Tparpp/lambda_Bp)/n1;
			int_parpm+=(phiK*Tparpm/lambda_Bm)/n1;
			
			int_DeltaFK+=(phiK/lambda_Bp/(1.-u))/n1;

			if(ie==0||ie==n1sav) n1=n1sav; 
		}
		
		
		/* Tau_P */		
		double complex TauP=xi_P*(-Cparp)+pi*pi/3.*param->f_B*param->f_K/mB*(-int_parpp-int_parpm);
		
		//double DeltaFK=8.*pi*pi/3.*param->f_B*param->f_K/mB*int_DeltaFK;
		
		//f0_fp=2.*E_K/mB*(1.+alphas_muf*4./3./2./pi*(1.-L)+alphas_muf*4./3./16./pi*mB*(mB-2.*E_K)/E_K/E_K*DeltaFK/xi_P);
		
		//fT_fp=(mK+mB)/mB*(1.+alphas_muf*4./3./2./pi*(-logb+L)-alphas_muf*4./3./8./pi*mB/E_K*DeltaFK/xi_P);

		
		double lambda=pow(mB,4.)+pow(mK,4.)+q2*q2-2.*(mB*mB*mK*mK+mB*mB*q2+mK*mK*q2);
	
 		//FT=xi_P*2.*sqrt(lambda)*beta_l/(mB+mK)*CT*fT_fp;
 		//FT5=xi_P*2.*sqrt(lambda)*beta_l/(mB+mK)*CT5*fT_fp;

		/**********************Full FF**********************************/
			FA=(C10+C10p)*fp;
			FV=(C9+Y+C9p)*fp+2.*(mb+alphas_mub/3./pi*DeltaM)/(mB+mK)*(C7eff+C7effp+4.*ml/mb*CT)*fT+2.*(mb+alphas_mub/3./pi*DeltaM)/(mB+mK)*TauP;
			FS=(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*(CQ1+CQ1p)*f0;
			FP=(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*(CQ2+CQ2p)*f0-ml*(C10+C10p)*(fp-(mB*mB-mK*mK)/q2*(f0-fp));

			//FT=2.*sqrt(lambda)*beta_l/(mB+mK)*CT*fT;
			//FT5=2.*sqrt(lambda)*beta_l/(mB+mK)*CT5*fT;
		/**********************Full FF**********************************/

		/* hadronic uncertainties */		
		FV*=1.+param->BtoKlow_FV_err_noq2+q2/6.*param->BtoKlow_FV_err_q2;
		FA*=1.+param->BtoKlow_FA_err_noq2+q2/6.*param->BtoKlow_FA_err_q2;
		FS*=1.+param->BtoKlow_FS_err_noq2+q2/6.*param->BtoKlow_FS_err_q2;
		FP*=1.+param->BtoKlow_FP_err_noq2+q2/6.*param->BtoKlow_FP_err_q2;


		double Gamma_0=pow(cabs(param->Vtb*conj(param->Vts)),2.)*param->Gfermi*param->Gfermi*alpha_em*alpha_em/512./pow(pi,5.)/pow(mB,3.);

		double C_K=Gamma_0*sqrt(lambda)*beta_l;
		
		a_l=C_K*(q2*(beta_l*beta_l*cabs(FS*conj(FS))+cabs(FP*conj(FP)))+lambda/4.*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*(mB*mB-mK*mK+q2)*creal(FP*conj(FA))+4.*ml*ml*mB*mB*cabs(FA*conj(FA))); 
		
		b_l=2.*C_K*(q2*(beta_l*beta_l*creal(FS*conj(FT))+creal(FP*conj(FT5)))+ml*(sqrt(lambda)*beta_l*creal(FS*conj(FV))+(mB*mB-mK*mK+q2)*creal(FT5*conj(FA)))); 
		
		c_l=C_K*(q2*(beta_l*beta_l*cabs(FT*conj(FT))+cabs(FT5*conj(FT5)))-lambda/4.*beta_l*beta_l*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*sqrt(lambda)*beta_l*creal(FT*conj(FV))); 
	}
	

	if(q2>6.)
	{
		double complex C10=Cmub[10];
	
		double complex C7effp=Cpb[7];
		double complex C9p=Cpb[9];
		double complex C10p=Cpb[10];
	
		double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double complex CQ1p=CQpb[1];
		double complex CQ2p=CQpb[2];
	
		
		double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
		double shat=q2/mB/mB;
		double z=4.*mb*mb/q2;
		
		double complex x1=0.5+0.5*I*csqrt(z-1.);
		double complex x2=0.5-0.5*I*csqrt(z-1.);
		double complex x3=0.5+0.5*I/csqrt(z-1.);
		double complex x4=0.5-0.5*I/csqrt(z-1.);
		
		double complex A=
		-104./243.*log(mb*mb/mu_b/mu_b)+4.*shat/27./(1.-shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*log(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(csqrt(z-1.)))+9.*shat*shat*shat*log(shat)*log(shat)+18.*pi*I*shat*(1.-2.*shat)*log(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
		
		double complex B=
		8./243./shat*((4.-34.*shat-17.*pi*I*shat)*log(mb*mb/mu_b/mu_b)+8.*shat*pow(log(mb*mb/mu_b/mu_b),2.)+17.*shat*log(shat)*log(mb*mb/mu_b/mu_b))
		+(2.+shat)*csqrt(z-1.)/729./shat*(-48.*log(mb*mb/mu_b/mu_b)*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
		-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
		-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
		-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
		+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
		-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
		+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
		double complex C=-16./81.*log(q2/mu_b/mu_b)+428./243.-64./27.*zeta3+16./81.*pi*I;
		
		double complex C9eff=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
		+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
		+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(shat))
		+8.*mc*mc/q2*((4./9.*Cmub[1]+1./3.*Cmub[2])*(1.+param->Vub*conj(param->Vus)/param->Vtb/conj(param->Vts))+2.*Cmub[3]+20.*Cmub[5]);
		
		double complex C7eff=Cmub[7]
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(shat,log(mu_b/mb)));
		
		//double kappa=1.-2.*alphas_mub/3./pi*log(mu_b/mb);
		
		FA=(C10+C10p)*fp;
		FV=(C9eff+C9p)*fp+2.*mb/(mB+mK)*(C7eff+C7effp+4.*ml/mb*CT)*fT;
		FS=(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*(CQ1+CQ1p)*f0;
		FP=(mB*mB-mK*mK)/2./(mb_mub-param->mass_s)*(CQ2+CQ2p)*f0-ml*(C10+C10p)*(fp-(mB*mB-mK*mK)/q2*(f0-fp));

		//FT=2.*sqrt(lambda)*beta_l/(mB+mK)*CT*fT;
		//FT5=2.*sqrt(lambda)*beta_l/(mB+mK)*CT5*fT;
	
		/* hadronic uncertainties */		
		FV*=1.+param->BtoKhigh_FV_err;
		FA*=1.+param->BtoKhigh_FA_err;
		FS*=1.+param->BtoKhigh_FS_err;
		FP*=1.+param->BtoKhigh_FP_err;

		double lambda=pow(mB,4.)+pow(mK,4.)+q2*q2-2.*(mB*mB*mK*mK+mB*mB*q2+mK*mK*q2);
	
		double Gamma_0=pow(cabs(param->Vtb*conj(param->Vts)),2.)*param->Gfermi*param->Gfermi*alpha_em*alpha_em/512./pow(pi,5.)/pow(mB,3.);
	
		double C_K=Gamma_0*sqrt(lambda)*beta_l;
	
	
		double a_l_high=C_K*(q2*(beta_l*beta_l*cabs(FS*conj(FS))+cabs(FP*conj(FP)))+lambda/4.*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*(mB*mB-mK*mK+q2)*creal(FP*conj(FA))+4.*ml*ml*mB*mB*cabs(FA*conj(FA))); 
		
		double b_l_high=2.*C_K*(q2*(beta_l*beta_l*creal(FS*conj(FT))+creal(FP*conj(FT5)))+ml*(sqrt(lambda)*beta_l*creal(FS*conj(FV))+(mB*mB-mK*mK+q2)*creal(FT5*conj(FA)))); 
		
		double c_l_high=C_K*(q2*(beta_l*beta_l*cabs(FT*conj(FT))+cabs(FT5*conj(FT5)))-lambda/4.*beta_l*beta_l*(cabs(FA*conj(FA))+cabs(FV*conj(FV)))+2.*ml*sqrt(lambda)*beta_l*creal(FT*conj(FV))); 
	
		if(q2>14.)
		{
			a_l=a_l_high;
			b_l=b_l_high;
			c_l=c_l_high;
		}
		else
		{
			a_l=a_l*(14.-q2)/8.+a_l_high*(q2-6.)/8.;
			b_l=b_l*(14.-q2)/8.+b_l_high*(q2-6.)/8.;
			c_l=c_l*(14.-q2)/8.+c_l_high*(q2-6.)/8.;
		}

	}


	double dGamma_BKll_dq2=2.*(a_l+c_l/3.);

	double AFB[3],FH[3];

	AFB[0]=b_l/dGamma_BKll_dq2;
	AFB[1]=b_l;
	AFB[2]=dGamma_BKll_dq2;

	FH[0]=2.*(a_l+c_l)/dGamma_BKll_dq2;
	FH[1]=2.*(a_l+c_l);
	FH[2]=dGamma_BKll_dq2;

	for(je=0;je<=Nobs_BKll;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=FH[ie];
	}

	return dGamma_BKll_dq2;
}

/*----------------------------------------------------------------------*/

double dGamma_BKll_dq2_full_calculator(int gen, int charge, double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(B->K mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(gen,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(gen,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(gen,Cpb,CQpb,mu_W,mu_b,&param);

	return dGamma_BKll_dq2_full(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_full(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=10;
	if((smin<0.099)||(smax-smin>10.)) nmax=100;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_BKll+1],obs_den[Nobs_BKll+1];
	for(je=0;je<=Nobs_BKll;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated FH */
	
	double dobs[Nobs_BKll+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=s0m=0.;
	s0p=1.;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGamma_BKll_dq2_full(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	
	for(ie=1;ie<nmax;ie++)	
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s+=h;
		s0p=s;
		
		Gamma+=4.*dGamma_BKll_dq2_full(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKll;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		Gamma+=2.*dGamma_BKll_dq2_full(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		for(je=1;je<=Nobs_BKll;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}

		if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	
	s=smax;
	Gamma+=4.*dGamma_BKll_dq2_full(gen,charge,s-h/2.,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGamma_BKll_dq2_full(gen,charge,s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_BKll;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	
	obs[0]=s0;
	for(je=1;je<=Nobs_BKll;je++) obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_BKll;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	if(charge==0) return param->life_Bd/hbar*Gamma;
	else return param->life_B/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double RK_BKll_full(int charge, double smin, double smax, double complex C0be[], double complex C1be[], double complex C2be[], double complex CQ0be[], double complex CQ1be[], double complex Cpbe[], double complex CQpbe[], double complex C0bmu[], double complex C1bmu[], double complex C2bmu[], double complex CQ0bmu[], double complex CQ1bmu[], double complex Cpbmu[], double complex CQpbmu[], struct parameters* param, double mu_b)
{
	double obs[Nobs_BKll+1];
	
	return BRBKll_full(2,charge,smin,smax,obs,C0bmu,C1bmu,C2bmu,CQ0bmu,CQ1bmu,Cpbmu,CQpbmu,param,mu_b)/BRBKll_full(1,charge,smin,smax,obs,C0be,C1be,C2be,CQ0be,CQ1be,Cpbe,CQpbe,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_lowq2_full(int gen,int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKll_full(gen,charge,1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_highq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKll_full(gen,charge,14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKmumu_lowq2_full_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->K mu+ mu-) and all the other observables */

	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKll_lowq2_full(2,1,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKmumu_highq2_full_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->K mu+ mu-) and all the other observables */
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKll_highq2_full(2,1,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double RK_BKll_full_calculator(int charge, double smin, double smax, char name[])
{
	double complex C0bmu[11],C1bmu[11],C2bmu[11],C0wmu[11],C1wmu[11],C2wmu[11],Cpbmu[11],CQ0bmu[3],CQ1bmu[3],CQpbmu[3];
	double complex C0be[11],C1be[11],C2be[11],C0we[11],C1we[11],C2we[11],Cpbe[11],CQ0be[3],CQ1be[3],CQpbe[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(2,C0wmu,C1wmu,C2wmu,mu_W,&param);
	C_calculator_base1(C0wmu,C1wmu,C2wmu,mu_W,C0bmu,C1bmu,C2bmu,mu_b,&param);
	CQ_calculator(2,CQ0bmu,CQ1bmu,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpbmu,CQpbmu,mu_W,mu_b,&param);

	CW_calculator(1,C0we,C1we,C2we,mu_W,&param);
	C_calculator_base1(C0we,C1we,C2we,mu_W,C0be,C1be,C2be,mu_b,&param);
	CQ_calculator(1,CQ0be,CQ1be,mu_W,mu_b,&param);
	Cprime_calculator(1,Cpbe,CQpbe,mu_W,mu_b,&param);
	
	return RK_BKll_full(charge,smin,smax,C0be,C1be,C2be,CQ0be,CQ1be,Cpbe,CQpbe,C0bmu,C1bmu,C2bmu,CQ0bmu,CQ1bmu,Cpbmu,CQpbmu,&param,mu_b);
}

/*----------------------------------------------------------------------*/
/*---------------------------- WRAPPER ---------------------------------*/
/*----------------------------------------------------------------------*/

double dGamma_BKll_dq2(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return dGamma_BKll_dq2_full(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return dGamma_BKll_dq2_soft(gen,charge,q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/
	
double dGamma_BKll_dq2_calculator(int gen, int charge, double q2, double obs[][3], char name[])
{
	return dGamma_BKll_dq2_full_calculator(gen,charge,q2,obs,name);
}

/*----------------------------------------------------------------------*/

double BRBKll(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBKll_full(gen,charge,smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBKll_soft(gen,charge,smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double RK_BKll(int charge, double smin, double smax, double complex C0be[], double complex C1be[], double complex C2be[], double complex CQ0be[], double complex CQ1be[], double complex Cpbe[], double complex CQpbe[], double complex C0bmu[], double complex C1bmu[], double complex C2bmu[], double complex CQ0bmu[], double complex CQ1bmu[], double complex Cpbmu[], double complex CQpbmu[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return RK_BKll_full(charge,smin,smax,C0be,C1be,C2be,CQ0be,CQ1be,Cpbe,CQpbe,C0bmu,C1bmu,C2bmu,CQ0bmu,CQ1bmu,Cpbmu,CQpbmu,param,mu_b);
	else return RK_BKll_soft(charge,smin,smax,C0be,C1be,C2be,CQ0be,CQ1be,Cpbe,CQpbe,C0bmu,C1bmu,C2bmu,CQ0bmu,CQ1bmu,Cpbmu,CQpbmu,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_lowq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBKll_lowq2_full(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBKll_lowq2_soft(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKll_highq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	if(param->fullFF==1) return BRBKll_highq2_full(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else return BRBKll_highq2_soft(gen,charge,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKmumu_lowq2_calculator(char name[], double obs[])
{
	return BRobs_BKmumu_lowq2_full_calculator(name,obs);
}

/*----------------------------------------------------------------------*/

double BRobs_BKmumu_highq2_calculator(char name[], double obs[])
{
		return BRobs_BKmumu_highq2_full_calculator(name,obs);
}

/*----------------------------------------------------------------------*/

double RK_BKll_calculator(int charge, double smin, double smax,char name[])
{
		return RK_BKll_full_calculator(charge,smin,smax,name);
}
