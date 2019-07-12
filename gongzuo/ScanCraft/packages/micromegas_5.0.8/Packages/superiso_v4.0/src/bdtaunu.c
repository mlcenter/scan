#include "include.h"

/* From 1309.0301 & 1212.1878 */

/*--------------------------------------------------------------------*/

double dGammaBDlnu_dq2(int gen, int charge, double q2, double obs[][3], struct parameters* param)
{
	double ml,lambda_l;
	if(gen==1) {ml=param->mass_e; lambda_l=param->lambda_l[1][1];}
	else if(gen==2) {ml=param->mass_mu; lambda_l=param->lambda_l[2][2];}
	else {ml=param->mass_tau; lambda_l=param->lambda_l[3][3];}
	
	double mB,mD;
	
	if(charge==0)
	{
		mB=param->m_B;
		mD=param->m_D0;
	}
	else
	{
		mB=param->m_Bd;
		mD=param->m_D;
	}

	double mc=running_mass(param->mass_c,param->mass_c,mB,param->mass_top_pole,param->mass_b,param);
	double mb=running_mass(param->mass_b,param->mass_b,mB,param->mass_top_pole,param->mass_b,param);
	
	double complex C_V1,C_V2,C_T,C_S1,C_S2;
	
	C_V1=C_V2=C_T=0.;
	
	if(param->SM==1) C_S1=C_S2=0.;
	
	else if(param->THDM_model==0)
	{
		C_S1=-ml*mb/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta);
		C_S2=-ml*mc/param->mass_H/param->mass_H/(1.+epsilon_0(param)*param->tan_beta);
	}
	else 
	{	
		C_S1=-ml*mb/param->mass_H/param->mass_H*param->lambda_d[3][3]*lambda_l;
		C_S2=ml*mc/param->mass_H/param->mass_H*param->lambda_u[2][2]*lambda_l;
	}
	 
	double lambda_D=((mB-mD)*(mB-mD)-q2)*((mB+mD)*(mB+mD)-q2);

	double Delta=param->Delta_BD;
	double r_D=mD/mB;
	double w=(mB*mB+mD*mD-q2)/2./mB/mD;
	double z=(sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
	double rho_D2=param->rho_D2_BD;
	double V1_1=param->V1_1_BD;

	double V1=V1_1*(1.-8.*rho_D2*z+(51.*rho_D2-10.)*z*z-(252.*rho_D2-84.)*z*z*z);
	double S1=V1*(1.+Delta*(-0.019+0.041*(w-1.)-0.015*(w-1.)*(w-1.)));

	double hp=1./2./(1.+r_D*r_D-2.*r_D*w)*(-(1.+r_D)*(1.+r_D)*(w-1.)*V1+(1.-r_D)*(1.-r_D)*(w+1.)*S1);
	double hm=(1.-r_D*r_D)*(w+1.)/2./(1.+r_D*r_D-2.*r_D*w)*(S1-V1);
	double hT=(mb+mc)/(mB+mD)*(hp-(1.+r_D)/(1.-r_D)*hm);

	double F_1=1./2./sqrt(mB*mD)*((mB+mD)*hp-(mB-mD)*hm);
	double F_0=1./2./sqrt(mB*mD)*(((mB+mD)*(mB+mD)-q2)/(mB+mD)*hp-((mB-mD)*(mB-mD)-q2)/(mB-mD)*hm);
	double F_T=(mB+mD)/2./sqrt(mB*mD)*hT;
	
	double Hs_V0=sqrt(lambda_D/q2)*F_1;
	double Hs_Vt=(mB*mB-mD*mD)/sqrt(q2)*F_0;
	double Hs_S=(mB*mB-mD*mD)/(mb-mc)*F_0;
	double Hs_T=-sqrt(lambda_D)/(mB+mD)*F_T;
	
	int ie,je;
	for(je=0;je<=Nobs_BDlnu;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;

	double dGamma_dq2=pow(param->Gfermi*cabs(param->Vcb),2.)/192./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_D)*pow(1.-ml*ml/q2,2.)*
	(pow(cabs(1.+C_V1+C_V2),2.)*((1.+ml*ml/2./q2)*Hs_V0*Hs_V0+3./2.*ml*ml/q2*Hs_Vt*Hs_Vt)
	+3./2.*pow(cabs(C_S1+C_S2),2.)*Hs_S*Hs_S+8.*pow(cabs(C_T),2.)*(1.+2.*ml*ml/q2)*Hs_T*Hs_T
	+3.*creal((1.+C_V1+C_V2)*conj(C_S1+C_S2))*ml/sqrt(q2)*Hs_S*Hs_Vt
	-12.*creal((1.+C_V1+C_V2)*conj(C_T))*ml/sqrt(q2)*Hs_T*Hs_V0);
	
	double b_theta=pow(param->Gfermi*cabs(param->Vcb),2.)/128./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_D)*pow(1.-ml*ml/q2,2.)*
	(pow(cabs(1.+C_V1+C_V2),2.)*ml*ml/q2*Hs_V0*Hs_Vt
	+creal((1.+C_V1+C_V2)*conj(C_S1+C_S2))*ml/sqrt(q2)*Hs_S*Hs_V0
	-4.*creal((1.+C_V1+C_V2)*conj(C_T))*ml/sqrt(q2)*Hs_T*Hs_Vt
	-4.*creal((C_S1+C_S2)*conj(C_T))*Hs_T*Hs_S);
	
	double dGamma_dq2_ltau_m12=pow(param->Gfermi*cabs(param->Vcb),2.)/192./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_D)*pow(1.-ml*ml/q2,2.)*
	(pow(cabs(1.+C_V1+C_V2),2.)*Hs_V0*Hs_V0+16.*pow(cabs(C_T),2.)*ml*ml/q2*Hs_T*Hs_T
	-8.*creal((1.+C_V1+C_V2)*conj(C_T))*ml/sqrt(q2)*Hs_T*Hs_V0);
	
	double AFB[3],Ptau[3];
	
	AFB[0]=b_theta/dGamma_dq2;
	AFB[1]=b_theta;
	AFB[2]=dGamma_dq2;
	
	Ptau[0]=(1.-2.*dGamma_dq2_ltau_m12)/dGamma_dq2;
	Ptau[1]=dGamma_dq2-2.*dGamma_dq2_ltau_m12;
	Ptau[2]=dGamma_dq2;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=Ptau[ie];
	}
	
	return dGamma_dq2;
}

/*--------------------------------------------------------------------*/

double BRBDlnu(int gen, int charge, double smin, double smax, double obs[], struct parameters* param)
{
	int ie,je;
	int nmax=10;
	double Gamma=0.;
	double s;

	double obs_num[Nobs_BDlnu+1],obs_den[Nobs_BDlnu+1];
	for(je=0;je<=Nobs_BDlnu;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated Ptau */

	double dobs[Nobs_BDlnu+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=s0m=0.;
	s0p=1.;
	
	/*Gamma=dGammaBDlnu_dq2(gen,charge,smin,dobs,param);/2.;
	for(je=1;je<=Nobs_BDlnu;je++) 
	{
		obs_num[je]+=dobs[je][1]/2.;
		obs_den[je]+=dobs[je][2]/2.;
	}
	Gamma+=dGammaBDlnu_dq2(gen,charge,smax,dobs,param);/2.;
	for(je=1;je<=Nobs_BDlnu;je++) 
	{
		obs_num[je]+=dobs[je][1]/2.;
		obs_den[je]+=dobs[je][2]/2.;
	}
	
	for(ie=1;ie<nmax;ie++)
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s=smin+(smax-smin)*ie/nmax;
		s0p=s;	

		Gamma+=dGammaBDlnu_dq2(gen,charge,s,dobs,param);
		
		for(je=1;je<=Nobs_BDlnu;je++) 
		{
			obs_num[je]+=dobs[je][1];
			obs_den[je]+=dobs[je][2];
		}
		
		if(ie>1){if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);}
	}
	Gamma*=(smax-smin)/nmax;*/
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGammaBDlnu_dq2(gen,charge,s,dobs,param);;
	for(je=1;je<=Nobs_BDlnu;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	
	for(ie=1;ie<nmax;ie++)	
	{
		s+=h;

		Gamma+=4.*dGammaBDlnu_dq2(gen,charge,s-h/2.,dobs,param);;
		for(je=1;je<=Nobs_BDlnu;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		if(ie>1){if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);}

		Gamma+=2.*dGammaBDlnu_dq2(gen,charge,s,dobs,param);;
		for(je=1;je<=Nobs_BDlnu;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}

		if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	
	s=smax;
	Gamma+=4.*dGammaBDlnu_dq2(gen,charge,s-h/2.,dobs,param);
	for(je=1;je<=Nobs_BDlnu;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGammaBDlnu_dq2(gen,charge,s,dobs,param);
	for(je=1;je<=Nobs_BDlnu;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_BDlnu;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	
	obs[0]=s0;
	for(je=1;je<=Nobs_BDlnu;je++) obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_BDlnu;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	if(charge==0) return param->life_B/hbar*Gamma;
	else return param->life_Bd/hbar*Gamma;	
}

/*--------------------------------------------------------------------*/

double BRBDlnu_calculator(int gen, int charge, double q2min, double q2max, double obs[], char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 ell nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BRBDlnu(gen,charge,q2min,q2max,obs,&param);
}

/*--------------------------------------------------------------------*/

double BRBDlnu_full(int gen, int charge, double obs[], struct parameters* param)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==2) ml=param->mass_mu;
	else ml=param->mass_tau;
	
	double mB,mD;
	if(charge==0)
	{
		mB=param->m_B;
		mD=param->m_D0;
	}
	else
	{
		mB=param->m_Bd;
		mD=param->m_D;
	}

	return BRBDlnu(gen,charge,ml*ml,0.999*pow(mB-mD,2.),obs,param);
}

/*--------------------------------------------------------------------*/

double BRBDlnu_full_calculator(int gen, int charge, double obs[], char name[])
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double ml;
	if(gen==1) ml=param.mass_e;
	else if(gen==2) ml=param.mass_mu;
	else ml=param.mass_tau;
	
	double mB,mD;
	if(charge==0)
	{
		mB=param.m_B;
		mD=param.m_D0;
	}
	else
	{
		mB=param.m_Bd;
		mD=param.m_D;
	}

	return BRBDlnu(gen,charge,ml*ml,0.999*pow(mB-mD,2.),obs,&param);
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

double dGammaBDstarlnu_dq2(int gen, int charge, double q2, double obs[][3], struct parameters* param)
{
	double ml,lambda_l;
	if(gen==1) {ml=param->mass_e; lambda_l=param->lambda_l[1][1];}
	else if(gen==2) {ml=param->mass_mu; lambda_l=param->lambda_l[2][2];}
	else {ml=param->mass_tau; lambda_l=param->lambda_l[3][3];}
	
	double mB,mDs;
	if(charge==0)
	{
		mB=param->m_B;
		mDs=param->m_Dstar0;
	}
	else
	{
		mB=param->m_Bd;
		mDs=param->m_Dstar;
	}

	double mc=running_mass(param->mass_c,param->mass_c,mB,param->mass_top_pole,param->mass_b,param);
	double mb=running_mass(param->mass_b,param->mass_b,mB,param->mass_top_pole,param->mass_b,param);
	
	double complex C_V1,C_V2,C_T,C_S1,C_S2;
	
	C_V1=C_V2=C_T=0.;
	
	if(param->SM==1) C_S1=C_S2=0.;
	
	else if(param->THDM_model==0)
	{
		C_S1=-ml*mb/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta);
		C_S2=-ml*mc/param->mass_H/param->mass_H/(1.+epsilon_0(param)*param->tan_beta);
	}
	else 
	{	
		C_S1=-ml*mb/param->mass_H/param->mass_H*param->lambda_d[3][3]*lambda_l;
		C_S2=ml*mc/param->mass_H/param->mass_H*param->lambda_u[2][2]*lambda_l;
	}
	 
	double lambda_Ds=((mB-mDs)*(mB-mDs)-q2)*((mB+mDs)*(mB+mDs)-q2);

	//double Delta=param->Delta_BDstar;
	double r_Dstar=mDs/mB;
	double w=(mB*mB+mDs*mDs-q2)/2./mB/mDs;
	double z=(sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
	double rho_Dstar2=param->rho_Dstar2_BDstar;
	double R1_1=param->R1_1_BDstar;
	double R2_1=param->R2_1_BDstar;
	double R3_1=param->R3_1_BDstar;
	//double V1_1=param->V1_1_BDstar;
	double hA1_1=param->hA1_1_BDstar;
	
	//double V1=V1_1*(1.-8.*rho_Dstar2*z+(51.*rho_Dstar2-10.)*z*z-(252.*rho_Dstar2-84.)*z*z*z);
	double hA1=hA1_1*(1.-8.*rho_Dstar2*z+(53.*rho_Dstar2-15.)*z*z-(231.*rho_Dstar2-91.)*z*z*z);
	double R1=R1_1-0.12*(w-1.)+0.05*(w-1.)*(w-1.);
	double R2=R2_1-0.11*(w-1.)-0.06*(w-1.)*(w-1.);
	double R3=R3_1-0.052*(w-1.)+0.026*(w-1.)*(w-1.);
	//double S1=V1*(1.+Delta*(-0.019+0.041*(w-1.)-0.015*(w-1.)*(w-1.)));

	//double hp=1./2./(1.+r_Dstar*r_Dstar-2.*r_Dstar*w)*(-(1.+r_Dstar)*(1.+r_Dstar)*(w-1.)*V1+(1.-r_Dstar)*(1.-r_Dstar)*(w+1.)*S1);
	//double hm=(1.-r_Dstar*r_Dstar)*(w+1.)/2./(1.+r_Dstar*r_Dstar-2.*r_Dstar*w)*(S1-V1);
	double hV=R1*hA1;
	double hA2=(R2-R3)/2./r_Dstar*hA1;
	double hA3=(R2+R3)/2.*hA1;
	//double hT=(mb+mc)/(mB+mDs)*(hp-(1.+r_Dstar)/(1.-r_Dstar)*hm);
	double hT1=1./2./(1.+r_Dstar*r_Dstar-2.*r_Dstar*w)*((mb-mc)/(mB-mDs)*(1.-r_Dstar)*(1.-r_Dstar)*(w+1.)*hA1-(mb+mc)/(mB+mDs)*(1.+r_Dstar)*(1.+r_Dstar)*(w-1.)*hV);
	double hT2=(1.-r_Dstar*r_Dstar)*(w+1.)/2./(1.+r_Dstar*r_Dstar-2.*r_Dstar*w)*((mb-mc)/(mB-mDs)*hA1-(mb+mc)/(mB+mDs)*hV);
	double hT3=-1./2./(1.+r_Dstar)/(1.+r_Dstar*r_Dstar-2.*r_Dstar*w)*(2.*(mb-mc)/(mB-mDs)*r_Dstar*(w+1.)*hA1-(mb-mc)/(mB-mDs)*(1.+r_Dstar*r_Dstar-2.*r_Dstar*w)*(hA3-r_Dstar*hA2)-(mb+mc)/(mB+mDs)*(1.+r_Dstar)*(1.+r_Dstar)*hV);
 
	double V=(mB+mDs)/2./sqrt(mB*mDs)*hV;
	double A_1=((mB+mDs)*(mB+mDs)-q2)/2./sqrt(mB*mDs)/(mB+mDs)*hA1;
	double A_2=(mB+mDs)/2./sqrt(mB*mDs)*(hA3+mDs/mB*hA2);
	double A_0=1./2./sqrt(mB*mDs)*(((mB+mDs)*(mB+mDs)-q2)/2./mDs*hA1-(mB*mB-mDs*mDs+q2)/2./mB*hA2-(mB*mB-mDs*mDs-q2)/2./mDs*hA3);	
	double T_1=1./2./sqrt(mB*mDs)*((mB+mDs)*hT1-(mB-mDs)*hT2);
	double T_2=1./2./sqrt(mB*mDs)*(((mB+mDs)*(mB+mDs)-q2)/(mB+mDs)*hT1-((mB-mDs)*(mB-mDs)-q2)/(mB-mDs)*hT2);
	double T_3=1./2./sqrt(mB*mDs)*((mB-mDs)*hT1-(mB+mDs)*hT2-2.*(mB*mB-mDs*mDs)/mB*hT3);
	
	double H_Vp=(mB+mDs)*A_1-sqrt(lambda_Ds)/(mB+mDs)*V;
	double H_Vm=(mB+mDs)*A_1+sqrt(lambda_Ds)/(mB+mDs)*V;
	double H_V0=(mB+mDs)/2./mDs/sqrt(q2)*(-(mB*mB-mDs*mDs-q2)*A_1+lambda_Ds/pow(mB+mDs,2)*A_2);
	double H_Vt=-sqrt(lambda_Ds/q2)*A_0;
	double H_S=-sqrt(lambda_Ds)/(mb+mc)*A_0;
	double H_Tp=1./sqrt(q2)*((mB*mB-mDs*mDs)*T_2+sqrt(lambda_Ds)*T_1);
	double H_Tm=1./sqrt(q2)*(-(mB*mB-mDs*mDs)*T_2+sqrt(lambda_Ds)*T_1);
	double H_T0=1./2./mDs*(-(mB*mB+3.*mDs*mDs-q2)*T_2+lambda_Ds/(mB*mB-mDs*mDs)*T_3);

	int ie,je;
	for(je=0;je<=Nobs_BDstarlnu;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;

	double dGamma_dq2=pow(param->Gfermi*cabs(param->Vcb),2.)/192./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_Ds)*pow(1.-ml*ml/q2,2.)*
	((pow(cabs(1.+C_V1),2.)+pow(cabs(C_V2),2.))*((1.+ml*ml/2./q2)*(H_Vp*H_Vp+H_Vm*H_Vm+H_V0*H_V0)+3./2.*ml*ml/q2*H_Vt*H_Vt)
	-2.*creal((1.+C_V1)*conj(C_V2))*((1.+ml*ml/2./q2)*(H_V0*H_V0+2.*H_Vp*H_Vm)+3./2.*ml*ml/q2*H_Vt*H_Vt)
	+3./2.*pow(cabs(C_S1-C_S2),2.)*H_S*H_S+8.*pow(cabs(C_T),2.)*(1.+2.*ml*ml/q2)*(H_Tp*H_Tp+H_Tm*H_Tm+H_T0*H_T0)
	+3.*creal((1.+C_V1-C_V2)*(conj(C_S1-C_S2)))*ml/sqrt(q2)*H_S*H_Vt
	-12.*creal((1.+C_V1)*conj(C_T))*ml/sqrt(q2)*(H_T0*H_V0+H_Tp*H_Vp-H_Tm*H_Vm)
	+12.*creal(C_V2*conj(C_T))*ml/sqrt(q2)*(H_T0*H_V0+H_Tp*H_Vm-H_Tm*H_Vp));
	
	double b_theta=pow(param->Gfermi*cabs(param->Vcb),2.)/128./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_Ds)*pow(1.-ml*ml/q2,2.)*
	(1./2.*(pow(cabs(1.+C_V1),2.)-pow(cabs(C_V2),2.))*(H_Vp*H_Vp-H_Vm*H_Vm)+pow(cabs(1.+C_V1-C_V2),2.)*ml*ml/q2*H_V0*H_Vt
	+8.*pow(cabs(C_T),2.)*ml*ml/q2*(H_Tp*H_Tp-H_Tm*H_Tm)
	+creal((1.+C_V1-C_V2)*(conj(C_S1-C_S2)))*ml/sqrt(q2)*H_S*H_V0
	-4.*creal((1.+C_V1)*conj(C_T))*ml/sqrt(q2)*(H_T0*H_Vt+H_Tp*H_Vp+H_Tm*H_Vm)
	+4.*creal(C_V2*conj(C_T))*ml/sqrt(q2)*(H_T0*H_Vt+H_Tp*H_Vm+H_Tm*H_Vp)
	-4.*creal((C_S1-C_S2)*conj(C_T))*H_T0*H_S);
	
	double dGamma_dq2_ltau_m12=pow(param->Gfermi*cabs(param->Vcb),2.)/192./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_Ds)*pow(1.-ml*ml/q2,2.)*
	((pow(cabs(1.+C_V1),2.)+pow(cabs(C_V2),2.))*(H_Vp*H_Vp+H_Vm*H_Vm+H_V0*H_V0)
	-2.*creal((1.+C_V1)*conj(C_V2))*(H_V0*H_V0+2.*H_Vp*H_Vm)
	+16.*pow(cabs(C_T),2.)*ml*ml/q2*(H_Tp*H_Tp+H_Tm*H_Tm+H_T0*H_T0)
	-8.*creal((1.+C_V1)*conj(C_T))*ml/sqrt(q2)*(H_T0*H_V0+H_Tp*H_Vp-H_Tm*H_Vm)
	+8.*creal(C_V2*conj(C_T))*ml/sqrt(q2)*(H_T0*H_V0+H_Tp*H_Vm-H_Tm*H_Vp));
	
	double dGamma_dq2_lDstar_0=pow(param->Gfermi*cabs(param->Vcb),2.)/192./pow(pi,3.)/pow(mB,3.)*q2*sqrt(lambda_Ds)*pow(1.-ml*ml/q2,2.)*
	((pow(cabs(1.+C_V1-C_V2),2.))*((1.+ml*ml/2./q2)*H_V0*H_V0+3./2.*ml*ml/q2*H_Vt*H_Vt)
	+3./2.*pow(cabs(C_S1-C_S2),2.)*H_S*H_S+8.*pow(cabs(C_T),2.)*(1.+2.*ml*ml/q2)*H_T0*H_T0
	+3.*creal((1.+C_V1-C_V2)*(conj(C_S1-C_S2)))*ml/sqrt(q2)*H_S*H_Vt
	-12.*creal((1.+C_V1-C_V2)*conj(C_T))*ml/sqrt(q2)*H_T0*H_V0);

	double AFB[3],Ptau[3],PDstar[3];
	
	AFB[0]=b_theta/dGamma_dq2;
	AFB[1]=b_theta;
	AFB[2]=dGamma_dq2;
	
	Ptau[0]=(1.-2.*dGamma_dq2_ltau_m12)/dGamma_dq2;
	Ptau[1]=dGamma_dq2-2.*dGamma_dq2_ltau_m12;
	Ptau[2]=dGamma_dq2;

	PDstar[0]=dGamma_dq2_lDstar_0/dGamma_dq2;
	PDstar[1]=dGamma_dq2_lDstar_0;
	PDstar[2]=dGamma_dq2;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=Ptau[ie];
		obs[3][ie]=PDstar[ie];
	}
	
	return dGamma_dq2;
}

/*--------------------------------------------------------------------*/

double BRBDstarlnu(int gen, int charge, double smin, double smax, double obs[], struct parameters* param)
{
	int ie,je;
	int nmax=10;
	double Gamma=0.;
	double s;

	double obs_num[Nobs_BDstarlnu+1],obs_den[Nobs_BDstarlnu+1];
	for(je=0;je<=Nobs_BDstarlnu;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated Ptau */
	obs[3]=0.; /* integrated PD* */

	double dobs[Nobs_BDstarlnu+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=s0m=0.;
	s0p=1.;
	
	/*Gamma=dGammaBDstarlnu_dq2(gen,charge,smin,dobs,param);/2.;
	for(je=1;je<=Nobs_BDstarlnu;je++) 
	{
		obs_num[je]+=dobs[je][1]/2.;
		obs_den[je]+=dobs[je][2]/2.;
	}
	Gamma+=dGammaBDstarlnu_dq2(gen,charge,smax,dobs,param);/2.;
	for(je=1;je<=Nobs_BDstarlnu;je++) 
	{
		obs_num[je]+=dobs[je][1]/2.;
		obs_den[je]+=dobs[je][2]/2.;
	}
	
	for(ie=1;ie<nmax;ie++)
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s=smin+(smax-smin)*ie/nmax;
		s0p=s;	

		Gamma+=dGammaBDstarlnu_dq2(gen,charge,s,dobs,param);
		
		for(je=1;je<=Nobs_BDstarlnu;je++) 
		{
			obs_num[je]+=dobs[je][1];
			obs_den[je]+=dobs[je][2];
		}
		
		if(ie>1){if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);}
	}
	Gamma*=(smax-smin)/nmax;*/
	
	double h=(smax-smin)/nmax;	
	s=smin;
	Gamma=dGammaBDstarlnu_dq2(gen,charge,s,dobs,param);;
	for(je=1;je<=Nobs_BDstarlnu;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}
	
	for(ie=1;ie<nmax;ie++)	
	{
		s+=h;

		Gamma+=4.*dGammaBDstarlnu_dq2(gen,charge,s-h/2.,dobs,param);;
		for(je=1;je<=Nobs_BDstarlnu;je++) 
		{
			obs_num[je]+=4.*dobs[je][1];
			obs_den[je]+=4.*dobs[je][2];
		}

		if(ie>1){if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);}

		Gamma+=2.*dGammaBDstarlnu_dq2(gen,charge,s,dobs,param);;
		for(je=1;je<=Nobs_BDstarlnu;je++) 
		{
			obs_num[je]+=2.*dobs[je][1];
			obs_den[je]+=2.*dobs[je][2];
		}

		if(dAFBtmp/dobs[1][0]<0.) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	
	s=smax;
	Gamma+=4.*dGammaBDstarlnu_dq2(gen,charge,s-h/2.,dobs,param);
	for(je=1;je<=Nobs_BDstarlnu;je++) 
	{
		obs_num[je]+=4.*dobs[je][1];
		obs_den[je]+=4.*dobs[je][2];
	}	
	Gamma+=dGammaBDstarlnu_dq2(gen,charge,s,dobs,param);
	for(je=1;je<=Nobs_BDstarlnu;je++) 
	{
		obs_num[je]+=dobs[je][1];
		obs_den[je]+=dobs[je][2];
	}	
	
	Gamma*=h/6.;
	for(je=1;je<=Nobs_BDstarlnu;je++) 
	{
		obs_num[je]*=h/6.;
		obs_den[je]*=h/6.;
	}
	
	obs[0]=s0;
	for(je=1;je<=Nobs_BDstarlnu;je++) obs[je]=obs_num[je]/obs_den[je];
	
	for(je=1;je<=Nobs_BDstarlnu;je++) if(fabs(obs[je])<1.e-15) obs[je]=0.;

	if(charge==0) return param->life_B/hbar*Gamma;
	else return param->life_Bd/hbar*Gamma;	}

/*--------------------------------------------------------------------*/

double BRBDstarlnu_calculator(int gen, int charge, double q2min, double q2max, double obs[], char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 ell nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BRBDstarlnu(gen,charge,q2min,q2max,obs,&param);
}

/*--------------------------------------------------------------------*/

double BRBDstarlnu_full(int gen, int charge, double obs[], struct parameters* param)
{	
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==2) ml=param->mass_mu;
	else ml=param->mass_tau;
	
	double mB,mDs;
	if(charge==0)
	{
		mB=param->m_B;
		mDs=param->m_Dstar0;
	}
	else
	{
		mB=param->m_Bd;
		mDs=param->m_Dstar;
	}

	return BRBDstarlnu(gen,charge,ml*ml,0.999*pow(mB-mDs,2.),obs,param);
}

/*--------------------------------------------------------------------*/

double BRBDstarlnu_full_calculator(int gen, int charge, double obs[], char name[])
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double ml;
	if(gen==1) ml=param.mass_e;
	else if(gen==2) ml=param.mass_mu;
	else ml=param.mass_tau;
	
	double mB,mDs;
	if(charge==0)
	{
		mB=param.m_B;
		mDs=param.m_Dstar0;
	}
	else
	{
		mB=param.m_Bd;
		mDs=param.m_Dstar;
	}

	return BRBDstarlnu(gen,charge,ml*ml,0.999*pow(mB-mDs,2.),obs,&param);
}

/*--------------------------------------------------------------------*/

double BDtaunu(struct parameters* param)
/* computes the branching ratio of B-> D+ tau nu */
{
	double obs[Nobs_BDlnu+1];
	return BRBDlnu_full(3,1,obs,param);
}

/*--------------------------------------------------------------------*/

double BDtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D+ tau nu) */
{
	double obs[Nobs_BDlnu+1];
	return BRBDlnu_full_calculator(3,1,obs,name);
}

/*--------------------------------------------------------------------*/

double BDtaunu_BDenu(struct parameters* param)
/* computes the ratio BR(B-> D+ tau nu)/BR(B-> D+ e nu) */
{
	double obs_tau[Nobs_BDlnu+1],obs_e[Nobs_BDlnu+1];
	return BRBDlnu_full(3,1,obs_tau,param)/BRBDlnu_full(1,1,obs_e,param);
}

/*--------------------------------------------------------------------*/

double BDtaunu_BDenu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D+ tau nu)/BR(B-> D+ e nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BDtaunu_BDenu(&param);
}

/*--------------------------------------------------------------------*/

double BDstartaunu(struct parameters* param)
/* computes the branching ratio of B-> Dstar+ tau nu */
{
	double obs[Nobs_BDstarlnu+1];
	return BRBDstarlnu_full(3,1,obs,param);
}

/*--------------------------------------------------------------------*/

double BDstartaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> Dstar+ tau nu) */
{
	double obs[Nobs_BDstarlnu+1];	
	return BRBDstarlnu_full_calculator(3,1,obs,name);
}

/*--------------------------------------------------------------------*/

double BDstartaunu_BDstarenu(struct parameters* param)
/* computes the ratio BR(B-> Dstar+ tau nu)/BR(B-> Dstar+ e nu) */
{
	double obs_tau[Nobs_BDstarlnu+1],obs_e[Nobs_BDstarlnu+1];
	return BRBDstarlnu_full(3,1,obs_tau,param)/BRBDstarlnu_full(1,1,obs_e,param);
}

/*--------------------------------------------------------------------*/

double BDstartaunu_BDstarenu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> Dstar+ tau nu)/BR(B-> Dstar+ e nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BDstartaunu_BDstarenu(&param);
}
