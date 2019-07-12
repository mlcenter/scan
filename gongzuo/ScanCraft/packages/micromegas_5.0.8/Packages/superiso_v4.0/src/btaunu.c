#include "include.h"

double Blnu(int gen, struct parameters* param)
/* computes the branching ratio of B-> ell nu */
{
	double ml,lambda_l;
	if(gen==1) {ml=param->mass_e; lambda_l=param->lambda_l[1][1];}
	else if(gen==2) {ml=param->mass_mu; lambda_l=param->lambda_l[2][2];}
	else {ml=param->mass_tau; lambda_l=param->lambda_l[3][3];}

	double Vub=cabs(param->Vub);

	if(param->SM==1) return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*ml*param->f_B*(1.-ml*ml/param->m_B/param->m_B),2.);

	if(param->THDM_model==0) return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*ml*param->f_B*(1.-ml*ml/param->m_B/param->m_B)	
	*(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta)),2.);

	else return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*ml*param->f_B*(1.-ml*ml/param->m_B/param->m_B)	
	*(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*lambda_l*param->lambda_d[3][3]),2.);
}

/*--------------------------------------------------------------------*/

double Blnu_calculator(int gen, char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> ell nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Blnu(gen,&param);
}

/*--------------------------------------------------------------------*/

double RBlnu(int gen, struct parameters* param)
/* computes the ratio of BR(B-> ell nu)_MSSM/BR(B-> ell nu)_SM */
{
	if(param->SM==1) return 1.;
	
	double lambda_l;
	if(gen==1) lambda_l=param->lambda_l[1][1];
	else if(gen==2) lambda_l=param->lambda_l[2][2];
	else lambda_l=param->lambda_l[3][3];

	
	if(param->THDM_model==0) return pow(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta),2.);

	else return pow(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*lambda_l*param->lambda_d[3][3],2.);
}

/*--------------------------------------------------------------------*/

double RBlnu_calculator(int gen, char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> ell nu)_MSSM/BR(B-> ell nu)_SM */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return RBlnu(gen,&param);
}

/*--------------------------------------------------------------------*/

double Btaunu(struct parameters* param)
/* computes the branching ratio of B-> tau nu */
{
	return Blnu(3,param);
}

/*--------------------------------------------------------------------*/

double Btaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu) */
{
	return Blnu_calculator(3,name);
}

/*--------------------------------------------------------------------*/

double RBtaunu(struct parameters* param)
/* computes the ratio of BR(B-> tau nu)_MSSM/BR(B-> tau nu)_SM */
{
	return RBlnu(3,param);
}

/*--------------------------------------------------------------------*/

double RBtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu)_MSSM/BR(B-> tau nu)_SM */
{
	return RBlnu_calculator(3,name);
}
