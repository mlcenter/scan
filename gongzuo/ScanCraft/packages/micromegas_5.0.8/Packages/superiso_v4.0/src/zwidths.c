#include "include.h"

double Zwidth_neu1neu1(struct parameters* param)
/* computes the partial Z width to neutralino 1 + neutralino 1 (in MeV) */
{
	if(param->mass_Z<2.*fabs(param->mass_neut[1])) return 0.;

	double alpha=1./137.03599011;

	double sw=sqrt(1.-param->mass_W*param->mass_W/param->mass_Z/param->mass_Z);
	double cw=sqrt(1.-sw*sw);
	
	double complex gA=(2.*param->neut_mix[1][4]*conj(param->neut_mix[1][4])-2.*param->neut_mix[1][3]*conj(param->neut_mix[1][3]))/4./sw/cw;

	double RA=pow(1.-4.*param->mass_neut[1]*param->mass_neut[1]/param->mass_Z/param->mass_Z,1.5);
	
	return alpha/3.*param->mass_Z*pow(cabs(gA),2.)*RA*1000.;
}

/*--------------------------------------------------------------------*/

double Zwidth_invisible(struct parameters* param)
/* computes the partial Z width to neutrinos (in MeV) */
{
	double sw=sin(atan(param->gp/param->g2));
	double cw=sqrt(1.-sw*sw);
	
	double gV=0.5/2./cw; /* I3(neutrino) = 0.5 */
	double gA=0.5/2./cw; /* Q(neutrino) = 0 */
	
	double RV=1.; /* neutrino masses negligible */
	double RA=1.;

	double width=param->g2*param->g2/12./pi*param->mass_Z*(gV*gV*RV+gA*gA*RA);

	return 3.*width*1000.+Zwidth_neu1neu1(param);
}

/*--------------------------------------------------------------------*/

double Zwidth_b1b1(struct parameters* param)
/* computes the partial Z width to sbottom 1 + anti-bottom 1 (in MeV) */
{
	if(param->mass_Z<2.*fabs(param->mass_b1)) return 0.;

	double alpha=1./137.03599011;
	
	double sw=sin(atan(param->gp/param->g2));
	double cw=sqrt(1.-sw*sw);
	
	return param->mass_Z/16.*alpha/sw/sw/cw/cw*pow(param->sbot_mix[1][1]*conj(param->sbot_mix[1][1])-2./3.*sw*sw,2.)*pow(1.-4.*param->mass_b1*param->mass_b1/param->mass_Z/param->mass_Z,1.5)*1000.;
}

/*--------------------------------------------------------------------*/

double Zwidth_neu1neu1_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the partial Z width to neutralino 1 + neutralino 1 (in MeV) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Zwidth_neu1neu1(&param);
}

/*--------------------------------------------------------------------*/

double Zwidth_invisible_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the total partial Z width to neutrinos (in MeV) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Zwidth_invisible(&param);
}

/*--------------------------------------------------------------------*/

double Zwidth_b1b1_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the partial Z width to sbottom 1 + anti-sbottom 1 (in MeV) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Zwidth_b1b1(&param);
}
