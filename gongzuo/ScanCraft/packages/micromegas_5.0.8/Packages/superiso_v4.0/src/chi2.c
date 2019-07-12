#include "include.h"

//#define DEBUG

/*--------------------------------------------------------------------*/

int observables(int ke, double complex deltaC[], double complex deltaCp[], double complex deltaCQ[], double complex deltaCQp[], char* names[], double values[], int onoff[], double values_ref[], struct parameters* param0)
/* Contains all the observables with experimental data */
{		
	int nbobs=0;
	int nbobsall=0;
	
	int ie;
	
	int method; /* Remark: myrand(X,-1)=1, myrand(X,0)=0, myrand(X,1)=random number */
	if(ke==0) method=0;
	else if(ke>0) method=-1;
	else if(ke<0) method=1;

	if(ke>0) {for(ie=0;ie<NBOBSMAX;ie++) values[ie]=values_ref[ie];}
	else {for(ie=0;ie<NBOBSMAX;ie++) values[ie]=0.;}

	double complex C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3];
	double complex C0eb[11],C1eb[11],C0ew[11],C1ew[11],C2ew[11],C2eb[11],Cpeb[11],CQpeb[3],CQ0eb[3],CQ1eb[3];
	double obsKstar[Nobs_BKsll+1];
	double obsK[Nobs_BKll+1];
	double obsphi[Nobs_Bsphill+1];

	double lambda_h=0.5;

	struct parameters param=*param0;
		
	double sigma_FFKstar_full[21]={0.0289419,0.256849,1.63438,0.026356,0.187894,1.02531,0.0208033,0.128772,0.656273,0.0332944,0.261268,1.53102,0.027496,0.189575,1.63965,0.027496,0.166219,0.803783,0.0633333,0.222129,2.20338};
	double sigma_FFKstar_soft[9]={sigma_FFKstar_full[3],sigma_FFKstar_full[4],sigma_FFKstar_full[5],sigma_FFKstar_full[6],sigma_FFKstar_full[7],sigma_FFKstar_full[8],sigma_FFKstar_full[9],sigma_FFKstar_full[10],sigma_FFKstar_full[11]};
	double sigma_FFphi_full[21]={0.0240514,0.23757,1.35909,0.0105759,0.103813,0.790297,0.0152447,0.125879,0.47881,0.0141353,0.164176,1.72682,0.0120688,0.0835241,1.00324,0.0120688,0.0803563,0.608872,0.0357469,0.330082,1.79526};
	double sigma_FFphi_soft[9]={sigma_FFphi_full[3],sigma_FFphi_full[4],sigma_FFphi_full[5],sigma_FFphi_full[6],sigma_FFphi_full[7],sigma_FFphi_full[8],sigma_FFphi_full[9],sigma_FFphi_full[10],sigma_FFphi_full[11]};
	double sigma_FFK[10]={0.03,0.1,1.07,2.74,0.02,0.09,0.76,0.02,0.13,1.00};


///* nuisance for all observables */
	param.alphas_MZ=0.1181+0.0011*myrand(1,method)*(method==1||ke==1);

	param.mass_b=4.18+0.04*myrand(1,method)*(method==1||ke==2);
	param.mass_c=1.27+0.003*myrand(1,method)*(method==1||ke==3);
	param.mass_s=0.096+0.008*myrand(1,method)*(method==1||ke==4);
	param.mass_top_pole=173.34+sqrt(0.27*0.27+0.71*0.71)*myrand(1,method)*(method==1||ke==5);

	param.mass_h0=125.09+0.24*myrand(1,method)*(method==1||ke==6);

	param.CKM_lambda=0.22506+0.00050*myrand(1,method)*(method==1||ke==7);
	param.CKM_A=0.811+0.026*myrand(1,method)*(method==1||ke==8);
	param.CKM_rhobar=0.124+0.019*myrand(1,method)*(method==1||ke==9);
	param.CKM_etabar=0.356+0.011*myrand(1,method)*(method==1||ke==10);

	slha_adjust(&param);

	double mu_W=param.mass_W*(1.25+0.75*myrand(2,method)*(method==1||ke==11));
	double mu_b=param.mass_b_pole*(1.25+0.75*myrand(2,method)*(method==1||ke==12));

/* inclusive b -> s */
	param.BR_BXclnu_exp=0.1065+0.0016*myrand(1,method)*(method==1||ke==13);
/* b -> s gamma */ 
	param.mu_G2_bsg=0.336+0.064*myrand(1,method)*(method==1||ke==14);
	param.rho_D3_bsg=0.153+0.045*myrand(1,method)*(method==1||ke==15);
	param.rho_LS3_bsg=-0.145+0.098*myrand(1,method)*(method==1||ke==16);
	double bsgamma_rand=0.+0.04*myrand(1,method)*(method==1||ke==17);
	param.mu_c_bsg=2.45+1.55*myrand(2,method)*(method==1||ke==18);
/* b -> s mu mu */
	double BRBXsmumu_lowq2_rand=0.+0.05*myrand(1,method)*(method==1||ke==19);
	double BRBXsmumu_highq2_rand=0.+0.05*myrand(1,method)*(method==1||ke==20);
	double BRBXsmumu_full_rand=0.+0.05*myrand(1,method)*(method==1||ke==21);
/* b -> s e e */
	double BRBXsee_lowq2_rand=0.+0.05*myrand(1,method)*(method==1||ke==22);
	double BRBXsee_highq2_rand=0.+0.05*myrand(1,method)*(method==1||ke==23);
	double BRBXsee_full_rand=0.+0.05*myrand(1,method)*(method==1||ke==24);

/* B */	
	param.f_B=0.1905+0.0042*myrand(1,method)*(method==1||ke==25);
	param.lambda_Bp=0.46+0.11*myrand(1,method)*(method==1||ke==26);
/* B -> K* */
	param.f_Kstar_par=0.204+0.007*myrand(1,method)*(method==1||ke==27);
	param.f_Kstar_perp=0.159+0.006*myrand(1,method)*(method==1||ke==28);
	param.a1perp=0.04+0.03*myrand(1,method)*(method==1||ke==29);
	param.a2perp=0.10+0.08*myrand(1,method)*(method==1||ke==30);
	param.a1par=0.06+0.04*myrand(1,method)*(method==1||ke==31);
	param.a2par=0.16+0.09*myrand(1,method)*(method==1||ke==32);
/* B -> K* gamma */
	param.T1_BKstar=0.312+0.027*myrand(1,method)*(method==1||ke==33);
	double mu_spec=sqrt(lambda_h*param.mass_b)*(1.25+0.75*myrand(2,method)*(method==1||ke==34));
/* low */
	switch(param.BKstar_implementation)
	{
		default:
		{
			param.BtoKstarlow_ALperp_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==35)*cexp(I*pi/4.);
			param.BtoKstarlow_ARperp_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==36)*cexp(I*pi/4.);
			param.BtoKstarlow_ALpar_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==37)*cexp(I*pi/4.);
			param.BtoKstarlow_ARpar_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==38)*cexp(I*pi/4.);
			param.BtoKstarlow_AL0_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==39)*cexp(I*pi/4.);
			param.BtoKstarlow_AR0_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==40)*cexp(I*pi/4.);
			param.BtoKstarlow_At_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==41)*cexp(I*pi/4.);
			param.BtoKstarlow_AS_err_noq2=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==42)*cexp(I*pi/4.);

			param.BtoKstarlow_ALperp_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==43)*cexp(I*pi/4.);
			param.BtoKstarlow_ARperp_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==44)*cexp(I*pi/4.);
			param.BtoKstarlow_ALpar_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==45)*cexp(I*pi/4.);
			param.BtoKstarlow_ARpar_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==46)*cexp(I*pi/4.);
			param.BtoKstarlow_AL0_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==47)*cexp(I*pi/4.);
			param.BtoKstarlow_AR0_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==48)*cexp(I*pi/4.);
			param.BtoKstarlow_At_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==49)*cexp(I*pi/4.);
			param.BtoKstarlow_AS_err_q2=2.5*param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==50)*cexp(I*pi/4.);
			
			break;
		}
		
		case 2:
		{
			param.real_alpha_perp0=-0.06e-4+0.21e-4*myrand(1,method)*(method==1||ke==35);
			param.real_alpha_perp1=-6.77e-4+0.27e-4*myrand(1,method)*(method==1||ke==36);
			param.real_alpha_perp2=18.96e-4+0.59e-4*myrand(1,method)*(method==1||ke==37);
			param.real_alpha_par0=-0.35e-4+0.62e-4*myrand(1,method)*(method==1||ke==38);
			param.real_alpha_par1=-3.13e-4+0.41e-4*myrand(1,method)*(method==1||ke==39);
			param.real_alpha_par2=12.20e-4+1.34e-4*myrand(1,method)*(method==1||ke==40);
			param.real_alpha_zero0=0.05e-4+1.52e-4*myrand(1,method)*(method==1||ke==41);
			param.real_alpha_zero1=17.26e-4+1.64e-4*myrand(1,method)*(method==1||ke==42);
			param.imag_alpha_perp0=-0.21e-4+2.25e-4*myrand(1,method)*(method==1||ke==43);
			param.imag_alpha_perp1=1.17e-4+3.58e-4*myrand(1,method)*(method==1||ke==44);
			param.imag_alpha_perp2=-0.08e-4+2.24e-4*myrand(1,method)*(method==1||ke==45);
			param.imag_alpha_par0=-0.04e-4+3.67e-4*myrand(1,method)*(method==1||ke==46);
			param.imag_alpha_par1=-2.14e-4+2.46e-4*myrand(1,method)*(method==1||ke==47);
			param.imag_alpha_par2=6.03e-4+2.50e-4*myrand(1,method)*(method==1||ke==48);
			param.imag_alpha_zero0=-0.05e-4+4.99e-4*myrand(1,method)*(method==1||ke==49);
			param.imag_alpha_zero1=4.29e-4+3.14e-4*myrand(1,method)*(method==1||ke==50);
			
			break;
		}
		
		case 3:
		{
			param.DeltaC9_M1_q2bar=0.72+0.47*myrand(1,method)*(method==1||ke==35);
			param.r1_M1=0.10+0.01*myrand(1,method)*(method==1||ke==36);
			param.r2_M1=1.13+0.01*myrand(1,method)*(method==1||ke==37);
			param.DeltaC9_M2_q2bar=0.76+0.56*myrand(1,method)*(method==1||ke==38);
			param.r1_M2=0.09+0.01*myrand(1,method)*(method==1||ke==39);
			param.r2_M2=1.12+0.01*myrand(1,method)*(method==1||ke==40);
			param.DeltaC9_M3_q2bar=1.11+0.92*myrand(1,method)*(method==1||ke==41);
			param.r1_M3=0.06+0.07*myrand(1,method)*(method==1||ke==42);
			param.r2_M3=1.05+0.05*myrand(1,method)*(method==1||ke==43);
			
			break;
		}
	}
	
/* high */
	param.BtoKstarhigh_ALperp_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==51)*cexp(I*pi/4.);
	param.BtoKstarhigh_ARperp_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==52)*cexp(I*pi/4.);
	param.BtoKstarhigh_ALpar_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==53)*cexp(I*pi/4.);
	param.BtoKstarhigh_ARpar_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==54)*cexp(I*pi/4.);
	param.BtoKstarhigh_AL0_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==55)*cexp(I*pi/4.);
	param.BtoKstarhigh_AR0_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==56)*cexp(I*pi/4.);
	param.BtoKstarhigh_At_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==57)*cexp(I*pi/4.);
	param.BtoKstarhigh_AS_err=param.hadrerrBKstar/100.*myrand(1,method)*(method==1||ke==58)*cexp(I*pi/4.);		

/* B -> K */
	param.f_K=0.156+0.005*myrand(1,method)*(method==1||ke==59);
	param.a1K=0.06+0.03*myrand(1,method)*(method==1||ke==60);
	param.a2K=0.25+0.15*myrand(1,method)*(method==1||ke==61);
	/* Form factors B->K ll */
	param.a00_BK=0.54+sigma_FFK[0]*myrand(1,method)*(method==1||ke==62);
	param.a10_BK=-1.91+sigma_FFK[1]*myrand(1,method)*(method==1||ke==63);
	param.a20_BK=1.83+sigma_FFK[2]*myrand(1,method)*(method==1||ke==64);
	param.a30_BK=-0.02+sigma_FFK[3]*myrand(1,method)*(method==1||ke==65);
	
	param.a0p_BK=0.43+sigma_FFK[4]*myrand(1,method)*(method==1||ke==66);
	param.a1p_BK=-0.67+sigma_FFK[5]*myrand(1,method)*(method==1||ke==67);
	param.a2p_BK=-1.12+sigma_FFK[6]*myrand(1,method)*(method==1||ke==68);

	param.a0T_BK=0.4+sigma_FFK[7]*myrand(1,method)*(method==1||ke==69);
	param.a1T_BK=-0.53+sigma_FFK[8]*myrand(1,method)*(method==1||ke==70);
	param.a2T_BK=-0.29+sigma_FFK[9]*myrand(1,method)*(method==1||ke==71);
	
/* low */
	param.BtoKlow_FV_err_noq2=param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==72)*cexp(I*pi/4.);
	param.BtoKlow_FA_err_noq2=param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==73)*cexp(I*pi/4.);
	param.BtoKlow_FS_err_noq2=param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==74)*cexp(I*pi/4.);
	param.BtoKlow_FP_err_noq2=param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==75)*cexp(I*pi/4.);

	param.BtoKlow_FV_err_q2=2.5*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==76)*cexp(I*pi/4.);
	param.BtoKlow_FA_err_q2=2.5*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==77)*cexp(I*pi/4.);
	param.BtoKlow_FS_err_q2=2.5*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==78)*cexp(I*pi/4.);
	param.BtoKlow_FP_err_q2=2.5*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==79)*cexp(I*pi/4.);

/* high */
	param.BtoKhigh_FV_err=param.hadrerrBK/100.*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==80)*cexp(I*pi/4.);
	param.BtoKhigh_FA_err=param.hadrerrBK/100.*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==81)*cexp(I*pi/4.);
	param.BtoKhigh_FS_err=param.hadrerrBK/100.*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==82)*cexp(I*pi/4.);
	param.BtoKhigh_FP_err=param.hadrerrBK/100.*param.hadrerrBK/100.*myrand(1,method)*(method==1||ke==83)*cexp(I*pi/4.);

/* Bs */
	param.life_Bs=1.511e-12+0.014e-12*myrand(1,method)*(method==1||ke==84);
	param.f_Bs=0.2277+0.0045*myrand(1,method)*(method==1||ke==85);
	param.lambda_Bsp=0.46+0.11*myrand(1,method)*(method==1||ke==86);
	
/* Bs -> phi */
	param.f_phi_par=0.233+0.004*myrand(1,method)*(method==1||ke==87);
	param.f_phi_perp=0.191+0.004*myrand(1,method)*(method==1||ke==88);
	param.a1phi_perp=0.;
	param.a1phi_par=0.;
	param.a2phi_perp=0.14+0.07*myrand(1,method)*(method==1||ke==89);
	param.a2phi_par=0.23+0.08*myrand(1,method)*(method==1||ke==90);

/* low */
	param.Bstophilow_ALperp_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==91)*cexp(I*pi/4.);
	param.Bstophilow_ARperp_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==92)*cexp(I*pi/4.);
	param.Bstophilow_ALpar_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==93)*cexp(I*pi/4.);
	param.Bstophilow_ARpar_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==94)*cexp(I*pi/4.);
	param.Bstophilow_AL0_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==95)*cexp(I*pi/4.);
	param.Bstophilow_AR0_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==96)*cexp(I*pi/4.);
	param.Bstophilow_At_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==97)*cexp(I*pi/4.);
	param.Bstophilow_AS_err_noq2=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==98)*cexp(I*pi/4.);

	param.Bstophilow_ALperp_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==99)*cexp(I*pi/4.);
	param.Bstophilow_ARperp_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==100)*cexp(I*pi/4.);
	param.Bstophilow_ALpar_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==101)*cexp(I*pi/4.);
	param.Bstophilow_ARpar_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==102)*cexp(I*pi/4.);
	param.Bstophilow_AL0_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==103)*cexp(I*pi/4.);
	param.Bstophilow_AR0_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==104)*cexp(I*pi/4.);
	param.Bstophilow_At_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==105)*cexp(I*pi/4.);
	param.Bstophilow_AS_err_q2=2.5*param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==106)*cexp(I*pi/4.);	
/* high */
	param.Bstophihigh_ALperp_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==107)*cexp(I*pi/4.);
	param.Bstophihigh_ARperp_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==108)*cexp(I*pi/4.);
	param.Bstophihigh_ALpar_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==109)*cexp(I*pi/4.);
	param.Bstophihigh_ARpar_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==110)*cexp(I*pi/4.);
	param.Bstophihigh_AL0_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==111)*cexp(I*pi/4.);
	param.Bstophihigh_AR0_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==112)*cexp(I*pi/4.);
	param.Bstophihigh_At_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==113)*cexp(I*pi/4.);
	param.Bstophihigh_AS_err=param.hadrerrBsphi/100.*myrand(1,method)*(method==1||ke==114)*cexp(I*pi/4.);
			
	if(param.fullFF)
	{
		/* Form factors B->K* ll */	
		param.a0A0_BKstar=0.369196+sigma_FFKstar_full[0]*myrand(1,method)*(method==1||ke==115);
		param.a1A0_BKstar=-1.36584+sigma_FFKstar_full[1]*myrand(1,method)*(method==1||ke==116);
		param.a2A0_BKstar=0.128191+sigma_FFKstar_full[2]*myrand(1,method)*(method==1||ke==117);
		param.a0A1_BKstar=0.29725+sigma_FFKstar_full[3]*myrand(1,method)*(method==1||ke==118);
		param.a1A1_BKstar=0.392378+sigma_FFKstar_full[4]*myrand(1,method)*(method==1||ke==119);
		param.a2A1_BKstar= 1.18916+sigma_FFKstar_full[5]*myrand(1,method)*(method==1||ke==120);
		param.a0A12_BKstar=0.265375+sigma_FFKstar_full[6]*myrand(1,method)*(method==1||ke==121);
		param.a1A12_BKstar=0.533638+sigma_FFKstar_full[7]*myrand(1,method)*(method==1||ke==122);
		param.a2A12_BKstar=0.483166+sigma_FFKstar_full[8]*myrand(1,method)*(method==1||ke==123);
		param.a0V_BKstar=0.376313+sigma_FFKstar_full[9]*myrand(1,method)*(method==1||ke==124);
		param.a1V_BKstar=-1.16597+sigma_FFKstar_full[10]*myrand(1,method)*(method==1||ke==125);
		param.a2V_BKstar=2.42443+sigma_FFKstar_full[11]*myrand(1,method)*(method==1||ke==126);
		param.a0T1_BKstar=0.312055+sigma_FFKstar_full[12]*myrand(1,method)*(method==1||ke==127);
		param.a1T1_BKstar=-1.00893+sigma_FFKstar_full[13]*myrand(1,method)*(method==1||ke==128);
		param.a2T1_BKstar=1.5272+sigma_FFKstar_full[14]*myrand(1,method)*(method==1||ke==129);
		param.a0T2_BKstar=0.312055+sigma_FFKstar_full[15]*myrand(1,method)*(method==1||ke==130);
		param.a1T2_BKstar=0.496846+sigma_FFKstar_full[16]*myrand(1,method)*(method==1||ke==131);
		param.a2T2_BKstar=1.61431+sigma_FFKstar_full[17]*myrand(1,method)*(method==1||ke==132);
		param.a0T23_BKstar=0.667412+sigma_FFKstar_full[18]*myrand(1,method)*(method==1||ke==133);
		param.a1T23_BKstar=1.31812+sigma_FFKstar_full[19]*myrand(1,method)*(method==1||ke==134);
		param.a2T23_BKstar=3.82334+sigma_FFKstar_full[20]*myrand(1,method)*(method==1||ke==135);
		
		/* Form factors Bs->phi ll */	
		param.a0A0_Bsphi=0.421328+sigma_FFphi_full[0]*myrand(1,method)*(method==1||ke==136);
		param.a1A0_Bsphi=-0.976454+sigma_FFphi_full[1]*myrand(1,method)*(method==1||ke==137);
		param.a2A0_Bsphi=3.2714+sigma_FFphi_full[2]*myrand(1,method)*(method==1||ke==138);
		param.a0A1_Bsphi=0.288007+sigma_FFphi_full[3]*myrand(1,method)*(method==1||ke==139);
		param.a1A1_Bsphi=0.350826+sigma_FFphi_full[4]*myrand(1,method)*(method==1||ke==140);
		param.a2A1_Bsphi=1.69688+sigma_FFphi_full[5]*myrand(1,method)*(method==1||ke==141);
		param.a0A12_Bsphi=0.267053+sigma_FFphi_full[6]*myrand(1,method)*(method==1||ke==142);
		param.a1A12_Bsphi=0.954402+sigma_FFphi_full[7]*myrand(1,method)*(method==1||ke==143);
		param.a2A12_Bsphi=2.15263+sigma_FFphi_full[8]*myrand(1,method)*(method==1||ke==144);
		param.a0V_Bsphi=0.364478+sigma_FFphi_full[9]*myrand(1,method)*(method==1||ke==145);
		param.a1V_Bsphi=-1.22389+sigma_FFphi_full[10]*myrand(1,method)*(method==1||ke==146);
		param.a2V_Bsphi=3.74061+sigma_FFphi_full[11]*myrand(1,method)*(method==1||ke==147);
		param.a0T1_Bsphi=0.299475+sigma_FFphi_full[12]*myrand(1,method)*(method==1||ke==148);
		param.a1T1_Bsphi=-1.1013+sigma_FFphi_full[13]*myrand(1,method)*(method==1||ke==149);
		param.a2T1_Bsphi=0.58459+sigma_FFphi_full[14]*myrand(1,method)*(method==1||ke==150);
		param.a0T2_Bsphi=0.299475+sigma_FFphi_full[15]*myrand(1,method)*(method==1||ke==151);
		param.a1T2_Bsphi=0.403564+sigma_FFphi_full[16]*myrand(1,method)*(method==1||ke==152);
		param.a2T2_Bsphi=1.03987+sigma_FFphi_full[17]*myrand(1,method)*(method==1||ke==153);
		param.a0T23_Bsphi=0.65233+sigma_FFphi_full[18]*myrand(1,method)*(method==1||ke==154);
		param.a1T23_Bsphi=2.09622+sigma_FFphi_full[19]*myrand(1,method)*(method==1||ke==155);
		param.a2T23_Bsphi=6.73572+sigma_FFphi_full[20]*myrand(1,method)*(method==1||ke==156);
	}
	else
	{
		/* Form factors B->K* ll */	
		param.a0A1_BKstar=0.29725+sigma_FFKstar_soft[0]*myrand(1,method)*(method==1||ke==115);
		param.a1A1_BKstar=0.392378+sigma_FFKstar_soft[1]*myrand(1,method)*(method==1||ke==116);
		param.a2A1_BKstar= 1.18916+sigma_FFKstar_soft[2]*myrand(1,method)*(method==1||ke==117);
		param.a0A12_BKstar=0.265375+sigma_FFKstar_soft[3]*myrand(1,method)*(method==1||ke==118);
		param.a1A12_BKstar=0.533638+sigma_FFKstar_soft[4]*myrand(1,method)*(method==1||ke==119);
		param.a2A12_BKstar=0.483166+sigma_FFKstar_soft[5]*myrand(1,method)*(method==1||ke==120);
		param.a0V_BKstar=0.376313+sigma_FFKstar_soft[6]*myrand(1,method)*(method==1||ke==121);
		param.a1V_BKstar=-1.16597+sigma_FFKstar_soft[7]*myrand(1,method)*(method==1||ke==122);
		param.a2V_BKstar=2.42443+sigma_FFKstar_soft[8]*myrand(1,method)*(method==1||ke==123);
		
		/* Form factors Bs->phi ll */	
		param.a0A1_Bsphi=0.288007+sigma_FFphi_soft[0]*myrand(1,method)*(method==1||ke==124);
		param.a1A1_Bsphi=0.350826+sigma_FFphi_soft[1]*myrand(1,method)*(method==1||ke==125);
		param.a2A1_Bsphi=1.69688+sigma_FFphi_soft[2]*myrand(1,method)*(method==1||ke==126);
		param.a0A12_Bsphi=0.267053+sigma_FFphi_soft[3]*myrand(1,method)*(method==1||ke==127);
		param.a1A12_Bsphi=0.954402+sigma_FFphi_soft[4]*myrand(1,method)*(method==1||ke==128);
		param.a2A12_Bsphi=2.15263+sigma_FFphi_soft[5]*myrand(1,method)*(method==1||ke==129);
		param.a0V_Bsphi=0.364478+sigma_FFphi_soft[6]*myrand(1,method)*(method==1||ke==130);
		param.a1V_Bsphi=-1.22389+sigma_FFphi_soft[7]*myrand(1,method)*(method==1||ke==131);
		param.a2V_Bsphi=3.74061+sigma_FFphi_soft[8]*myrand(1,method)*(method==1||ke==132);
	}
	
	int test_delta0Bkstargamma=(ke<=12)||(ke>=25&&ke<=34);

	int test_BKstargamma=(ke<=12)||(ke>=25&&ke<=33);
	
	int test_bsgamma=(ke<=18);
	
	int test_Bsmumu=(ke<=12)||(ke>=84&&ke<=86);
	
	int test_BXsll=(ke<=13);
	
	int test_lowBKstar=(ke<=12)||(ke>=25&&ke<=32)||((ke>=35&&ke<=50&&(param.BKstar_implementation==1||param.BKstar_implementation==2))||(ke>=35&&ke<=43&&param.BKstar_implementation==3))||(param.fullFF&&(ke>=115&&ke<=135))||(!param.fullFF&&(ke>=115&&ke<=123));
	int test_intermBKstar=(ke<=12)||(ke>=25&&ke<=32)||((ke>=35&&ke<=50&&(param.BKstar_implementation==1||param.BKstar_implementation==2))||(ke>=35&&ke<=43&&param.BKstar_implementation==3))||(ke>=51&&ke<=58)||(param.fullFF&&(ke>=115&&ke<=135))||(!param.fullFF&&(ke>=115&&ke<=123));
	int test_highBKstar=(ke<=12)||(ke>=25&&ke<=32)||(ke>=51&&ke<=58)||(param.fullFF&&(ke>=115&&ke<=135))||(!param.fullFF&&(ke>=115&&ke<=123));

	int test_lowBK=(ke<=12)||(ke>=25&&ke<=26)||(ke>=59&&ke<=79);
	int test_intermBK=(ke<=12)||(ke>=25&&ke<=26)||(ke>=59&&ke<=83);
	int test_highBK=(ke<=12)||(ke>=25&&ke<=26)||(ke>=59&&ke<=71)||(ke>=80&&ke<=83);

	int test_lowBsphi=(ke<=12)||(ke>=84&&ke<=106)||(param.fullFF&&(ke>=136&&ke<=156))||(!param.fullFF&&(ke>=124&&ke<=132));
	int test_intermBsphi=(ke<=12)||(ke>=84&&ke<=114)||(param.fullFF&&(ke>=136&&ke<=156))||(!param.fullFF&&(ke>=124&&ke<=132));
	int test_highBsphi=(ke<=12)||(ke>=84&&ke<=98)||(ke>=107&&ke<=114)||(param.fullFF&&(ke>=136&&ke<=156))||(!param.fullFF&&(ke>=124&&ke<=132));

	CW_calculator(2,C0w,C1w,C2w,mu_W,&param); /* 2 = muon */
	CW_calculator(1,C0ew,C1ew,C2ew,mu_W,&param); /* 1 = electron */

	if(onoff[nbobsall]&&(test_delta0Bkstargamma))
	{
		C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
		C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
		Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);		
		
		if(deltaC!=NULL) for(ie=1;ie<=8;ie++) C0b[ie]+=deltaC[ie];
		if(deltaCp!=NULL) for(ie=1;ie<=8;ie++) Cpb[ie]+=deltaCp[ie];
		
		values[nbobs]=delta0(C0b,C0spec,C1b,C1spec,Cpb,&param,mu_b,mu_spec,lambda_h);
	}
	if(onoff[nbobsall++]) names[nbobs++]="delta0";
	
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);

	C_calculator_base1(C0ew,C1ew,C2ew,mu_W,C0eb,C1eb,C2eb,mu_b,&param);
	Cprime_calculator(1,Cpeb,CQpeb,mu_W,mu_b,&param);
	CQ_calculator(1,CQ0eb,CQ1eb,mu_W,mu_b,&param);

	if(deltaC!=NULL) for(ie=1;ie<=10;ie++)
	{
		if(ie>=9)
		{
			C0b[ie]+=deltaC[ie];
			C0eb[ie]+=deltaC[ie+10];
		}
		else
		{
			C0b[ie]+=deltaC[ie];
			C0eb[ie]+=deltaC[ie];
		}
	}
	if(deltaCp!=NULL) for(ie=1;ie<=10;ie++)
	{
		if(ie>=9)
		{
			Cpb[ie]+=deltaCp[ie];
			Cpeb[ie]+=deltaCp[ie+10];
		}
		else
		{
			Cpb[ie]+=deltaCp[ie];
			Cpeb[ie]+=deltaCp[ie];
		}
	}
	if(deltaCQ!=NULL) for(ie=1;ie<=2;ie++)
	{
		CQ0b[ie]+=deltaCQ[ie];
		CQ0eb[ie]+=deltaCQ[ie+2];
	}
	if(deltaCQp!=NULL) for(ie=1;ie<=2;ie++)
	{
		CQpb[ie]+=deltaCQp[ie];
		CQpeb[ie]+=deltaCQp[ie+2];
	}

	if(onoff[nbobsall]&&(test_bsgamma)) values[nbobs]=bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,&param)*(1.+bsgamma_rand);
	if(onoff[nbobsall++]) names[nbobs++]="bsgamma";
	
	if(onoff[nbobsall]&&(test_Bsmumu)) values[nbobs]=Bsmumu_untag(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
	if(onoff[nbobsall++]) names[nbobs++]="Bsmumu_untag";
	if(onoff[nbobsall]&&(test_Bsmumu)) values[nbobs]=Bsll_untag(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
	if(onoff[nbobsall++]) names[nbobs++]="Bsee_untag";
	if(onoff[nbobsall]&&(test_Bsmumu)) values[nbobs]=Bdmumu(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
	if(onoff[nbobsall++]) names[nbobs++]="Bdmumu";

	if(onoff[nbobsall]&&(test_BXsll||(ke==19))) values[nbobs]=BRBXsll_lowq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)*(1.+BRBXsmumu_lowq2_rand);
	if(onoff[nbobsall++]) names[nbobs++]="BRBXsmumu_lowq2";
	if(onoff[nbobsall]&&(test_BXsll||(ke==20))) values[nbobs]=BRBXsll_highq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)*(1.+BRBXsmumu_highq2_rand);
	if(onoff[nbobsall++]) names[nbobs++]="BRBXsmumu_highq2";
	if(onoff[nbobsall]&&(test_BXsll||(ke==21))) values[nbobs]=BRBXsll(2,0.045/param.mass_b_1S/param.mass_b_1S,0.999,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)*(1.+BRBXsmumu_full_rand);
	if(onoff[nbobsall++]) names[nbobs++]="BRBXsmumu_full";

	if(onoff[nbobsall]&&(test_BXsll||(ke==22))) values[nbobs]=BRBXsll_lowq2(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)*(1.+BRBXsee_lowq2_rand);
	if(onoff[nbobsall++]) names[nbobs++]="BRBXsee_lowq2";
	if(onoff[nbobsall]&&(test_BXsll||(ke==23))) values[nbobs]=BRBXsll_highq2(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)*(1.+BRBXsee_highq2_rand);
	if(onoff[nbobsall++]) names[nbobs++]="BRBXsee_highq2";
	if(onoff[nbobsall]&&(test_BXsll||(ke==24))) values[nbobs]=BRBXsll(1,0.045/param.mass_b_1S/param.mass_b_1S,0.999,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)*(1.+BRBXsee_full_rand);
	if(onoff[nbobsall++]) names[nbobs++]="BRBXsee_full";
	
	if(onoff[nbobsall]&&(test_BKstargamma)) values[nbobs]=BR_BKstargamma(0,C0b,C1b,C2b,Cpb,&param,mu_b);
	if(onoff[nbobsall++]) names[nbobs++]="BR_B0Kstar0gamma";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,0.1,0.98,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(0.98-0.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_0.1_0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	if(onoff[nbobsall++]) names[nbobs++]="ATRe_B0Kstar0mumu_0.1_0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_0.1-0.98";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,1.1,2.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.5-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_1.1-2.5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_1.1-2.5";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,2.5,4,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4-2.5);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_2.5-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_2.5-4";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,4.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_4-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_4-6";

	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,0,6.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_6-8";

	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,0,11.,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_11-12.5";

	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,15.,17.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_15-17";

	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,17.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-17.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_17-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_17-19";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	if(onoff[nbobsall++]) names[nbobs++]="ATRe_B0Kstar0mumu_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_1.1_6";
	
	
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,15.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_15-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_15-19";
	

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,1.1,2.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_1.1-2";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,2.,3.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(3.-2.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_2-3";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,3.,4.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-3.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_3-4";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,4.,5.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(5.-4.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_4-5";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,5.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-5.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_5-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_5-6";

	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,0,6.,7.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(7.-6.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_6-7";

	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,0,7.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-7.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_7-8";

	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,0,11.,11.75,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(11.75-11.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_11-11.75";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_11-11.75";

	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,0,11.75,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.75);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_11.75-12.5";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_11.75-12.5";

	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,15.,16.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(16.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_15-16";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_15-16";

	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,16,17,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17-16);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_16-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_16-17";

	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,17.,18.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(18.-17.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_17-18";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_17-18";
	
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,0,18.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-18.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[1];
	if(onoff[nbobsall++]) names[nbobs++]="AFB_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[25];
	if(onoff[nbobsall++]) names[nbobs++]="S3_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[26];
	if(onoff[nbobsall++]) names[nbobs++]="S4_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[27];
	if(onoff[nbobsall++]) names[nbobs++]="S5_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[28];
	if(onoff[nbobsall++]) names[nbobs++]="S7_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[29];
	if(onoff[nbobsall++]) names[nbobs++]="S8_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[30];
	if(onoff[nbobsall++]) names[nbobs++]="S9_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[18];
	if(onoff[nbobsall++]) names[nbobs++]="P5p_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[14];
	if(onoff[nbobsall++]) names[nbobs++]="P2_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[15];
	if(onoff[nbobsall++]) names[nbobs++]="P3_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[17];
	if(onoff[nbobsall++]) names[nbobs++]="P4p_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[19];
	if(onoff[nbobsall++]) names[nbobs++]="P6p_B0Kstar0mumu_18-19";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=obsKstar[21];
	if(onoff[nbobsall++]) names[nbobs++]="P8p_B0Kstar0mumu_18-19";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,1,0.1,2.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-0.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_0.1-2";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,1,2.,4.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-2.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_2-4";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,1,4.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_4-6";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,1,6.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(2,1,11.,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_11-12.5";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,1,15.,17.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_15-17";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,1,17.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-17.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_17-19";

	
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,1,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_1.1-6";
	if(onoff[nbobsall]&&(test_highBKstar)) values[nbobs]=BRBKstarll(2,1,15.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+Kstar+mumu_15-19";
	
	if(onoff[nbobsall]&&(test_intermBKstar)) values[nbobs]=BRBKstarll(1,0,0.1,pow(param.m_Bd-param.m_Kstar0,2.)*0.999,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
	if(onoff[nbobsall++]) names[nbobs++]="BR_B0Kstar0ee_full";


	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(1,0,0.0009,1.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(1.-0.0009);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0ee_0.0009-1";
			
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(1,0,0.002,1.12,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(1.12-0.002);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0ee_0.002-1.12";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0ee_0.002_1.12";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="AT2_B0Kstar0ee_0.002_1.12";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	if(onoff[nbobsall++]) names[nbobs++]="ATRe_B0Kstar0ee_0.002_1.12";
			
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(1,0,1.1,6.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(6.-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0ee_1.1-6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[2];
	if(onoff[nbobsall++]) names[nbobs++]="FL_B0Kstar0ee_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=obsKstar[5];
	if(onoff[nbobsall++]) names[nbobs++]="P1_B0Kstar0ee_1.1_6";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	if(onoff[nbobsall++]) names[nbobs++]="ATRe_B0Kstar0ee_1.1_6";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(1,0,1.,6.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(6.-1.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0Kstar0ee_1-6";

	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,0.045,1.1,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/BRBKstarll(1,0,0.045,1.1,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)-1.;
	if(onoff[nbobsall++]) names[nbobs++]="RKstar_0.045-1.1-1";
	if(onoff[nbobsall]&&(test_lowBKstar)) values[nbobs]=BRBKstarll(2,0,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/BRBKstarll(1,0,1.1,6.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)-1.;
	if(onoff[nbobsall++]) names[nbobs++]="RKstar_1.1-6-1";

	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,0,0.1,2.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-0.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_0.1-2";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,0,2.,4.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-2.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_2-4";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,0,4.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_4-6";
	if(onoff[nbobsall]&&(test_intermBK)) values[nbobs]=BRBKll(2,0,6.,8.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_6-8";
	if(onoff[nbobsall]&&(test_intermBK)) values[nbobs]=BRBKll(2,0,11.,12.5,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_11-12.5";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,0,15.,17.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_15-17";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,0,17.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-17.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_17-22";
	
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,0,1.1,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_1.1-6";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,0,15.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B0K0mumu_15-22";


	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,0.1,0.98,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(0.98-0.1);		
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_0.1-0.98";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,1.1,2.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-1.1);		
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_1.1-2";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,2.,3.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(3.-2.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_2-3";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,3,4.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-3.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_3-4";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,4.,5.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(5.-4.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_4-5";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,5.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-5.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_5-6";
	if(onoff[nbobsall]&&(test_intermBK)) values[nbobs]=BRBKll(2,1,6.,7.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(7.-6.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_6-7";
	if(onoff[nbobsall]&&(test_intermBK)) values[nbobs]=BRBKll(2,1,7.,8.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-7.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_7-8";
	if(onoff[nbobsall]&&(test_intermBK)) values[nbobs]=BRBKll(2,1,11.,11.8,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(11.8-11.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_11-11.8";
	if(onoff[nbobsall]&&(test_intermBK)) values[nbobs]=BRBKll(2,1,11.8,12.5,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.8);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_11.8-12.5";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,15.,16.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(16.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_15-16";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,16.,17.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-16.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_16-17";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,17.,18.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(18.-17.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_17-18";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,18.,19.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-18.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_18-19";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,19.,20.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(20.-19.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_19-20";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,20.,21.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(21.-20.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_20-21";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,21.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-21.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_21-22";

	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,1.1,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_1.1-6";
	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=obsK[2];
	if(onoff[nbobsall++]) names[nbobs++]="FH_B+K+mumu_1.1-6";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=BRBKll(2,1,15.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_15-22";
	if(onoff[nbobsall]&&(test_highBK)) values[nbobs]=obsK[2];
	if(onoff[nbobsall++]) names[nbobs++]="FH_B+K+mumu_15-22";

	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,1.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+mumu_1-6";

	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(1,1,1.,6.,obsK,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(6.-1.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_B+K+ee_1-6";

	if(onoff[nbobsall]&&(test_lowBK)) values[nbobs]=BRBKll(2,1,1.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/BRBKll(1,1,1.,6.,obsK,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)-1.;
	if(onoff[nbobsall++]) names[nbobs++]="RK_1-6-1";

	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=BRBsphill(2,0.1,2.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-0.1);
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_0.1-2";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_0.1-2";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_0.1-2";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_0.1-2";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_0.1-2";

	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=BRBsphill(2,2.,5.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(5.-2.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_2-5";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_2-5";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_2-5";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_2-5";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_2-5";

	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=BRBsphill(2,5.,8.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-5.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_5-8";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_5-8";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_5-8";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_5-8";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_5-8";

	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=BRBsphill(2,11.,12.5,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_11-12.5";
	if(onoff[nbobsall]&&(test_intermBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_11-12.5";

	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=BRBsphill(2,15.,17.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_15-17";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_15-17";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_15-17";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_15-17";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_15-17";

	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=BRBsphill(2,17.,pow(param.m_Bs-param.m_phi,2.)*0.999,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(pow(param.m_Bs-param.m_phi,2.)*0.999-17.); // endpoint problem!
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_17-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_17-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_17-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_17-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_17-19";

	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=BRBsphill(2,1.,6.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.);
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_1-6";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_1-6";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_1-6";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_1-6";
	if(onoff[nbobsall]&&(test_lowBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_1-6";

	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=BRBsphill(2,15.,pow(param.m_Bs-param.m_phi,2.)*0.999,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(pow(param.m_Bs-param.m_phi,2.)*0.999-15.); // endpoint problem!
	if(onoff[nbobsall++]) names[nbobs++]="dG_Bsphimumu_15-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[1];
	if(onoff[nbobsall++]) names[nbobs++]="FL_Bsphimumu_15-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[3];
	if(onoff[nbobsall++]) names[nbobs++]="S3_Bsphimumu_15-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[4];
	if(onoff[nbobsall++]) names[nbobs++]="S4_Bsphimumu_15-19";
	if(onoff[nbobsall]&&(test_highBsphi)) values[nbobs]=obsphi[6];
	if(onoff[nbobsall++]) names[nbobs++]="S7_Bsphimumu_15-19";
	
	return nbobs;
}

/*---------------------------------------------------------------------*/

void activate_observables(int nbobs0, char* names0[], int nbobsin, char** namesin, int onoff[])
{	
	int ie,je;
	for(ie=0;ie<NBOBSMAX;ie++) onoff[ie]=0;
	
	for(ie=0;ie<nbobsin;ie++)
	{
		je=0;
		while(strcmp(namesin[ie],names0[je])&&je<nbobs0) je++;
		if(je<nbobs0) {if(!strcmp(namesin[ie],names0[je])) onoff[je]=1;}
		else printf("%s experimental measurement not available, removed from fit!\n",namesin[ie]);
	}
	
	return;
}

/*---------------------------------------------------------------------*/

void correlation_matrices(int fullFF, double **correlation_matrix_FFKstar, double **correlation_matrix_FFphi, double **correlation_matrix_FFK)
{
	int ie,je;
	
	double correlation_matrix_FFKstar_full[21][21];
	/* correlation matrix for B-> Kstar ll from 1503.95534 - updated to v2 */
	for(ie=0;ie<21;ie++) for(je=0;je<21;je++) correlation_matrix_FFKstar_full[ie][je]=(ie==je);
	correlation_matrix_FFKstar_full[0][0]=1.0;
	correlation_matrix_FFKstar_full[0][1]=0.633689;
	correlation_matrix_FFKstar_full[0][2]=0.0575305;
	correlation_matrix_FFKstar_full[1][1]=1.0;
	correlation_matrix_FFKstar_full[1][2]=0.486314;
	correlation_matrix_FFKstar_full[2][2]=1.0;
	correlation_matrix_FFKstar_full[0][3]=0.0615373;
	correlation_matrix_FFKstar_full[0][4]=0.169044;
	correlation_matrix_FFKstar_full[0][5]=0.240647;
	correlation_matrix_FFKstar_full[1][3]=0.300568;
	correlation_matrix_FFKstar_full[1][4]=0.472167;
	correlation_matrix_FFKstar_full[1][5]=0.321569;
	correlation_matrix_FFKstar_full[2][3]=0.515475;
	correlation_matrix_FFKstar_full[2][4]=0.743538;
	correlation_matrix_FFKstar_full[2][5]=0.661689;
	correlation_matrix_FFKstar_full[0][6]=0.999; /* was 1 */
	correlation_matrix_FFKstar_full[0][7]=0.609495;
	correlation_matrix_FFKstar_full[0][8]=0.0355029;
	correlation_matrix_FFKstar_full[1][6]=0.621073; 
	correlation_matrix_FFKstar_full[1][7]=0.756564;
	correlation_matrix_FFKstar_full[1][8]=0.46913;
	correlation_matrix_FFKstar_full[2][6]=0.05503;
	correlation_matrix_FFKstar_full[2][7]=0.326536;
	correlation_matrix_FFKstar_full[2][8]=0.665149;
	correlation_matrix_FFKstar_full[0][9]=0.0393316;
	correlation_matrix_FFKstar_full[0][10]=0.178168;
	correlation_matrix_FFKstar_full[0][11]=0.195421;
	correlation_matrix_FFKstar_full[1][9]=0.280475;
	correlation_matrix_FFKstar_full[1][10]=0.553833;
	correlation_matrix_FFKstar_full[1][11]=0.0752289;
	correlation_matrix_FFKstar_full[2][9]=0.545574;
	correlation_matrix_FFKstar_full[2][10]=0.731129;
	correlation_matrix_FFKstar_full[2][11]=0.00860747;
	correlation_matrix_FFKstar_full[0][12]=0.0164883;
	correlation_matrix_FFKstar_full[0][13]=0.0452644;
	correlation_matrix_FFKstar_full[0][14]=-0.158457;
	correlation_matrix_FFKstar_full[1][12]=0.241261;
	correlation_matrix_FFKstar_full[1][13]=0.448732;
	correlation_matrix_FFKstar_full[1][14]=-0.182566;
	correlation_matrix_FFKstar_full[2][12]=0.475401;
	correlation_matrix_FFKstar_full[2][13]=0.674459;
	correlation_matrix_FFKstar_full[2][14]=-0.125613;
	correlation_matrix_FFKstar_full[0][15]=0.0217598;
	correlation_matrix_FFKstar_full[0][16]=0.150264;
	correlation_matrix_FFKstar_full[0][17]=0.0311554;
	correlation_matrix_FFKstar_full[1][15]=0.24405;
	correlation_matrix_FFKstar_full[1][16]=0.495428;
	correlation_matrix_FFKstar_full[1][17]=0.181197;
	correlation_matrix_FFKstar_full[2][15]=0.474921;
	correlation_matrix_FFKstar_full[2][16]=0.701034;
	correlation_matrix_FFKstar_full[2][17]=0.273928;
	correlation_matrix_FFKstar_full[0][18]=0.451075;
	correlation_matrix_FFKstar_full[0][19]=0.104911;
	correlation_matrix_FFKstar_full[0][20]=-0.347452;
	correlation_matrix_FFKstar_full[1][18]=0.573536;
	correlation_matrix_FFKstar_full[1][19]=0.53356;
	correlation_matrix_FFKstar_full[1][20]=-0.19402;
	correlation_matrix_FFKstar_full[2][18]=0.663092;
	correlation_matrix_FFKstar_full[2][19]=0.426721;
	correlation_matrix_FFKstar_full[2][20]=-0.318126;
	correlation_matrix_FFKstar_full[3][3]=1.0;
	correlation_matrix_FFKstar_full[3][4]=0.702494;
	correlation_matrix_FFKstar_full[3][5]=0.153604;
	correlation_matrix_FFKstar_full[4][4]=1.0;
	correlation_matrix_FFKstar_full[4][5]=0.747682;
	correlation_matrix_FFKstar_full[5][5]=1.0;
	correlation_matrix_FFKstar_full[3][6]=0.0555181;
	correlation_matrix_FFKstar_full[3][7]=0.0105278;
	correlation_matrix_FFKstar_full[3][8]=0.205325;
	correlation_matrix_FFKstar_full[4][6]=0.166407;
	correlation_matrix_FFKstar_full[4][7]=0.284455;
	correlation_matrix_FFKstar_full[4][8]=0.480939;
	correlation_matrix_FFKstar_full[5][6]=0.245093;
	correlation_matrix_FFKstar_full[5][7]=0.316222;
	correlation_matrix_FFKstar_full[5][8]=0.535592;
	correlation_matrix_FFKstar_full[3][9]=0.922982;
	correlation_matrix_FFKstar_full[3][10]=0.731071;
	correlation_matrix_FFKstar_full[3][11]=-0.35833;
	correlation_matrix_FFKstar_full[4][9]=0.678427;
	correlation_matrix_FFKstar_full[4][10]=0.862475;
	correlation_matrix_FFKstar_full[4][11]=-0.0633373;
	correlation_matrix_FFKstar_full[5][9]=0.194572;
	correlation_matrix_FFKstar_full[5][10]=0.517205;
	correlation_matrix_FFKstar_full[5][11]=0.253395;
	correlation_matrix_FFKstar_full[3][12]=0.899902;
	correlation_matrix_FFKstar_full[3][13]=0.68101;
	correlation_matrix_FFKstar_full[3][14]=-0.468536;
	correlation_matrix_FFKstar_full[4][12]=0.630965;
	correlation_matrix_FFKstar_full[4][13]=0.760031;
	correlation_matrix_FFKstar_full[4][14]=-0.373063;
	correlation_matrix_FFKstar_full[5][12]=0.135433;
	correlation_matrix_FFKstar_full[5][13]=0.400555;
	correlation_matrix_FFKstar_full[5][14]=-0.127614;
	correlation_matrix_FFKstar_full[3][15]=0.899616;
	correlation_matrix_FFKstar_full[3][16]=0.692637;
	correlation_matrix_FFKstar_full[3][17]=-0.202092;
	correlation_matrix_FFKstar_full[4][15]=0.634299;
	correlation_matrix_FFKstar_full[4][16]=0.915195;
	correlation_matrix_FFKstar_full[4][17]=0.226359;
	correlation_matrix_FFKstar_full[5][15]=0.140073;
	correlation_matrix_FFKstar_full[5][16]=0.59715;
	correlation_matrix_FFKstar_full[5][17]=0.402339;
	correlation_matrix_FFKstar_full[3][18]=0.564524;
	correlation_matrix_FFKstar_full[3][19]=0.30614;
	correlation_matrix_FFKstar_full[3][20]=-0.362387;
	correlation_matrix_FFKstar_full[4][18]=0.695868;
	correlation_matrix_FFKstar_full[4][19]=0.457217;
	correlation_matrix_FFKstar_full[4][20]=-0.440385;
	correlation_matrix_FFKstar_full[5][18]=0.553618;
	correlation_matrix_FFKstar_full[5][19]=0.252273;
	correlation_matrix_FFKstar_full[5][20]=-0.40495;
	correlation_matrix_FFKstar_full[6][6]=1.0;
	correlation_matrix_FFKstar_full[6][7]=0.617726;
	correlation_matrix_FFKstar_full[6][8]=0.0443495;
	correlation_matrix_FFKstar_full[7][7]=1.0;
	correlation_matrix_FFKstar_full[7][8]=0.700096;
	correlation_matrix_FFKstar_full[8][8]=1.0;
	correlation_matrix_FFKstar_full[6][9]=0.0335119;
	correlation_matrix_FFKstar_full[6][10]=0.172545;
	correlation_matrix_FFKstar_full[6][11]=0.202654;
	correlation_matrix_FFKstar_full[7][9]=-0.0140921;
	correlation_matrix_FFKstar_full[7][10]=0.350603;
	correlation_matrix_FFKstar_full[7][11]=0.234661;
	correlation_matrix_FFKstar_full[8][9]=0.196602;
	correlation_matrix_FFKstar_full[8][10]=0.485841;
	correlation_matrix_FFKstar_full[8][11]=0.139392;
	correlation_matrix_FFKstar_full[6][12]=0.010524;
	correlation_matrix_FFKstar_full[6][13]=0.0391092;
	correlation_matrix_FFKstar_full[6][14]=-0.154503;
	correlation_matrix_FFKstar_full[7][12]=-0.0243359;
	correlation_matrix_FFKstar_full[7][13]=0.234436;
	correlation_matrix_FFKstar_full[7][14]=-0.0609604;
	correlation_matrix_FFKstar_full[8][12]=0.173373;
	correlation_matrix_FFKstar_full[8][13]=0.402833;
	correlation_matrix_FFKstar_full[8][14]=-0.0779755;
	correlation_matrix_FFKstar_full[6][15]=0.0157345;
	correlation_matrix_FFKstar_full[6][16]=0.145566;
	correlation_matrix_FFKstar_full[6][17]=0.0343108;
	correlation_matrix_FFKstar_full[7][15]=-0.0214341;
	correlation_matrix_FFKstar_full[7][16]=0.28577;
	correlation_matrix_FFKstar_full[7][17]=0.248484;
	correlation_matrix_FFKstar_full[8][15]=0.173504;
	correlation_matrix_FFKstar_full[8][16]=0.422339;
	correlation_matrix_FFKstar_full[8][17]=0.253397;
	correlation_matrix_FFKstar_full[6][18]=0.448529;
	correlation_matrix_FFKstar_full[6][19]=0.100428;
	correlation_matrix_FFKstar_full[6][20]=-0.348189;
	correlation_matrix_FFKstar_full[7][18]=0.386552;
	correlation_matrix_FFKstar_full[7][19]=0.537224;
	correlation_matrix_FFKstar_full[7][20]=-0.0377692;
	correlation_matrix_FFKstar_full[8][18]=0.403006;
	correlation_matrix_FFKstar_full[8][19]=0.571118;
	correlation_matrix_FFKstar_full[8][20]=-0.0276853;
	correlation_matrix_FFKstar_full[9][9]=1.0;
	correlation_matrix_FFKstar_full[9][10]=0.757379;
	correlation_matrix_FFKstar_full[9][11]=-0.397005;
	correlation_matrix_FFKstar_full[10][10]=1.0;
	correlation_matrix_FFKstar_full[10][11]=0.0346143;
	correlation_matrix_FFKstar_full[11][11]=1.0;
	correlation_matrix_FFKstar_full[9][12]=0.901557;
	correlation_matrix_FFKstar_full[9][13]=0.716437;
	correlation_matrix_FFKstar_full[9][14]=-0.400738;
	correlation_matrix_FFKstar_full[10][12]=0.709376;
	correlation_matrix_FFKstar_full[10][13]=0.906137;
	correlation_matrix_FFKstar_full[10][14]=-0.236675;
	correlation_matrix_FFKstar_full[11][12]=-0.354268;
	correlation_matrix_FFKstar_full[11][13]=-0.0783557;
	correlation_matrix_FFKstar_full[11][14]=0.338328;
	correlation_matrix_FFKstar_full[9][15]=0.899723;
	correlation_matrix_FFKstar_full[9][16]=0.702845;
	correlation_matrix_FFKstar_full[9][17]=-0.147198;
	correlation_matrix_FFKstar_full[10][15]=0.709471;
	correlation_matrix_FFKstar_full[10][16]=0.916209;
	correlation_matrix_FFKstar_full[10][17]=0.261192;
	correlation_matrix_FFKstar_full[11][15]=-0.352355;
	correlation_matrix_FFKstar_full[11][16]=-0.0541063;
	correlation_matrix_FFKstar_full[11][17]=0.38493;
	correlation_matrix_FFKstar_full[9][18]=0.573778;
	correlation_matrix_FFKstar_full[9][19]=0.277866;
	correlation_matrix_FFKstar_full[9][20]=-0.360324;
	correlation_matrix_FFKstar_full[10][18]=0.70266;
	correlation_matrix_FFKstar_full[10][19]=0.480664;
	correlation_matrix_FFKstar_full[10][20]=-0.355976;
	correlation_matrix_FFKstar_full[11][18]=0.016202;
	correlation_matrix_FFKstar_full[11][19]=6.71602e-5;
	correlation_matrix_FFKstar_full[11][20]=0.0569126;
	correlation_matrix_FFKstar_full[12][12]=1.0;
	correlation_matrix_FFKstar_full[12][13]=0.713877;
	correlation_matrix_FFKstar_full[12][14]=-0.473872;
	correlation_matrix_FFKstar_full[13][13]=1.0;
	correlation_matrix_FFKstar_full[13][14]=0.00093065;
	correlation_matrix_FFKstar_full[14][14]=1.0;
	correlation_matrix_FFKstar_full[12][15]=0.999; /* was 1 */
	correlation_matrix_FFKstar_full[12][16]=0.686349;
	correlation_matrix_FFKstar_full[12][17]=-0.229853;
	correlation_matrix_FFKstar_full[13][15]=0.706907;
	correlation_matrix_FFKstar_full[13][16]=0.881232;
	correlation_matrix_FFKstar_full[13][17]=0.365115;
	correlation_matrix_FFKstar_full[14][15]=-0.48143;
	correlation_matrix_FFKstar_full[14][16]=-0.239308;
	correlation_matrix_FFKstar_full[14][17]=0.626511;
	correlation_matrix_FFKstar_full[12][18]=0.515437;
	correlation_matrix_FFKstar_full[12][19]=0.236189;
	correlation_matrix_FFKstar_full[12][20]=-0.310316;
	correlation_matrix_FFKstar_full[13][18]=0.601021;
	correlation_matrix_FFKstar_full[13][19]=0.427022;
	correlation_matrix_FFKstar_full[13][20]=-0.210129;
	correlation_matrix_FFKstar_full[14][18]=-0.27088;
	correlation_matrix_FFKstar_full[14][19]=-0.105814;
	correlation_matrix_FFKstar_full[14][20]=0.447129;
	correlation_matrix_FFKstar_full[15][15]=1.0;
	correlation_matrix_FFKstar_full[15][16]=0.692809;
	correlation_matrix_FFKstar_full[15][17]=-0.227854;
	correlation_matrix_FFKstar_full[16][16]=1.0;
	correlation_matrix_FFKstar_full[16][17]=0.382807;
	correlation_matrix_FFKstar_full[17][17]=1.0;
	correlation_matrix_FFKstar_full[15][18]=0.517204;
	correlation_matrix_FFKstar_full[15][19]=0.237634;
	correlation_matrix_FFKstar_full[15][20]=-0.313381;
	correlation_matrix_FFKstar_full[16][18]=0.658722;
	correlation_matrix_FFKstar_full[16][19]=0.46875;
	correlation_matrix_FFKstar_full[16][20]=-0.326765;
	correlation_matrix_FFKstar_full[17][18]=0.0969361;
	correlation_matrix_FFKstar_full[17][19]=0.216438;
	correlation_matrix_FFKstar_full[17][20]=0.27039;
	correlation_matrix_FFKstar_full[18][18]=1.0;
	correlation_matrix_FFKstar_full[18][19]=0.327243;
	correlation_matrix_FFKstar_full[18][20]=-0.711856;
	correlation_matrix_FFKstar_full[19][19]=1.0;
	correlation_matrix_FFKstar_full[19][20]=0.287454;
	correlation_matrix_FFKstar_full[20][20]=1.0;
    for(ie=0;ie<21;ie++) for(je=ie+1;je<21;je++) correlation_matrix_FFKstar_full[je][ie]=correlation_matrix_FFKstar_full[ie][je];
    
    double correlation_matrix_FFKstar_soft[9][9];
	for(ie=0;ie<9;ie++) for(je=0;je<9;je++) correlation_matrix_FFKstar_soft[ie][je]=(ie==je);
	correlation_matrix_FFKstar_soft[3-3][3-3]=correlation_matrix_FFKstar_full[3][3];
	correlation_matrix_FFKstar_soft[3-3][4-3]=correlation_matrix_FFKstar_full[3][4];
	correlation_matrix_FFKstar_soft[3-3][5-3]=correlation_matrix_FFKstar_full[3][5];
	correlation_matrix_FFKstar_soft[4-3][4-3]=correlation_matrix_FFKstar_full[4][4];
	correlation_matrix_FFKstar_soft[4-3][5-3]=correlation_matrix_FFKstar_full[4][5];
	correlation_matrix_FFKstar_soft[5-3][5-3]=correlation_matrix_FFKstar_full[5][5];
	correlation_matrix_FFKstar_soft[3-3][6-3]=correlation_matrix_FFKstar_full[3][6];
	correlation_matrix_FFKstar_soft[3-3][7-3]=correlation_matrix_FFKstar_full[3][7];
	correlation_matrix_FFKstar_soft[3-3][8-3]=correlation_matrix_FFKstar_full[3][8];
	correlation_matrix_FFKstar_soft[4-3][6-3]=correlation_matrix_FFKstar_full[4][6];
	correlation_matrix_FFKstar_soft[4-3][7-3]=correlation_matrix_FFKstar_full[4][7];
	correlation_matrix_FFKstar_soft[4-3][8-3]=correlation_matrix_FFKstar_full[4][8];
	correlation_matrix_FFKstar_soft[5-3][6-3]=correlation_matrix_FFKstar_full[5][6];
	correlation_matrix_FFKstar_soft[5-3][7-3]=correlation_matrix_FFKstar_full[5][7];
	correlation_matrix_FFKstar_soft[5-3][8-3]=correlation_matrix_FFKstar_full[5][8];
	correlation_matrix_FFKstar_soft[3-3][9-3]=correlation_matrix_FFKstar_full[3][9];
	correlation_matrix_FFKstar_soft[3-3][10-3]=correlation_matrix_FFKstar_full[3][10];
	correlation_matrix_FFKstar_soft[3-3][11-3]=correlation_matrix_FFKstar_full[3][11];
	correlation_matrix_FFKstar_soft[4-3][9-3]=correlation_matrix_FFKstar_full[4][9];
	correlation_matrix_FFKstar_soft[4-3][10-3]=correlation_matrix_FFKstar_full[4][10];
	correlation_matrix_FFKstar_soft[4-3][11-3]=correlation_matrix_FFKstar_full[4][11];
	correlation_matrix_FFKstar_soft[5-3][9-3]=correlation_matrix_FFKstar_full[5][9];
	correlation_matrix_FFKstar_soft[5-3][10-3]=correlation_matrix_FFKstar_full[5][10];
	correlation_matrix_FFKstar_soft[5-3][11-3]=correlation_matrix_FFKstar_full[5][11];
	correlation_matrix_FFKstar_soft[6-3][6-3]=correlation_matrix_FFKstar_full[6][6];
	correlation_matrix_FFKstar_soft[6-3][7-3]=correlation_matrix_FFKstar_full[6][7];
	correlation_matrix_FFKstar_soft[6-3][8-3]=correlation_matrix_FFKstar_full[6][8];
	correlation_matrix_FFKstar_soft[7-3][7-3]=correlation_matrix_FFKstar_full[7][7];
	correlation_matrix_FFKstar_soft[7-3][8-3]=correlation_matrix_FFKstar_full[7][8];
	correlation_matrix_FFKstar_soft[8-3][8-3]=correlation_matrix_FFKstar_full[8][8];
	correlation_matrix_FFKstar_soft[6-3][9-3]=correlation_matrix_FFKstar_full[6][9];
	correlation_matrix_FFKstar_soft[6-3][10-3]=correlation_matrix_FFKstar_full[6][10];
	correlation_matrix_FFKstar_soft[6-3][11-3]=correlation_matrix_FFKstar_full[6][11];
	correlation_matrix_FFKstar_soft[7-3][9-3]=correlation_matrix_FFKstar_full[7][9];
	correlation_matrix_FFKstar_soft[7-3][10-3]=correlation_matrix_FFKstar_full[7][10];
	correlation_matrix_FFKstar_soft[7-3][11-3]=correlation_matrix_FFKstar_full[7][11];
	correlation_matrix_FFKstar_soft[8-3][9-3]=correlation_matrix_FFKstar_full[8][9];
	correlation_matrix_FFKstar_soft[8-3][10-3]=correlation_matrix_FFKstar_full[8][10];
	correlation_matrix_FFKstar_soft[8-3][11-3]=correlation_matrix_FFKstar_full[8][11];
	correlation_matrix_FFKstar_soft[9-3][9-3]=correlation_matrix_FFKstar_full[9][9];
	correlation_matrix_FFKstar_soft[9-3][10-3]=correlation_matrix_FFKstar_full[9][10];
	correlation_matrix_FFKstar_soft[9-3][11-3]=correlation_matrix_FFKstar_full[9][11];
	correlation_matrix_FFKstar_soft[10-3][10-3]=correlation_matrix_FFKstar_full[10][10];
	correlation_matrix_FFKstar_soft[10-3][11-3]=correlation_matrix_FFKstar_full[10][11];
	correlation_matrix_FFKstar_soft[11-3][11-3]=correlation_matrix_FFKstar_full[11][11];

    for(ie=0;ie<9;ie++) for(je=ie+1;je<9;je++) correlation_matrix_FFKstar_soft[je][ie]=correlation_matrix_FFKstar_soft[ie][je];
    
    
    /* correlation matrix for Bs-> phi ll from 1503.95534 - updated to v2 */
    double correlation_matrix_FFphi_full[21][21];
	for(ie=0;ie<21;ie++) for(je=0;je<21;je++) correlation_matrix_FFphi_full[ie][je]=(ie==je);
	correlation_matrix_FFphi_full[0][0]=1.0;
	correlation_matrix_FFphi_full[0][1]=0.687348;
	correlation_matrix_FFphi_full[0][2]=-0.349657;
	correlation_matrix_FFphi_full[1][1]=1.0;
	correlation_matrix_FFphi_full[1][2]=0.121069;
	correlation_matrix_FFphi_full[2][2]=1.0;
	correlation_matrix_FFphi_full[0][3]=0.206736;
	correlation_matrix_FFphi_full[0][4]=-0.554941;
	correlation_matrix_FFphi_full[0][5]=-0.597622;
	correlation_matrix_FFphi_full[1][3]=0.0064726;
	correlation_matrix_FFphi_full[1][4]=-0.609677;
	correlation_matrix_FFphi_full[1][5]=-0.668983;
	correlation_matrix_FFphi_full[2][3]=0.091877;
	correlation_matrix_FFphi_full[2][4]=0.0390177;
	correlation_matrix_FFphi_full[2][5]=0.272613;
	correlation_matrix_FFphi_full[0][6]=0.996; /* was 1 */
	correlation_matrix_FFphi_full[0][7]=0.772589;
	correlation_matrix_FFphi_full[0][8]=0.42463;
	correlation_matrix_FFphi_full[1][6]=0.668299;
	correlation_matrix_FFphi_full[1][7]=0.791667;
	correlation_matrix_FFphi_full[1][8]=0.646183;
	correlation_matrix_FFphi_full[2][6]=-0.367399;
	correlation_matrix_FFphi_full[2][7]=-0.290877;
	correlation_matrix_FFphi_full[2][8]=0.166971;
	correlation_matrix_FFphi_full[0][9]=0.187815;
	correlation_matrix_FFphi_full[0][10]=-0.422447;
	correlation_matrix_FFphi_full[0][11]=-0.456553;
	correlation_matrix_FFphi_full[1][9]=0.0174825;
	correlation_matrix_FFphi_full[1][10]=-0.52219;
	correlation_matrix_FFphi_full[1][11]=-0.662488;
	correlation_matrix_FFphi_full[2][9]=0.207248;
	correlation_matrix_FFphi_full[2][10]=-0.0243012;
	correlation_matrix_FFphi_full[2][11]=-0.0754036;
	correlation_matrix_FFphi_full[0][12]=0.0865627;
	correlation_matrix_FFphi_full[0][13]=-0.297008;
	correlation_matrix_FFphi_full[0][14]=0.0671143;
	correlation_matrix_FFphi_full[1][12]=0.0545029;
	correlation_matrix_FFphi_full[1][13]=-0.247063;
	correlation_matrix_FFphi_full[1][14]=-0.0128991;
	correlation_matrix_FFphi_full[2][12]=0.182472;
	correlation_matrix_FFphi_full[2][13]=-0.00649649;
	correlation_matrix_FFphi_full[2][14]=-0.220386;
	correlation_matrix_FFphi_full[0][15]=0.0885951;
	correlation_matrix_FFphi_full[0][16]=-0.457242;
	correlation_matrix_FFphi_full[0][17]=-0.196984;
	correlation_matrix_FFphi_full[1][15]=0.0486156;
	correlation_matrix_FFphi_full[1][16]=-0.525849;
	correlation_matrix_FFphi_full[1][17]=-0.32988;
	correlation_matrix_FFphi_full[2][15]=0.174415;
	correlation_matrix_FFphi_full[2][16]=-0.0908376;
	correlation_matrix_FFphi_full[2][17]=-0.260704;
	correlation_matrix_FFphi_full[0][18]=0.703702;
	correlation_matrix_FFphi_full[0][19]=0.728003;
	correlation_matrix_FFphi_full[0][20]=0.602114;
	correlation_matrix_FFphi_full[1][18]=0.640094;
	correlation_matrix_FFphi_full[1][19]=0.775898;
	correlation_matrix_FFphi_full[1][20]=0.645661;
	correlation_matrix_FFphi_full[2][18]=-0.255547;
	correlation_matrix_FFphi_full[2][19]=-0.197701;
	correlation_matrix_FFphi_full[2][20]=-0.0842799;
	correlation_matrix_FFphi_full[3][3]=1.0;
	correlation_matrix_FFphi_full[3][4]=0.250919;
	correlation_matrix_FFphi_full[3][5]=-0.105449;
	correlation_matrix_FFphi_full[4][4]=1.0;
	correlation_matrix_FFphi_full[4][5]=0.78406;
	correlation_matrix_FFphi_full[5][5]=1.0;
	correlation_matrix_FFphi_full[3][6]=0.200928;
	correlation_matrix_FFphi_full[3][7]=-0.13586;
	correlation_matrix_FFphi_full[3][8]=-0.146239;
	correlation_matrix_FFphi_full[4][6]=-0.549042;
	correlation_matrix_FFphi_full[4][7]=-0.662113;
	correlation_matrix_FFphi_full[4][8]=-0.589264;
	correlation_matrix_FFphi_full[5][6]=-0.58909;
	correlation_matrix_FFphi_full[5][7]=-0.716176;
	correlation_matrix_FFphi_full[5][8]=-0.453332;
	correlation_matrix_FFphi_full[3][9]=0.406354;
	correlation_matrix_FFphi_full[3][10]=0.180876;
	correlation_matrix_FFphi_full[3][11]=0.109717;
	correlation_matrix_FFphi_full[4][9]=-0.00746984;
	correlation_matrix_FFphi_full[4][10]=0.642019;
	correlation_matrix_FFphi_full[4][11]=0.552109;
	correlation_matrix_FFphi_full[5][9]=0.0402144;
	correlation_matrix_FFphi_full[5][10]=0.47354;
	correlation_matrix_FFphi_full[5][11]=0.497547;
	correlation_matrix_FFphi_full[3][12]=0.267591;
	correlation_matrix_FFphi_full[3][13]=0.101937;
	correlation_matrix_FFphi_full[3][14]=-0.164558;
	correlation_matrix_FFphi_full[4][12]=0.142301;
	correlation_matrix_FFphi_full[4][13]=0.571088;
	correlation_matrix_FFphi_full[4][14]=-0.22105;
	correlation_matrix_FFphi_full[5][12]=0.156373;
	correlation_matrix_FFphi_full[5][13]=0.388671;
	correlation_matrix_FFphi_full[5][14]=-0.169048;
	correlation_matrix_FFphi_full[3][15]=0.275557;
	correlation_matrix_FFphi_full[3][16]=0.190196;
	correlation_matrix_FFphi_full[3][17]=-0.0562319;
	correlation_matrix_FFphi_full[4][15]=0.139241;
	correlation_matrix_FFphi_full[4][16]=0.794425;
	correlation_matrix_FFphi_full[4][17]=0.195856;
	correlation_matrix_FFphi_full[5][15]=0.148962;
	correlation_matrix_FFphi_full[5][16]=0.548463;
	correlation_matrix_FFphi_full[5][17]=0.124908;
	correlation_matrix_FFphi_full[3][18]=0.0931432;
	correlation_matrix_FFphi_full[3][19]=-0.0724199;
	correlation_matrix_FFphi_full[3][20]=-0.155631;
	correlation_matrix_FFphi_full[4][18]=-0.601549;
	correlation_matrix_FFphi_full[4][19]=-0.673956;
	correlation_matrix_FFphi_full[4][20]=-0.681441;
	correlation_matrix_FFphi_full[5][18]=-0.652782;
	correlation_matrix_FFphi_full[5][19]=-0.693798;
	correlation_matrix_FFphi_full[5][20]=-0.569813;
	correlation_matrix_FFphi_full[6][6]=1.0;
	correlation_matrix_FFphi_full[6][7]=0.775667;
	correlation_matrix_FFphi_full[6][8]=0.424398;
	correlation_matrix_FFphi_full[7][7]=1.0;
	correlation_matrix_FFphi_full[7][8]=0.803423;
	correlation_matrix_FFphi_full[8][8]=1.0;
	correlation_matrix_FFphi_full[6][9]=0.184513;
	correlation_matrix_FFphi_full[6][10]=-0.420811;
	correlation_matrix_FFphi_full[6][11]=-0.447059;
	correlation_matrix_FFphi_full[7][9]=-0.0661659;
	correlation_matrix_FFphi_full[7][10]=-0.502237;
	correlation_matrix_FFphi_full[7][11]=-0.553914;
	correlation_matrix_FFphi_full[8][9]=0.0196031;
	correlation_matrix_FFphi_full[8][10]=-0.450477;
	correlation_matrix_FFphi_full[8][11]=-0.475719;
	correlation_matrix_FFphi_full[6][12]=0.0773811;
	correlation_matrix_FFphi_full[6][13]=-0.304431;
	correlation_matrix_FFphi_full[6][14]=0.0722582;
	correlation_matrix_FFphi_full[7][12]=-0.148005;
	correlation_matrix_FFphi_full[7][13]=-0.406181;
	correlation_matrix_FFphi_full[7][14]= 0.0686011;
	correlation_matrix_FFphi_full[8][12]=-0.124325;
	correlation_matrix_FFphi_full[8][13]=-0.459704;
	correlation_matrix_FFphi_full[8][14]=-0.0415268;
	correlation_matrix_FFphi_full[6][15]=0.0795224;
	correlation_matrix_FFphi_full[6][16]=-0.455068;
	correlation_matrix_FFphi_full[6][17]=-0.187699;
	correlation_matrix_FFphi_full[7][15]=-0.149814;
	correlation_matrix_FFphi_full[7][16]=-0.597757;
	correlation_matrix_FFphi_full[7][17]=-0.236414;
	correlation_matrix_FFphi_full[8][15]=-0.129832;
	correlation_matrix_FFphi_full[8][16]=-0.651114;
	correlation_matrix_FFphi_full[8][17]=-0.338525;
	correlation_matrix_FFphi_full[6][18]=0.698961;
	correlation_matrix_FFphi_full[6][19]=0.723823;
	correlation_matrix_FFphi_full[6][20]=0.599831;
	correlation_matrix_FFphi_full[7][18]=0.629096;
	correlation_matrix_FFphi_full[7][19]=0.8456;
	correlation_matrix_FFphi_full[7][20]=0.733289;
	correlation_matrix_FFphi_full[8][18]=0.369588;
	correlation_matrix_FFphi_full[8][19]=0.656874;
	correlation_matrix_FFphi_full[8][20]=0.656771;
	correlation_matrix_FFphi_full[9][9]=1.0;
	correlation_matrix_FFphi_full[9][10]=0.247761;
	correlation_matrix_FFphi_full[9][11]=-0.123124;
	correlation_matrix_FFphi_full[10][10]=1.0;
	correlation_matrix_FFphi_full[10][11]=0.739804;
	correlation_matrix_FFphi_full[11][11]=1.0;
	correlation_matrix_FFphi_full[9][12]=0.526268;
	correlation_matrix_FFphi_full[9][13]=-0.0393931;
	correlation_matrix_FFphi_full[9][14]=-0.415068;
	correlation_matrix_FFphi_full[10][12]=0.218671;
	correlation_matrix_FFphi_full[10][13]=0.571275;
	correlation_matrix_FFphi_full[10][14]=-0.444945;
	correlation_matrix_FFphi_full[11][12]=-0.099119;
	correlation_matrix_FFphi_full[11][13]=0.323867;
	correlation_matrix_FFphi_full[11][14]=-0.0726924;
	correlation_matrix_FFphi_full[9][15]=0.535066;
	correlation_matrix_FFphi_full[9][16]=-0.0392003;
	correlation_matrix_FFphi_full[9][17]=-0.368212;
	correlation_matrix_FFphi_full[10][15]=0.223869;
	correlation_matrix_FFphi_full[10][16]=0.633854;
	correlation_matrix_FFphi_full[10][17]=-0.0725463;
	correlation_matrix_FFphi_full[11][15]=-0.0931479;
	correlation_matrix_FFphi_full[11][16]=0.50992;
	correlation_matrix_FFphi_full[11][17]=0.213584;
	correlation_matrix_FFphi_full[9][18]=0.160071;
	correlation_matrix_FFphi_full[9][19]=0.0458389;
	correlation_matrix_FFphi_full[9][20]=0.0186265;
	correlation_matrix_FFphi_full[10][18]=-0.617296;
	correlation_matrix_FFphi_full[10][19]=-0.617176;
	correlation_matrix_FFphi_full[10][20]=-0.629361;
	correlation_matrix_FFphi_full[11][18]=-0.677613;
	correlation_matrix_FFphi_full[11][19]=-0.651525;
	correlation_matrix_FFphi_full[11][20]=-0.547041;
	correlation_matrix_FFphi_full[12][12]=1.0;
	correlation_matrix_FFphi_full[12][13]=0.448903;
	correlation_matrix_FFphi_full[12][14]=-0.56142;
	correlation_matrix_FFphi_full[13][13]=1.0;
	correlation_matrix_FFphi_full[13][14]=-0.169757;
	correlation_matrix_FFphi_full[14][14]=1.0;
	correlation_matrix_FFphi_full[12][15]=0.996; /* was 1 */
	correlation_matrix_FFphi_full[12][16]=0.189067;
	correlation_matrix_FFphi_full[12][17]=-0.51613;
	correlation_matrix_FFphi_full[13][15]=0.423782;
	correlation_matrix_FFphi_full[13][16]=0.660917;
	correlation_matrix_FFphi_full[13][17]=-0.00861968;
	correlation_matrix_FFphi_full[14][15]=-0.561865;
	correlation_matrix_FFphi_full[14][16]=-0.137407;
	correlation_matrix_FFphi_full[14][17]=0.711765;
	correlation_matrix_FFphi_full[12][18]=0.106616;
	correlation_matrix_FFphi_full[12][19]=-0.176712;
	correlation_matrix_FFphi_full[12][20]=-0.210211;
	correlation_matrix_FFphi_full[13][18]=-0.321698;
	correlation_matrix_FFphi_full[13][19]=-0.487256;
	correlation_matrix_FFphi_full[13][20]=-0.54514;
	correlation_matrix_FFphi_full[14][18]=0.220799;
	correlation_matrix_FFphi_full[14][19]=0.232593;
	correlation_matrix_FFphi_full[14][20]=0.364997;
	correlation_matrix_FFphi_full[15][15]=1.0;
	correlation_matrix_FFphi_full[15][16]=0.212618;
	correlation_matrix_FFphi_full[15][17]=-0.495299;
	correlation_matrix_FFphi_full[16][16]=1.0;
	correlation_matrix_FFphi_full[16][17]=0.410351;
	correlation_matrix_FFphi_full[17][17]=1.0;
	correlation_matrix_FFphi_full[15][18]=0.109431;
	correlation_matrix_FFphi_full[15][19]=-0.176517;
	correlation_matrix_FFphi_full[15][20]=-0.208143;
	correlation_matrix_FFphi_full[16][18]=-0.455203;
	correlation_matrix_FFphi_full[16][19]=-0.62091;
	correlation_matrix_FFphi_full[16][20]=-0.661853;
	correlation_matrix_FFphi_full[17][18]=-0.0306923;
	correlation_matrix_FFphi_full[17][19]=-0.114741;
	correlation_matrix_FFphi_full[17][20]=0.00905559;
	correlation_matrix_FFphi_full[18][18]=1.0;
	correlation_matrix_FFphi_full[18][19]=0.777787;
	correlation_matrix_FFphi_full[18][20]=0.573289;
	correlation_matrix_FFphi_full[19][19]=1.0;
	correlation_matrix_FFphi_full[19][20]=0.899178;
	correlation_matrix_FFphi_full[20][20]=1.0;
    for(ie=0;ie<21;ie++) for(je=ie+1;je<21;je++) correlation_matrix_FFphi_full[je][ie]=correlation_matrix_FFphi_full[ie][je];    

	double correlation_matrix_FFphi_soft[9][9];
	for(ie=0;ie<9;ie++) for(je=0;je<9;je++) correlation_matrix_FFphi_soft[ie][je]=(ie==je);
	correlation_matrix_FFphi_soft[3-3][3-3]=correlation_matrix_FFphi_full[3][3];
	correlation_matrix_FFphi_soft[3-3][4-3]=correlation_matrix_FFphi_full[3][4];
	correlation_matrix_FFphi_soft[3-3][5-3]=correlation_matrix_FFphi_full[3][5];
	correlation_matrix_FFphi_soft[4-3][4-3]=correlation_matrix_FFphi_full[4][4];
	correlation_matrix_FFphi_soft[4-3][5-3]=correlation_matrix_FFphi_full[4][5];
	correlation_matrix_FFphi_soft[5-3][5-3]=correlation_matrix_FFphi_full[5][5];
	correlation_matrix_FFphi_soft[3-3][6-3]=correlation_matrix_FFphi_full[3][6];
	correlation_matrix_FFphi_soft[3-3][7-3]=correlation_matrix_FFphi_full[3][7];
	correlation_matrix_FFphi_soft[3-3][8-3]=correlation_matrix_FFphi_full[3][8];
	correlation_matrix_FFphi_soft[4-3][6-3]=correlation_matrix_FFphi_full[4][6];
	correlation_matrix_FFphi_soft[4-3][7-3]=correlation_matrix_FFphi_full[4][7];
	correlation_matrix_FFphi_soft[4-3][8-3]=correlation_matrix_FFphi_full[4][8];
	correlation_matrix_FFphi_soft[5-3][6-3]=correlation_matrix_FFphi_full[5][6];
	correlation_matrix_FFphi_soft[5-3][7-3]=correlation_matrix_FFphi_full[5][7];
	correlation_matrix_FFphi_soft[5-3][8-3]=correlation_matrix_FFphi_full[5][8];
	correlation_matrix_FFphi_soft[3-3][9-3]=correlation_matrix_FFphi_full[3][9];
	correlation_matrix_FFphi_soft[3-3][10-3]=correlation_matrix_FFphi_full[3][10];
	correlation_matrix_FFphi_soft[3-3][11-3]=correlation_matrix_FFphi_full[3][11];
	correlation_matrix_FFphi_soft[4-3][9-3]=correlation_matrix_FFphi_full[4][9];
	correlation_matrix_FFphi_soft[4-3][10-3]=correlation_matrix_FFphi_full[4][10];
	correlation_matrix_FFphi_soft[4-3][11-3]=correlation_matrix_FFphi_full[4][11];
	correlation_matrix_FFphi_soft[5-3][9-3]=correlation_matrix_FFphi_full[5][9];
	correlation_matrix_FFphi_soft[5-3][10-3]=correlation_matrix_FFphi_full[5][10];
	correlation_matrix_FFphi_soft[5-3][11-3]=correlation_matrix_FFphi_full[5][11];
	correlation_matrix_FFphi_soft[6-3][6-3]=correlation_matrix_FFphi_full[6][6];
	correlation_matrix_FFphi_soft[6-3][7-3]=correlation_matrix_FFphi_full[6][7];
	correlation_matrix_FFphi_soft[6-3][8-3]=correlation_matrix_FFphi_full[6][8];
	correlation_matrix_FFphi_soft[7-3][7-3]=correlation_matrix_FFphi_full[7][7];
	correlation_matrix_FFphi_soft[7-3][8-3]=correlation_matrix_FFphi_full[7][8];
	correlation_matrix_FFphi_soft[8-3][8-3]=correlation_matrix_FFphi_full[8][8];
	correlation_matrix_FFphi_soft[6-3][9-3]=correlation_matrix_FFphi_full[6][9];
	correlation_matrix_FFphi_soft[6-3][10-3]=correlation_matrix_FFphi_full[6][10];
	correlation_matrix_FFphi_soft[6-3][11-3]=correlation_matrix_FFphi_full[6][11];
	correlation_matrix_FFphi_soft[7-3][9-3]=correlation_matrix_FFphi_full[7][9];
	correlation_matrix_FFphi_soft[7-3][10-3]=correlation_matrix_FFphi_full[7][10];
	correlation_matrix_FFphi_soft[7-3][11-3]=correlation_matrix_FFphi_full[7][11];
	correlation_matrix_FFphi_soft[8-3][9-3]=correlation_matrix_FFphi_full[8][9];
	correlation_matrix_FFphi_soft[8-3][10-3]=correlation_matrix_FFphi_full[8][10];
	correlation_matrix_FFphi_soft[8-3][11-3]=correlation_matrix_FFphi_full[8][11];
	correlation_matrix_FFphi_soft[9-3][9-3]=correlation_matrix_FFphi_full[9][9];
	correlation_matrix_FFphi_soft[9-3][10-3]=correlation_matrix_FFphi_full[9][10];
	correlation_matrix_FFphi_soft[9-3][11-3]=correlation_matrix_FFphi_full[9][11];
	correlation_matrix_FFphi_soft[10-3][10-3]=correlation_matrix_FFphi_full[10][10];
	correlation_matrix_FFphi_soft[10-3][11-3]=correlation_matrix_FFphi_full[10][11];
	correlation_matrix_FFphi_soft[11-3][11-3]=correlation_matrix_FFphi_full[11][11];
    for(ie=0;ie<9;ie++) for(je=ie+1;je<9;je++) correlation_matrix_FFphi_soft[je][ie]=correlation_matrix_FFphi_soft[ie][je];
	
    /* correlation matrix for B-> K ll from 1411.3161 */
	double correlation_matrix_FFK_tmp[10][10]=
	{{1.0,-0.39,-0.71,-0.63,0.49,-0.03,-0.22,0.16,-0.08,-0.09},
	{-0.39,1.0,0.66,0.26,0.05,0.72,0.48,-0.08,0.03,0.01},
	{-0.71,0.66,1.0,0.54,-0.17,0.51,0.59,-0.16,0.05,0.09},
	{-0.63,0.26,0.54,1.0,0.05,0.14,0.05,0.0,0.03,-0.01},
	{0.49,0.05,-0.17,0.05,1.0,0.09,-0.47,0.34,-0.06,-0.28},
	{-0.03,0.72,0.51,0.14,0.09,1.0,0.43,-0.06,0.11,-0.04},
	{-0.22,0.48,0.59,0.05,-0.47,0.43,1.0,-0.32,-0.05,0.29},
	{0.16,-0.08,-0.16,0.0,0.34,-0.06,-0.32,1.0,0.0,-0.35},
	{-0.08,0.03,0.05,0.03,-0.06,0.11,-0.05,0.0,1.0,0.21},
	{-0.09,0.01,0.09,-0.01,-0.28,-0.04,0.29,-0.35,0.21,1.0}};
	for(ie=0;ie<10;ie++) for(je=0;je<10;je++) correlation_matrix_FFK[ie][je]=correlation_matrix_FFK_tmp[ie][je];

	if(fullFF)
	{
			for(ie=0;ie<21;ie++) for(je=0;je<21;je++)
			{
				correlation_matrix_FFKstar[ie][je]=correlation_matrix_FFKstar_full[ie][je];
				correlation_matrix_FFphi[ie][je]=correlation_matrix_FFphi_full[ie][je];
			}
	}
	else
	{
			for(ie=0;ie<9;ie++) for(je=0;je<9;je++)
			{
				correlation_matrix_FFKstar[ie][je]=correlation_matrix_FFKstar_soft[ie][je];
				correlation_matrix_FFphi[ie][je]=correlation_matrix_FFphi_soft[ie][je];
			}
	}

	return;
}

/*---------------------------------------------------------------------*/

void experimental_covariance(int likelihood, char* namesin[], int nbobsin, double* central_exp, double* errors_exp, double **correlations)
{	
	int ie,je;
	
	int translate[nbobsin];
		
	for(ie=0;ie<nbobsin;ie++) central_exp[ie]=errors_exp[ie]=0.;
	
	double stat2BRBsll,stat2BRB0Kstar0,stat2AngB0Kstar0,stat2BRBpKstarp,stat2BRBKstaree,stat2AngBKstaree,stat2BRBK,stat2BRBpKp,statBRBpKpee,statRK,statBRBsPhi,statAngBsPhi,stat2RKstar;
    double syst2BRBsll,syst2BRB0Kstar0,syst2AngB0Kstar0,syst2BRBpKstarp,syst2BRBKstaree,syst2AngBKstaree,syst2BRBK,syst2BRBpKp,systBRBpKpee,systRK,systBRBsPhi,systAngBsPhi,syst2RKstar;
    double other2BRB0Kstar0,other2BRBKstaree,otherBRBsPhi;

	stat2BRBsll=stat2BRB0Kstar0=stat2AngB0Kstar0=stat2BRBpKstarp=stat2BRBKstaree=stat2AngBKstaree=stat2BRBK=stat2BRBpKp=statBRBpKpee=statRK=statBRBsPhi=statAngBsPhi=stat2RKstar=1.;
	syst2BRBsll=syst2BRB0Kstar0=syst2AngB0Kstar0=syst2BRBpKstarp=syst2BRBKstaree=syst2AngBKstaree=syst2BRBK=syst2BRBpKp=systBRBpKpee=systRK=systBRBsPhi=systAngBsPhi=syst2RKstar=1.;
	other2BRB0Kstar0=other2BRBKstaree=otherBRBsPhi=1.;
	
	char* names[NBOBSMAX];
	double central[NBOBSMAX],errors[NBOBSMAX];

#include "experimental_input.h"
	
	for(ie=0;ie<nbobsin;ie++) translate[ie]=-1;

	for(ie=0;ie<nbobsin;ie++)
	{
		je=0;
		while(je<NBOBSMAX)
		{
			if(names[je]!=NULL)
			{
				if(!strcmp(namesin[ie],names[je]))
				{
					translate[ie]=je;
					je=NBOBSMAX;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsin;ie++)
	{
		if(translate[ie]>-1&&errors[translate[ie]]!=0.)
		{
			central_exp[ie]=central[translate[ie]];
			errors_exp[ie]=errors[translate[ie]];
		}
		else 
		{
			printf("Warning! %s not present (or experimental error set to 0) in experimental_input.h. Error set to 10^30...\n",namesin[ie]);  
			central_exp[ie]=0.;
			errors_exp[ie]=1.e30;
		}
	}

	for(ie=0;ie<nbobsin;ie++) for(je=0;je<nbobsin;je++) correlations[ie][je]=(ie==je);

	double correlations_exp[NBOBSMAX][NBOBSMAX];
	for(ie=0;ie<NBOBSMAX;ie++) for(je=0;je<NBOBSMAX;je++) correlations_exp[ie][je]=(ie==je);	
#include "correlations_exp_input.h"
	
	for(ie=0;ie<nbobsin;ie++)
	{
		for(je=0;je<nbobsin;je++)
		{
			if((translate[ie]==-1||translate[je]==-1))
			{
				correlations[ie][je]=(ie==je);
			}
			else
			{
				correlations[ie][je]=correlations_exp[translate[ie]][translate[je]];
			}
		}
	}

	return;
}

/*---------------------------------------------------------------------*/

int get_covariance(double ***covariance_th, double ***covariance_exp, double **central_exp, char* namesin[], int nbobsin, double complex deltaC[], double complex deltaCp[], double complex deltaCQ[], double complex deltaCQp[], struct parameters* param, int *nbobs_out)
{		
	double values_ref[NBOBSMAX];
	int onoff[NBOBSMAX];
	char *names0[NBOBSMAX],*namesout[NBOBSMAX];
	
	int ie,je,ke,le;
	
	for(ie=0;ie<NBOBSMAX;ie++) onoff[ie]=1;
	
	int nbobs0=observables(10000,NULL,NULL,NULL,NULL,names0,values_ref,onoff,values_ref,param); /* get the name and number of observables with experimental values and correlations */

	activate_observables(nbobs0,names0,nbobsin,namesin,onoff);
		
	int n; /* number of iterations */
	int nsoft=132;
	int nfull=156; /* nsoft + 2*(21-9) */
	if(param->fullFF) n=nfull; else n=nsoft;
	
	int nbobs=observables(0,deltaC,deltaCp,deltaCQ,deltaCQp,namesout,values_ref,onoff,NULL,param);
	
	if(*nbobs_out!=nbobs) return -1;
	
		
	double **values_mod;
	values_mod=(double **) malloc((n+1)*sizeof(double *));
	for(ie=0;ie<=n;ie++) values_mod[ie]=(double *) malloc(nbobs*sizeof(double));	

	for(ie=0;ie<nbobs;ie++) values_mod[0][ie]=values_ref[ie];		

	int nFF; /* size of correlation matrix for the FF */
	if(param->fullFF) nFF=21; else nFF=9;
	double **correlation_matrix_FFKstar,**correlation_matrix_FFphi,**correlation_matrix_FFK;
	correlation_matrix_FFKstar=(double **) malloc(nFF*sizeof(double *));
	correlation_matrix_FFphi=(double **) malloc(nFF*sizeof(double *));
	correlation_matrix_FFK=(double **) malloc(10*sizeof(double *));
	for(ie=0;ie<nFF;ie++) correlation_matrix_FFKstar[ie]=(double *) malloc(nFF*sizeof(double));	
	for(ie=0;ie<nFF;ie++) correlation_matrix_FFphi[ie]=(double *) malloc(nFF*sizeof(double));	
	for(ie=0;ie<10;ie++) correlation_matrix_FFK[ie]=(double *) malloc(10*sizeof(double));	

	correlation_matrices(param->fullFF,correlation_matrix_FFKstar,correlation_matrix_FFphi,correlation_matrix_FFK);
	
	double **corr; /* correlation matrix between the parameters */
	corr=(double **) malloc((n+1)*sizeof(double *));
	for(ie=0;ie<=n;ie++) corr[ie]=(double *) malloc((n+1)*sizeof(double));	
	
	for(ie=1;ie<=n;ie++) for(je=1;je<=n;je++) corr[ie][je]=(ie==je);

	for(ie=0;ie<nFF;ie++) for(je=0;je<nFF;je++) corr[114+ie][114+je]=correlation_matrix_FFKstar[ie][je];
	for(ie=0;ie<nFF;ie++) for(je=0;je<nFF;je++) corr[114+nFF+ie][114+nFF+je]=correlation_matrix_FFphi[ie][je];

#if defined(_OPENMP)	
#pragma omp parallel for private(ie)
#endif
	for(ke=1;ke<=n;ke++)
	{			
		double values[NBOBSMAX];
					
		observables(ke,deltaC,deltaCp,deltaCQ,deltaCQp,namesout,values,onoff,values_ref,param);

		for(ie=0;ie<nbobs;ie++) values_mod[ke][ie]=values[ie];
	}	
		
	*covariance_th=NULL;
	*covariance_th=(double **) malloc(nbobs*sizeof(double *));
	for(ie=0;ie<nbobs;ie++) (*covariance_th)[ie]=(double *) malloc(nbobs*sizeof(double));		

	for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++)
	{
		(*covariance_th)[ie][je]=0.;
		for(ke=1;ke<=n;ke++) for(le=1;le<=n;le++) (*covariance_th)[ie][je]+=corr[ke][le]*(values_mod[ke][ie]-values_ref[ie])*(values_mod[le][je]-values_ref[je]);
	}

	*covariance_exp=NULL;
	*covariance_exp=(double **) malloc(nbobs*sizeof(double *));
	for(ie=0;ie<nbobs;ie++) (*covariance_exp)[ie]=(double *) malloc(nbobs*sizeof(double));		

	*central_exp=NULL;
	*central_exp=(double *) malloc(nbobs*sizeof(double));		

	double *sigma_exp=NULL;	
	sigma_exp=(double *) malloc(nbobs*sizeof(double));		
 
	experimental_covariance(param->likelihoodBKstarmumu,namesout,nbobs,*central_exp,sigma_exp,*covariance_exp);
				
	for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++) (*covariance_exp)[ie][je]*=sigma_exp[ie]*sigma_exp[je];
		
	*nbobs_out=nbobs;
	
	return 1;
}

/*---------------------------------------------------------------------*/

void get_covtot(double ***covariance_th, double ***covariance_exp, double ***covariance_tot, int nbobs)
{		
	int ie,je;
	
	*covariance_tot=NULL;
	*covariance_tot=(double **) malloc(nbobs*sizeof(double *));
	for(ie=0;ie<nbobs;ie++) (*covariance_tot)[ie]=(double *) malloc(nbobs*sizeof(double));		

	
	for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++) (*covariance_tot)[ie][je]=(*covariance_exp)[ie][je]+(*covariance_th)[ie][je];	
		
	return;
}

/*---------------------------------------------------------------------*/

int reduce_covariance(double ***covariance_in, char* namesin[], int nbobsin, double ***covariance_out, char* namesout[], int nbobsout)
{		
	if(nbobsout>nbobsin) return 0;

	int ie,je;
	
	int translate[nbobsout];
	for(ie=0;ie<nbobsout;ie++) translate[ie]=-1;
	
	for(ie=0;ie<nbobsout;ie++)
	{
		je=0;
		while(je<nbobsin)
		{
			if(namesin[je]!=NULL)
			{
				if(!strcmp(namesout[ie],namesin[je]))
				{
					translate[ie]=je;
					je=nbobsin;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsout;ie++) if(translate[ie]==-1) return 0;
	
	*covariance_out=NULL;
	*covariance_out=(double **) malloc(nbobsout*sizeof(double *));
	for(ie=0;ie<nbobsout;ie++) (*covariance_out)[ie]=(double *) malloc(nbobsout*sizeof(double));		

	for(ie=0;ie<nbobsout;ie++)  for(je=0;je<nbobsout;je++) (*covariance_out)[ie][je]=(*covariance_in)[translate[ie]][translate[je]];
	return 1;
}

/*---------------------------------------------------------------------*/

int reduce_values(double **values_in, char* namesin[], int nbobsin, double **values_out, char* namesout[], int nbobsout)
{		
	if(nbobsout>nbobsin) return 0;

	int ie,je;
	
	int translate[nbobsout];
	for(ie=0;ie<nbobsout;ie++) translate[ie]=-1;
	
	for(ie=0;ie<nbobsout;ie++)
	{
		je=0;
		while(je<nbobsin)
		{
			if(namesin[je]!=NULL)
			{
				if(!strcmp(namesout[ie],namesin[je]))
				{
					translate[ie]=je;
					je=nbobsin;
				}
				else je++;
			} 
			else je++;
		}
	}

	for(ie=0;ie<nbobsout;ie++) if(translate[ie]==-1) return 0;

	*values_out=NULL;
	*values_out=(double *) malloc(nbobsout*sizeof(double));		

	for(ie=0;ie<nbobsout;ie++)  (*values_out)[ie]=(*values_in)[translate[ie]];
	
	return 1;
}

/*---------------------------------------------------------------------*/

int get_invcovtot(double ***covariance_tot, double ***inv_cov_tot, int nbobs)
{	
	int ie;
				
	*inv_cov_tot=NULL;
	*inv_cov_tot=(double **) malloc(nbobs*sizeof(double *));
	for(ie=0;ie<nbobs;ie++) (*inv_cov_tot)[ie]=(double *) malloc(nbobs*sizeof(double));		

	if(invert_matrix(nbobs,*covariance_tot,*inv_cov_tot)==0) return -1; /* inversion of the total covariance matrix */
	else return 1;
}

/*---------------------------------------------------------------------*/

void print_covariance(char filename[], double **covariance_in, double *central_in, char* namesin[], int nbobsin)
{		
	int ie,je;
	
	FILE *output;
	
	output=fopen(filename,"w");
	
	fprintf(output,"NOBS %d\n\n",nbobsin);
	
	for(ie=0;ie<nbobsin;ie++) fprintf(output,"%d\t%-25s\t%.4e\n",ie,namesin[ie],central_in[ie]);
	fprintf(output,"\n");
	
	for(ie=0;ie<nbobsin;ie++) for(je=0;je<=ie;je++) fprintf(output,"%d\t%d\t%.4e\n",ie,je,covariance_in[ie][je]);
	
	fclose(output);
	
	return;
}

/*---------------------------------------------------------------------*/

int read_nbobs(char filename[])
{		
	if(!test_file(filename)) return 0;
	
	int nbobs;
	char dummy[50];
	
	FILE *input;
	input=fopen(filename,"r");

	fscanf(input,"%s",dummy);	
	fscanf(input,"%d",&nbobs);
	
	fclose(input);
	
	return nbobs;
}

/*---------------------------------------------------------------------*/

int read_covariance(char filename[], double ***covariance, double **central, char* names[], int *nbobs)
{		
	if(!test_file(filename)) return 0;
	
	int ie,je,nbobstmp;
	char dummy[50];
	
	FILE *input;
	input=fopen(filename,"r");

	fscanf(input,"%s",dummy);	
	fscanf(input,"%d",&nbobstmp);
	
	*nbobs=nbobstmp;

	*central=NULL;
	*central=(double *) malloc(nbobstmp*sizeof(double));		
	
	for(ie=0;ie<nbobstmp;ie++) 
	{
		fscanf(input,"%s",dummy);
		
		fscanf(input,"%s",dummy);
		names[ie]=malloc(strlen(dummy)+1);
		strcpy(names[ie], dummy);
		
		fscanf(input,"%lf",&(*central)[ie]);
	}

	*covariance=NULL;
	*covariance=(double **) malloc(nbobstmp*sizeof(double *));
	for(ie=0;ie<nbobstmp;ie++) (*covariance)[ie]=(double *) malloc(nbobstmp*sizeof(double));		
		
	for(ie=0;ie<nbobstmp;ie++) for(je=0;je<=ie;je++)
	{
		fscanf(input,"%s",dummy);
		fscanf(input,"%s",dummy);
		fscanf(input,"%lf",&(*covariance)[ie][je]);
		if(ie!=je) (*covariance)[je][ie]=(*covariance)[ie][je];
	}	
	fclose(input);
	
	return 1;
}

/*---------------------------------------------------------------------*/

int get_predictions(char* namesin[], int nbobsin, double** predictions, double complex deltaC[], double complex deltaCp[], double complex deltaCQ[], double complex deltaCQp[], struct parameters* param)
{		
 	double values_ref[NBOBSMAX];
	int onoff[NBOBSMAX];
	char *names0[NBOBSMAX];
	
	int ie;
	
	for(ie=0;ie<NBOBSMAX;ie++) onoff[ie]=1;
	
	int nbobs0=observables(10000,NULL,NULL,NULL,NULL,names0,values_ref,onoff,values_ref,param); /* get the name and number of observables with experimental values and correlations */

	activate_observables(nbobs0,names0,nbobsin,namesin,onoff);
		
	int nbobs=observables(0,deltaC,deltaCp,deltaCQ,deltaCQp,namesin,values_ref,onoff,NULL,param); /* give number of observables */
		
	*predictions=NULL;
	*predictions=(double *) malloc(nbobs*sizeof(double));	
	
	for(ie=0;ie<nbobs;ie++) (*predictions)[ie]=values_ref[ie];

	return nbobs;
}
/*---------------------------------------------------------------------*/

double get_chi2(double **inv_cov_tot, double *predictions, double *central_exp, int nbobs)
{		
	int ie,je;
	double chi2=0.;
	
	for(ie=0;ie<nbobs;ie++) for(je=0;je<nbobs;je++) chi2+=(predictions[ie]-central_exp[ie])*inv_cov_tot[ie][je]*(predictions[je]-central_exp[je]);
		
	return chi2;
}

/*---------------------------------------------------------------------*/

double chi2(char* namesin[], int nbobsin, struct parameters* param, int *nbobs)
{
	double complex deltaC[21],deltaCp[21],deltaCQ[5],deltaCQp[5];
	int ie;
	double *predictions;
	
	for(ie=1;ie<=20;ie++)
	{
		deltaC[ie]=param->deltaC[ie];
		deltaCp[ie]=param->deltaCp[ie];
	}
	for(ie=1;ie<=4;ie++)
	{
		deltaCQ[ie]=param->deltaCQ[ie];
		deltaCQp[ie]=param->deltaCQp[ie];
	}
	
	*nbobs=get_predictions(namesin,nbobsin,&predictions,deltaC,deltaCp,deltaCQ,deltaCQp,param);

	double **covariance_th,**covariance_exp,*central_exp;
	if(get_covariance(&covariance_th,&covariance_exp,&central_exp,namesin,nbobsin,deltaC,deltaCp,deltaCQ,deltaCQp,param,nbobs)<0) return -2.;

#ifdef DEBUG
	for(ie=0;ie<*nbobs;ie++) printf("\t%-25s:\tprediction %.3e\tcentral %.3e\texperr %.3e\ttherr %.3e\n",namesin[ie],predictions[ie],central_exp[ie],sqrt(covariance_exp[ie][ie]),sqrt(covariance_th[ie][ie]));
#endif

	double **covariance_tot;
	get_covtot(&covariance_th,&covariance_exp,&covariance_tot,*nbobs);

	double **inv_cov_tot;	
	if(get_invcovtot(&covariance_tot,&inv_cov_tot,*nbobs)<0) return -3.;
			
	return get_chi2(inv_cov_tot,predictions,central_exp,nbobsin);	
}
