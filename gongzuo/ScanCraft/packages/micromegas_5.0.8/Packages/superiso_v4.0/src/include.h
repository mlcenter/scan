#ifndef H_GLOBAL
#define H_GLOBAL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <complex.h>
#include <string.h>
#include <strings.h>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/
/*#define SM_ChargedHiggs*/

/*--------------------------------------------------------------------*/

#define pi    3.1415926535897932385
#define zeta3 1.2020569031595942855
#define hbar  6.58211889e-25 /* in GeV.s */
#define Gn    6.67428e-8  /* in cm^3.g^-1.s^-2 */
#define Mplanck     1.2209102930946623e+19 /* in GeV, more precise definition than before */

#define MN 0.939 // nucleon mass in GeV
#define mp 0.9382720
#define mn 0.9395654

#define hbarc 197327 //keV fm


#define Nobs_BKll 2
#define Nobs_BKsll 30
#define Nobs_Bsphill 6
#define Nobs_BDlnu 2
#define Nobs_BDstarlnu 3

#define NBOBSMAX 500

/*--------------------------------------------------------------------*/

typedef struct parameters
/* structure containing all the scanned parameters from the SLHA file */
{
	int SM;
	int model; /* CMSSM=1, GMSB=2, AMSB=3 */
	int generator; /* ISAJET=1, SOFTSUSY=3, SPHENO=4, SUSPECT=5, NMSSMTOOLS=6 */
	double Q; /* Qmax ; default = M_EWSB = sqrt(m_stop1*mstop2) */

	double m0,m12,tan_beta,sign_mu,A0; /* CMSSM parameters */
	double Lambda,Mmess,N5,cgrav,m32; /* AMSB, GMSB parameters */
	double mass_Z,mass_W,mass_b,mass_top_pole,mass_tau_pole; /* SM parameters */
	double inv_alpha_em,alphas_MZ,Gfermi,GAUGE_Q; /* SM parameters */
	double charg_Umix[3][3],charg_Vmix[3][3],stop_mix[3][3],sbot_mix[3][3],stau_mix[3][3],neut_mix[6][6],mass_neut[6],alpha; /* mass mixing matrices */
	double Min,M1_Min,M2_Min,M3_Min,At_Min,Ab_Min,Atau_Min,M2H1_Min,M2H2_Min,mu_Min,M2A_Min,tb_Min,mA_Min; /* optional input parameters at scale Min */
	double MeL_Min,MmuL_Min,MtauL_Min,MeR_Min,MmuR_Min,MtauR_Min; /* optional input parameters at scale Min */
	double MqL1_Min,MqL2_Min,MqL3_Min,MuR_Min,McR_Min,MtR_Min,MdR_Min,MsR_Min,MbR_Min; /* optional input parameters at scale Min */
	double N51,N52,N53,M2H1_Q,M2H2_Q; /* optional input parameters (N51...3: GMSB) */
	double mass_d,mass_u,mass_s,mass_c,mass_t,mass_e,mass_nue,mass_mu,mass_num,mass_tau,mass_nut; /* SM masses */
	double mass_gluon,mass_photon,mass_Z0; /* SM masses */
	double mass_h0,mass_H0,mass_A0,mass_H,mass_dnl,mass_upl,mass_stl,mass_chl,mass_b1,mass_t1; /* Higgs & superparticle masses */
	double mass_el,mass_nuel,mass_mul,mass_numl,mass_tau1,mass_nutl,mass_gluino,mass_cha1,mass_cha2; /* superparticle masses */
	double mass_dnr,mass_upr,mass_str,mass_chr,mass_b2,mass_t2,mass_er,mass_mur,mass_tau2; /* superparticle masses */
	double mass_nuer,mass_numr,mass_nutr,mass_graviton,mass_gravitino; /* superparticle masses */
	double gp,g2,gp_Q,g2_Q,g3_Q,YU_Q,yut[4],YD_Q,yub[4],YE_Q,yutau[4]; /* couplings */
	double HMIX_Q,mu_Q,tanb_GUT,Higgs_VEV,mA2_Q,MSOFT_Q,M1_Q,M2_Q,M3_Q; /* parameters at scale Q */
	double MeL_Q,MmuL_Q,MtauL_Q,MeR_Q,MmuR_Q,MtauR_Q,MqL1_Q,MqL2_Q,MqL3_Q,MuR_Q,McR_Q,MtR_Q,MdR_Q,MsR_Q,MbR_Q; /* masses at scale Q */
	double AU_Q,A_u,A_c,A_t,AD_Q,A_d,A_s,A_b,AE_Q,A_e,A_mu,A_tau; /* trilinear couplings */
	
	/* SLHA2 */
	int NMSSM,RV,CPV,FV;
	double CKM_lambda,CKM_A,CKM_rhobar,CKM_etabar;
	double PMNS_theta12,PMNS_theta23,PMNS_theta13,PMNS_delta13,PMNS_alpha1,PMNS_alpha2;
	double lambdaNMSSM_Min,kappaNMSSM_Min,AlambdaNMSSM_Min,AkappaNMSSM_Min,lambdaSNMSSM_Min,xiFNMSSM_Min,xiSNMSSM_Min,mupNMSSM_Min,mSp2NMSSM_Min,mS2NMSSM_Min,mass_H03,mass_A02,NMSSMRUN_Q,lambdaNMSSM,kappaNMSSM,AlambdaNMSSM,AkappaNMSSM,lambdaSNMSSM,xiFNMSSM,xiSNMSSM,mupNMSSM,mSp2NMSSM,mS2NMSSM; /* NMSSM parameters */
	double PMNSU_Q,CKM_Q,IMCKM_Q,MSE2_Q,MSU2_Q,MSD2_Q,MSL2_Q,MSQ2_Q,TU_Q,TD_Q,TE_Q;
	double CKM[4][4],IMCKM[4][4]; /* CKM matrix */
	double H0_mix[4][4],A0_mix[4][4]; /* Higgs mixing matrices */
	double sU_mix[7][7],sD_mix[7][7],sE_mix[7][7], sNU_mix[4][4]; /* mixing matrices */
	double sCKM_msq2[4][4],sCKM_msl2[4][4],sCKM_msd2[4][4],sCKM_msu2[4][4],sCKM_mse2[4][4]; /* super CKM matrices */
	double PMNS_U[4][4]; /* PMNS mixing matrices */
	double TU[4][4],TD[4][4],TE[4][4]; /* trilinear couplings */
	
	/* non-SLHA*/
	double mass_c_pole,mass_b_1S,mass_b_pole,mtmt;
	int scheme_c_mass;
	double Lambda3,Lambda4,Lambda5,Lambda6; /* Lambda QCD */
	double alphasMZ_Lambda3,alphasMZ_Lambda4,alphasMZ_Lambda5,alphasMZ_Lambda6; /* Lambda QCD */
	
	/* Flavour constants */
	double f_B,f_Bs,f_Ds,f_D,fK_fpi,f_K;
	double f_Kstar_par,f_Kstar_perp;
	double f_phi_perp,f_phi_par,a1phi_perp,a2phi_perp,a1phi_par,a2phi_par;
	double m_B,m_Bs,m_Bd,m_pi,m_Ds,m_K0,m_K,m_Kstar0,m_Kstar,m_D0,m_D,m_Dstar0,m_Dstar,m_phi;
	double life_pi,life_K,life_B,life_Bs,life_Bd,life_D,life_Ds;
	double a1par,a2par,a1perp,a2perp,a1K,a2K;
	double zeta3A,zeta3V,wA10,deltatp,deltatm,deltatp_phi,deltatm_phi;
	double lambda_Bp,lambda_Bsp,rho1,lambda2;
	double BR_BXclnu_exp; /* Used in bsgamma.c and bsll.c */
	int fullFF; /* full or soft form factor approach */
	
	/* CKM matrix */
	double complex Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb;

	/* b -> s gamma */
	double mu_G2_bsg,rho_D3_bsg,rho_LS3_bsg,mu_c_bsg;
	
	/* B->K* gamma */
	double T1_BKstar;
	
	/* B -> Kstar FF */
	double a0V_BKstar,a1V_BKstar,a2V_BKstar,MV_BKstar;
	double a0A1_BKstar,a1A1_BKstar,a2A1_BKstar,MA1_BKstar;	
	double a0A12_BKstar,a1A12_BKstar,a2A12_BKstar,MA12_BKstar;
	double a0A0_BKstar,a1A0_BKstar,a2A0_BKstar,MA0_BKstar;
	double a0T1_BKstar,a1T1_BKstar,a2T1_BKstar,MT1_BKstar;
	double a0T2_BKstar,a1T2_BKstar,a2T2_BKstar,MT2_BKstar;
	double a0T23_BKstar,a1T23_BKstar,a2T23_BKstar,MT23_BKstar;

	/* Bs -> phi FF */
	double a0V_Bsphi,a1V_Bsphi,a2V_Bsphi,MV_Bsphi;
	double a0A1_Bsphi,a1A1_Bsphi,a2A1_Bsphi,MA1_Bsphi;	
	double a0A12_Bsphi,a1A12_Bsphi,a2A12_Bsphi,MA12_Bsphi;
	double a0A0_Bsphi,a1A0_Bsphi,a2A0_Bsphi,MA0_Bsphi;
	double a0T1_Bsphi,a1T1_Bsphi,a2T1_Bsphi,MT1_Bsphi;
	double a0T2_Bsphi,a1T2_Bsphi,a2T2_Bsphi,MT2_Bsphi;
	double a0T23_Bsphi,a1T23_Bsphi,a2T23_Bsphi,MT23_Bsphi;
	
	/* B -> K FF */
	double a00_BK,a10_BK,a20_BK,a30_BK;
	double a0p_BK,a1p_BK,a2p_BK,DmBp_BK;
	double a0T_BK,a1T_BK,a2T_BK,DmBT_BK;
	
	/* B -> Kstar hadronic uncertainties */
	double hadrerrBKstar;
	double complex BtoKstarlow_ALperp_err_noq2,BtoKstarlow_ARperp_err_noq2,BtoKstarlow_ALpar_err_noq2,BtoKstarlow_ARpar_err_noq2,BtoKstarlow_AL0_err_noq2,BtoKstarlow_AR0_err_noq2,BtoKstarlow_At_err_noq2,BtoKstarlow_AS_err_noq2;
	double complex BtoKstarlow_ALperp_err_q2,BtoKstarlow_ARperp_err_q2,BtoKstarlow_ALpar_err_q2,BtoKstarlow_ARpar_err_q2,BtoKstarlow_AL0_err_q2,BtoKstarlow_AR0_err_q2,BtoKstarlow_At_err_q2,BtoKstarlow_AS_err_q2;
	double complex BtoKstarhigh_ALperp_err,BtoKstarhigh_ARperp_err,BtoKstarhigh_ALpar_err,BtoKstarhigh_ARpar_err,BtoKstarhigh_AL0_err,BtoKstarhigh_AR0_err,BtoKstarhigh_At_err,BtoKstarhigh_AS_err;

	/* B -> K hadronic uncertainties */
	double hadrerrBK;
	double complex BtoKlow_FV_err_noq2,BtoKlow_FA_err_noq2,BtoKlow_FS_err_noq2,BtoKlow_FP_err_noq2;
	double complex BtoKlow_FV_err_q2,BtoKlow_FA_err_q2,BtoKlow_FS_err_q2,BtoKlow_FP_err_q2;
	double complex BtoKhigh_FV_err,BtoKhigh_FA_err,BtoKhigh_FS_err,BtoKhigh_FP_err;	

	/* Bs -> phi hadronic uncertainties */
	double hadrerrBsphi;
	double complex Bstophilow_ALperp_err_noq2,Bstophilow_ARperp_err_noq2,Bstophilow_ALpar_err_noq2,Bstophilow_ARpar_err_noq2,Bstophilow_AL0_err_noq2,Bstophilow_AR0_err_noq2,Bstophilow_At_err_noq2,Bstophilow_AS_err_noq2;
	double complex Bstophilow_ALperp_err_q2,Bstophilow_ARperp_err_q2,Bstophilow_ALpar_err_q2,Bstophilow_ARpar_err_q2,Bstophilow_AL0_err_q2,Bstophilow_AR0_err_q2,Bstophilow_At_err_q2,Bstophilow_AS_err_q2;
	double complex Bstophihigh_ALperp_err,Bstophihigh_ARperp_err,Bstophihigh_ALpar_err,Bstophihigh_ARpar_err,Bstophihigh_AL0_err,Bstophihigh_AR0_err,Bstophihigh_At_err,Bstophihigh_AS_err;

	/* B -> K* power correction implementation */
	int BKstar_implementation;
	
	double complex hplus0,hminus0,hplus1,hminus1,hplus2,hminus2,hzero0,hzero1,hzero2; /* hadronic parameters (nonfactorisable power corrections) */
	
	double real_alpha_perp0,real_alpha_perp1,real_alpha_perp2;
	double real_alpha_par0,real_alpha_par1,real_alpha_par2;
	double real_alpha_zero0,real_alpha_zero1;
	double imag_alpha_perp0,imag_alpha_perp1,imag_alpha_perp2;
	double imag_alpha_par0,imag_alpha_par1,imag_alpha_par2;
	double imag_alpha_zero0,imag_alpha_zero1;
	
	double DeltaC9_M1_q2bar,r1_M1,r2_M1;
	double DeltaC9_M2_q2bar,r1_M2,r2_M2;
	double DeltaC9_M3_q2bar,r1_M3,r2_M3;
	
	/* B->K* mu mu - likelihood or method of moments */
	int likelihoodBKstarmumu;

	/* B -> D parameters */
	double Delta_BD,rho_D2_BD,V1_1_BD;
	
	/* B -> D* parameters */
	double Delta_BDstar,rho_Dstar2_BDstar,R1_1_BDstar,R2_1_BDstar,R3_1_BDstar,V1_1_BDstar,hA1_1_BDstar;
			
	/* NP contributions to Wilson coefficients */
	double complex deltaC[21],deltaCp[21],deltaCQ[5],deltaCQp[5];
				
	/* Decay widths */
	int widthcalc; /* 0=none, 1=hdecay, 2=feynhiggs */
	double width_h0,width_H0,width_A0,width_H,width_Z,width_W,width_top,width_H03,width_A02;
	double width_gluino,width_t1,width_t2,width_b1,width_b2,width_ul,width_ur,width_dl,width_dr;
	double width_cl,width_cr,width_sl,width_sr,width_el,width_er,width_ml,width_mr,width_tau1,width_tau2,width_gravitino;
	double width_nuel,width_numl,width_nutaul,width_c1,width_c2,width_o1,width_o2,width_o3,width_o4,width_o5;
				
	/* 2HDM */
	int THDM_model;
	double lambda_u[4][4],lambda_d[4][4],lambda_l[4][4];

	/* NMSSMTools */
	int NMSSMcoll,NMSSMtheory,NMSSMups1S,NMSSMetab1S;
        
    /* SDECAY */
	double BRtbW,BRtbH,BRtt1o1,BRtt1o2,BRtt1o3,BRtt1o4,BRtt2o1,BRtt2o2,BRtt2o3,BRtt2o4;
	double BRgluinot1tbar,BRgluinot1bart,BRgluinodldbar,BRgluinodlbard,BRgluinodrdbar,BRgluinodrbard,BRgluinoulubar,BRgluinoulbaru,BRgluinourubar,BRgluinourbaru,BRgluinoslsbar,BRgluinoslbars,BRgluinosrsbar,BRgluinosrbars,BRgluinoclcbar,BRgluinoclbarc,BRgluinocrcbar,BRgluinocrbarc,BRgluinob1bbar,BRgluinob1barb,BRgluinob2bbar,BRgluinob2barb,BRgluinot2tbar,BRgluinot2bart,BRgluinoo1g,BRgluinoo2g,BRgluinoo3g,BRgluinoo4g,BRgluinoo1ddbar,BRgluinoo1uubar,BRgluinoo1ssbar,BRgluinoo1ccbar,BRgluinoo1bbbar,BRgluinoo1ttbar,BRgluinoo2ddbar,BRgluinoo2uubar,BRgluinoo2ssbar,BRgluinoo2ccbar,BRgluinoo2bbbar,BRgluinoo2ttbar,BRgluinoo3ddbar,BRgluinoo3uubar,BRgluinoo3ssbar,BRgluinoo3ccbar,BRgluinoo3bbbar,BRgluinoo3ttbar,BRgluinoo4ddbar,BRgluinoo4uubar,BRgluinoo4ssbar,BRgluinoo4ccbar,BRgluinoo4bbbar,BRgluinoo4ttbar,BRgluinoc1dubar,BRgluinoc1udbar,BRgluinoc1scbar,BRgluinoc1csbar,BRgluinoc1btbar,BRgluinoc1tbbar,BRgluinoc2dubar,BRgluinoc2udbar,BRgluinoc2scbar,BRgluinoc2csbar,BRgluinoc2btbar,BRgluinoc2tbbar,BRgluinot1barW,BRgluinot1W,BRgluinot1barH,BRgluinot1H;
	double BRt1o1t,BRt1o2t,BRt1o3t,BRt1o4t,BRt1c1b,BRt1c2b,BRt1o1c,BRt1o1u,BRt1gluinoc;
	double BRt2o1t,BRt2o2t,BRt2o3t,BRt2o4t,BRt2c1b,BRt2c2b,BRt2t1h,BRt2t1Z,BRt2b1W,BRt2o1c,BRt2o1u,BRt2gluinoc;
	double BRb1o1b,BRb1o2b,BRb1o3b,BRb1o4b,BRb1c1t,BRb1c2t,BRb1gluinob,BRb1t1W;
	double BRb2o1b,BRb2o2b,BRb2o3b,BRb2o4b,BRb2c1t,BRb2c2t,BRb2gluinob,BRb2b1h,BRb2b1Z,BRb2t1W,BRb2t2W;
	double BRulo1u,BRulo2u,BRulo3u,BRulo4u,BRulc1d,BRulc2d,BRulgluinou;
	double BRuro1u,BRuro2u,BRuro3u,BRuro4u,BRurc1d,BRurgluinou,BRurc2d;
	double BRdlo1d,BRdlo2d,BRdlo3d,BRdlo4d,BRdlc1u,BRdlc2u,BRdlgluinod;
	double BRdro1d,BRdro2d,BRdro3d,BRdro4d,BRdrgluinod,BRdrc1u,BRdrc2u;
	double BRclo1c,BRclo2c,BRclo3c,BRclo4c,BRclc1s,BRclc2s,BRclgluinoc;
	double BRcro1c,BRcro2c,BRcro3c,BRcro4c,BRcrc1s,BRcrgluinoc,BRcrc2s;
	double BRslo1s,BRslo2s,BRslo3s,BRslo4s,BRslc1c,BRslc2c,BRslgluinos;
	double BRsro1s,BRsro2s,BRsro3s,BRsro4s,BRsrgluinos,BRsrc1c,BRsrc2c;
	double BRelo1e,BRelo2e,BRelo3e,BRelo4e,BRelc1nue,BRelc2nue;
	double BRero1e,BRero2e,BRero3e,BRero4e,BRerc1nue,BRerc2nue;
	double BRmlo1m,BRmlo2m,BRmlo3m,BRmlo4m,BRmlc1num,BRmlc2num;
	double BRmro1m,BRmro2m,BRmro3m,BRmro4m,BRmrc1num,BRmrc2num;
	double BRtau1o1tau,BRtau1o2tau,BRtau1o3tau,BRtau1o4tau,BRtau1c1nutau,BRtau1c2nutau,BRtau1nutaulH,BRtau1nutaulW;
	double BRtau2o1tau,BRtau2o2tau,BRtau2o3tau,BRtau2o4tau,BRtau2c1nutau,BRtau2c2nutau,BRtau2tau1h,BRtau2tau1Z,BRtau2nutaulH,BRtau2nutaulW,BRtau2tau1H,BRtau2tau1A;
	double BRnuelo1nue,BRnuelo2nue,BRnuelo3nue,BRnuelo4nue,BRnuelc1e,BRnuelc2e;
	double BRnumlo1num,BRnumlo2num,BRnumlo3num,BRnumlo4num,BRnumlc1m,BRnumlc2m;
	double BRnutaulo1nutau,BRnutaulo2nutau,BRnutaulo3nutau,BRnutaulo4nutau,BRnutaulc1tau,BRnutaulc2tau,BRnutaultau1W,BRnutaultau1H,BRnutaultau2H,BRnutaultau2W;
	double BRc1o1W,BRc1tau1nutau,BRc1o1udbar,BRc1o1csbar,BRc1o1enue,BRc1o1mnum,BRc1o1taunutau,BRc1o2udbar,BRc1o2csbar,BRc1o2enue,BRc1o2mnum,BRc1o2taunutau,BRc1o3udbar,BRc1o3csbar,BRc1o3enue,BRc1o3mnum,BRc1o3taunutau,BRc1o4udbar,BRc1o4csbar,BRc1o4enue,BRc1o4mnum,BRc1o4taunutau,BRc1nuele,BRc1numlm,BRc1elnue,BRc1mlnum,BRc1tau2nutau,BRc1c1Z,BRc1c1h,BRc1nutaultau,BRc1o2W;
	double BRc2c1Z,BRc2o1W,BRc2o2W,BRc2c1h,BRc2nuele,BRc2numlm,BRc2nutaultau,BRc2elnue,BRc2mlnum,BRc2tau1nutau,BRc2tau2nutau;
	double BRo2o1Z,BRo2o1h,BRo2tau1taubar,BRo2tau1bartau,BRo2o1gamma,BRo2o1ubaru,BRo2o1dbard,BRo2o1cbarc,BRo2o1sbars,BRo2o1bbarb,BRo2o1tbart,BRo2o1ebare,BRo2o1mbarm,BRo2o1taubartau,BRo2o1nuebarnue,BRo2o1numbarnum,BRo2o1nutaubarnutau,BRo2c1ubard,BRo2c1dbaru,BRo2c1cbars,BRo2c1sbarc,BRo2c1tbarb,BRo2c1bbart,BRo2c1nuebare,BRo2c1nueebar,BRo2c1numbarm,BRo2c1nummbar,BRo2c1nutaubartau,BRo2c1nutautaubar,BRo2c2ubard,BRo2c2dbaru,BRo2c2cbars,BRo2c2sbarc,BRo2c2tbarb,BRo2c2bbart,BRo2c2nuebare,BRo2c2nueebar,BRo2c2numbarm,BRo2c2nummbar,BRo2c2nutaubartau,BRo2c2nutautaubar,BRo2elebar,BRo2elbare,BRo2erebar,BRo2erbare,BRo2mlmbar,BRo2mlbarm,BRo2mrmbar,BRo2mrbarm,BRo2tau2taubar,BRo2tau2bartau,BRo2nuelnuebar,BRo2nuelbarnue,BRo2numlnumbar,BRo2numlbarnum,BRo2nutaulnutaubar,BRo2nutaulbarnutau;
	double BRo3o1Z,BRo3o2Z,BRo3c1W,BRo3c1barW,BRo3o1h,BRo3o2h,BRo3elebar,BRo3elbare,BRo3erebar,BRo3erbare,BRo3mlmbar,BRo3mlbarm,BRo3mrmbar,BRo3mrbarm,BRo3tau1taubar,BRo3tau1bartau,BRo3tau2taubar,BRo3tau2bartau,BRo3nuelnuebar,BRo3nuelbarnue,BRo3numlnumbar,BRo3numlbarnum,BRo3nutaulnutaubar,BRo3nutaulbarnutau,BRo3o1gamma,BRo3o2gamma;
	double BRo4o1Z,BRo4o2Z,BRo4c1W,BRo4c1barW,BRo4o1h,BRo4o2h,BRo4elebar,BRo4elbare,BRo4erebar,BRo4erbare,BRo4mlmbar,BRo4mlbarm,BRo4mrmbar,BRo4mrbarm,BRo4tau1taubar,BRo4tau1bartau,BRo4tau2taubar,BRo4tau2bartau,BRo4nuelnuebar,BRo4nuelbarnue,BRo4numlnumbar,BRo4numlbarnum,BRo4nutaulnutaubar,BRo4nutaulbarnutau,BRo4o1gamma,BRo4o2gamma,BRo4o3gamma;
	double BRo5o1Z,BRo5o2Z,BRo5c1W,BRo5c1barW,BRo5o1h,BRo5o2h,BRo5elebar,BRo5elbare,BRo5erebar,BRo5erbare,BRo5mlmbar,BRo5mlbarm,BRo5mrmbar,BRo5mrbarm,BRo5tau1taubar,BRo5tau1bartau,BRo5tau2taubar,BRo5tau2bartau,BRo5nuelnuebar,BRo5nuelbarnue,BRo5numlnumbar,BRo5numlbarnum,BRo5nutaulnutaubar,BRo5nutaulbarnutau,BRo5o1gamma,BRo5o2gamma,BRo5o3gamma;

    /* HDECAY & FeynHiggs */
	double mass_h0SM,width_h0SM;
	double mass_H0SM,width_H0SM;
	double mass_A0SM,width_A0SM;
	double BRh0bb_SM,BRh0tautau_SM,BRh0WW_SM,BRh0gg_SM,BRh0gaga_SM,BRh0ZZ_SM;
	double BRH0bb_SM,BRH0tautau_SM,BRH0WW_SM,BRH0gg_SM,BRH0gaga_SM,BRH0ZZ_SM;
	double BRA0bb_SM,BRA0tautau_SM,BRA0WW_SM,BRA0gg_SM,BRA0gaga_SM,BRA0ZZ_SM;

	double BRh0bb,BRh0tautau,BRh0WW,BRh0gg,BRh0gaga,BRh0ZZ;
	double BRH0bb,BRH0tautau,BRH0WW,BRH0gg,BRH0gaga,BRH0ZZ;
	double BRA0bb,BRA0tautau,BRA0WW,BRA0gg,BRA0gaga,BRA0ZZ;
	double BRh0mumu,BRh0ss,BRh0cc,BRh0tt,BRh0gaZ;
	double BRh0n1n2,BRh0n1n3,BRh0n1n4,BRh0n2n3,BRh0n2n4,BRh0c1c1,BRh0c1c2,BRh0n1n1,BRh0n2n2;
	double BRH0mumu,BRH0ss,BRH0cc,BRH0tt,BRH0gaZ,BRH0hZ;
	double BRH0n1n2,BRH0n1n3,BRH0n1n4,BRH0n2n3,BRH0n2n4,BRH0c1c1,BRH0c1c2,BRH0n1n1,BRH0n2n2,BRH0hh;
	double BRA0mumu,BRA0ss,BRA0cc,BRA0tt,BRA0gaZ,BRA0hZ;
	double BRA0n1n2,BRA0n1n3,BRA0n1n4,BRA0n2n3,BRA0n2n4,BRA0c1c1,BRA0c1c2,BRA0n1n1,BRA0n2n2,BRA0hh;
	double BRHmunu,BRHtaunu,BRHub,BRHus,BRHcs,BRHcb,BRHtb,BRHWh,BRHWA,BRHc1n1,BRHc1n2,BRHc1n3,BRHc1n4,BRHc2n1,BRHc2n2;
	double BRh0mumu_SM,BRh0ss_SM,BRh0cc_SM,BRh0tt_SM,BRh0gaZ_SM;
	double BRh0stau1stau1,BRh0stau1stau2,BRh0stau2stau2;
	double BRH0stau1stau1,BRH0stau1stau2,BRH0stau2stau2;
	double BRA0stau1stau1,BRA0stau1stau2,BRA0stau2stau2;
	double BRh0b1b1,BRh0b1b2,BRh0b2b2;
	double BRH0b1b1,BRH0b1b2,BRH0b2b2;
	double BRA0b1b1,BRA0b1b2,BRA0b2b2;
	
	double BRH03bb,BRH03tautau,BRH03WW,BRH03gg,BRH03gaga,BRH03ZZ;
	double BRA02bb,BRA02tautau,BRA02WW,BRA02gg,BRA02gaga,BRA02ZZ;
	double BRH03mumu,BRH03ss,BRH03cc,BRH03tt,BRH03gaZ,BRH03hZ;
	double BRH03n1n2,BRH03n1n3,BRH03n1n4,BRH03n2n3,BRH03n2n4,BRH03c1c1,BRH03c1c2,BRH03n1n1,BRH03n2n2,BRH03hh;
	double BRA02mumu,BRA02ss,BRA02cc,BRA02tt,BRA02gaZ,BRA02hZ;
	double BRA02n1n2,BRA02n1n3,BRA02n1n4,BRA02n2n3,BRA02n2n4,BRA02c1c1,BRA02c1c2,BRA02n1n1,BRA02n2n2,BRA02hh;
	double BRH03stau1stau1,BRH03stau1stau2,BRH03stau2stau2;
	double BRA02stau1stau1,BRA02stau1stau2,BRA02stau2stau2;
	double BRH03b1b1,BRH03b1b2,BRH03b2b2;
	double BRA02b1b1,BRA02b1b2,BRA02b2b2;
}
parameters;

/*--------------------------------------------------------------------*/

typedef struct flhaparam
/* structure containing all the scanned parameters from the FLHA file */
{
	double complex dC7[3],dC8[3],dC9[3],dC10[3],dC9e[3],dC10e[3];
	double complex dC7p[3],dC8p[3],dC9p[3],dC10p[3],dC9pe[3],dC10pe[3];
	double complex C7[3],C8[3],C9[3],C10[3],C9e[3],C10e[3];
	double complex C7p[3],C8p[3],C9p[3],C10p[3],C9pe[3],C10pe[3];
	double complex C7SM[3],C8SM[3],C9SM[3],C10SM[3],C9eSM[3],C10eSM[3];
	double complex C7pSM[3],C8pSM[3],C9pSM[3],C10pSM[3],C9peSM[3],C10peSM[3];
	double Q;
}
flhaparam;

/*--------------------------------------------------------------------*/
/* Prototypes */

/* isajet.c */
int isajet_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, char name[]);
int isajet_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, char name[]);
int isajet_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, char name[]);
int isajet_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, char name[]);
int isajet_mmamsb(double alpha, double m32, double tanb, double sgnmu, double mtop, char name[]);
int isajet_hcamsb(double alpha, double m32, double tanb, double sgnmu, double mtop, char name[]);

/* softsusy.c */
int softsusy_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_mssm(double m1, double m2, double m3, double tanb, double mA, double at, double ab, double atau, double mu, double mer, double mel, double mstaul, double mstaur, double mql, double mq3l, double mqur, double mqtr, double mqdr, double mqbr, double Q, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_slhain(char name_in[], char name_out[]);

/* leshouches.c */ 
int Les_Houches_Reader(char name[], struct parameters* param);
void Init_param(struct parameters* param);
void slha_adjust(struct parameters* param);
int test_slha(char name[]);
int SM_Decays_Reader(char name[], struct parameters* param);
int Higgs_Decays_Reader(char name[], struct parameters* param);
int SUSY_Decays_Reader(char name[], struct parameters* param);

/* alphas.c */
double alphas_running(double Q, double mtop, double mbot, struct parameters* param);

/* quark_masses.c */
double running_mass(double quark_mass, double Qinit, double Qfin,  double mtop, double mbot, struct parameters* param);
double mb_pole(struct parameters* param);
double mb_pole_1loop(struct parameters* param);
double mc_pole(struct parameters* param);
double mc_pole_1loop(struct parameters* param);
double mc_pole_3loops(struct parameters* param);
double mb_1S(struct parameters* param);
double mt_mt(struct parameters* param);
double mcmc_from_pole(double mcpole, int loop, struct parameters* param);
double mbmb_from_pole(double mbpole, int loop, struct parameters* param);

/* general.c */
double max(double x, double y);
double min(double x, double y);
double Ei1(double x);
double Ei2(double x);
double Ei3(double x);
double Ei(double x);
double complex polylog(int n, int m, double x);
double complex cd(double x, double y);
double complex hpl_base1(int i, double complex x);
double complex hpl_base2(int i1, int i2, double complex x);
double complex hpl_base3(int i1, int i2, int i3, double complex x);
double complex hpl_base4(int i1, int i2, int i3, int i4, double complex x);
double complex hpl1(int i, double complex x);
double complex hpl2(int i1, int i2, double complex x);
double complex hpl3(int i1, int i2, int i3, double complex x);
double complex hpl4(int i1, int i2, int i3, int i4, double complex x);
double Li2(double x);
double Li3(double x);
double Li4(double x);
double complex CLi2(double complex x);
double complex CLi3(double complex x);
double complex CLi4(double complex x);
double Cl2(double x);
double Cl3(double x);
double I0(double x);
double I1(double x);
double K0(double x);
double K1(double x);
double K2(double x);
double K3(double x);
double K4(double x);
double Lbessel(double x);
double Mbessel(double x);
double Nbessel(double x);
double expcor(double x);
double K0exp(double x,double z);
double K1exp(double x,double z);
double K2exp(double x,double z);
double kron(int x, int y);
double rand_gauss(void);
int test_integer(char name[]);
int test_file(char *name);
int invert_matrix(int n, double** in, double** out);
double rand_gauss();
double myrand(int type, int on);
int cholesky(int n, double** in, double** out);
double interpol_fromtable(double val, double table[], int dim, int order);
double simpson_noerr(double* f, double*x, int n1, int n2);
void simpson_rec(double (*func)(double), double* f, double *x, int* n, int n1, int n2, double* restot, double err);
double simpson(double (*func)(double), double a, double b, double err);
void braket(int n,double (*func)(double[]), double* a0, double*x0, double*b0, double xinit[], double xi[], double xlimmin[], double xlimmax[], double factor);
int brentmethod(int n, double (*func)(double[]),double xmin[],double xi[], double xlimmin[], double xlimmax[], double *fmin, double factor);
double brentmethod1Dfunc(double (*func)(double[]), double x, double xmini[], double xi[], int n);
int brentmethod1D(double (*func)(double[]), double a0, double x0,double b0, double *fmin, double *xmin, double factor, double xmini[], double xi[], int n);
void brentzero(double (*func)(double), double a, double b, double *fmin, double *xmin, double eps);
int powellaux(int n, double (*func)(double[]),double x0[], double xlimmin[], double xlimmax[],  double* fmin,  double xmin[], double ftol, double factor);
int powell(int n, double (*func)(double[]), double xlimmin[], double xlimmax[],  double* fmin,  double xmin[], double ftol, char* option);

/* wilson.c */
double H2(double x, double y);
double D3(double x);
double h10(double x);
double h20(double x);
double h30(double x);
double h40(double x);
double h50(double x);
double h60(double x);
double f20(double x);
double f30(double x,double y);
double f40(double x,double y);
double f50(double x,double y,double z);
double f60(double x,double y,double z);
double f70(double x, double y);
double f80(double x);
double f90(double w,double x,double y,double z);
double f100(double w,double x,double y,double z);
double f110(double x,double y);
double h11(double x,double y);
double h21(double x,double y);
double h31(double x,double y);
double h41(double x,double y);
double h51(double x,double y);
double h61(double x,double y);
double h71(double x,double y);
double f31(double x, double y);
double f41(double x,double y);
double f51(double x,double y);
double f61(double x,double y);
double f71(double x,double y,double z);
double f81(double x,double y,double z);
double f91(double x,double y,double z);
double f111(double x, double y);
double f121(double x, double y, double z);
double f131(double x, double y, double z);
double f141(double x, double y);
double f151(double x);
double f161(double x);
double f171(double x, double y);
double f181(double x, double y);
double f191(double x, double y);
double q11(double x,double y);
double q21(double x,double y);
double q31(double x,double y);
double q41(double x,double y);
double q51(double x,double y);
double q61(double x,double y);
double A0t(double x);
double A1t(double x, double l);
double B0t(double x);
double C0t(double x);
double D0t(double x);
double B1t(double x, double l);
double C1t(double x, double l);
double D1t(double x, double l);
double F0t(double x);
double F1t(double x,double l);
double E0t(double x);
double G1t(double x, double l);
double E1t(double x, double l);
double T(double x);
double F7_1(double x);
double F7_2(double x);
double F8_1(double x);
double F8_2(double x);
double B(double m1, double m2, double Q);
double G7H(double x, double lu, double ld);
double Delta7H(double x, double lu, double ld);
double G8H(double x, double lu, double ld);
double Delta8H(double x, double lu, double ld);
double EH(double x, double lu);
double G4H(double x, double lu);
double Delta4H(double x, double lu);
double G3H(double x, double lu);
double Delta3H(double x, double lu);
double C9llH0(double x, double y, double lu);
double D9H0(double x, double lu);
double C9llH1(double x, double y, double lu, double L);
double D9H1(double x, double lu, double L);
double C7t2mt(double x);
double C7c2MW(double x);
double C8t2mt(double x);
double C8c2MW(double x);
double C7H2lulumt(double r);
double C7H2ldlumt(double r);
double C7H2(double yt, double lu, double ld, double L);
double C8H2lulumt(double r);
double C8H2ldlumt(double r);
double C8H2(double yt, double lu, double ld, double L);
double C10Wt2mt(double x);
double C10Wc2MW(double x);
double C10Zt2mt(double x);
double C10Z2tri(double x);
double F0SP(double xt);double F1SP(double xt, double xH);
double F2SP(double xt, double xH);
double F3SP(double xt, double xH);
double F4SP(double xt, double xH);
double F5SP(double xt, double xH);
double F6SP(double xt, double xH);
double F7SP(double xt, double xH);
double F8SP(double xt, double xH);
double F9SP(double xt, double xH);
double F10SP(double xt, double xH);
double F11SP(double xt, double xH);
double F12SP(double xt, double xH);
double CSc_2HDM(double xH, double xt, double lu, double ld, double ll);
double CPc_2HDM(double xH, double xt, double lu, double ld, double ll, double sw2);
double epsilon_0(struct parameters* param);
double epsilon_2(struct parameters* param);
double epsilon_b(struct parameters* param);
double epsilon_bp(struct parameters* param);
double epsilon_0p(struct parameters* param);
double epsilon_1p(struct parameters* param);
void CW_calculator(int gen, double complex C0w[], double complex C1w[], double complex C2w[], double mu_W, struct parameters* param); 
void C_calculator_base1(double complex C0w[], double complex C1w[], double complex C2w[], double mu_W, double complex C0b[], double complex C1b[], double complex C2b[], double mu, struct parameters* param); 
void C_calculator_base2(double complex C0w[], double complex C1w[], double mu_W, double complex C0b[], double complex C1b[], double mu, struct parameters* param); 
void Cprime_calculator(int gen, double complex Cpb[], double complex CQpb[], double mu_W, double mu, struct parameters* param);
void CQ_calculator(int gen, double complex CQ0b[], double complex CQ1b[], double mu_W, double mu, struct parameters* param);

/* bsgamma.c */
double phi77(double delta);
double phi78(double delta);
double phi88(double delta, double b);
double G1(double t);
double phi22(double delta, double z);
double phi11(double delta, double z);
double phi12(double delta, double z);
double G2(double t);
double phi27(double delta, double z);
double phi17(double delta, double z);
double phi18(double delta, double z);
double phi28(double delta, double z);
double phi47(double delta);
double phi48(double delta);
double F2_nf(double z);
double F2_a(double z);
double F2_na(double z);
double phi77_2beta(double delta,double mu, struct parameters* param);
double phi77_2rem(double delta, struct parameters* param);
double Re_a(double z);
double Re_b(double z);
double C_BXlnu(struct parameters* param);
double bsgamma_Ecut(double complex C0[], double complex C1[], double complex C2[], double complex Cp[], double mu, double mu_W, double ETcut, struct parameters* param);
double bsgamma(double complex C0[], double complex C1[], double complex C2[], double complex Cp[], double mu, double mu_W, struct parameters* param);
double bsgamma_calculator(char name[]);

/* isospin.c */
double F_perp(double a1_perp, double a2_perp);
double X_perp1(double a1_perp, double a2_perp);
double X_perp2(double a1_perp, double a2_perp);
double complex G(double s, double x);
double complex G_perp(double s, double a1_perp, double a2_perp);
double complex H_perp(double s, double a1_par, double a2_par, double zeta3A, double zeta3V, double wA10, double deltatp, double deltatm);
double H8_perp(double a1_perp, double a2_perp);
double complex h(double u,double s);
double complex H2_perp(double s, double a1_perp, double a2_perp);
double delta0(double complex C0[],double complex C0_spec[],double complex C1[],double complex C1_spec[],double complex Cp[],struct parameters* param,double mub,double muspec, double lambda_h);
double delta0_calculator(char name[]);

/* excluded_masses.c */
int excluded_Higgs_mass_calculator(char name[]);
int excluded_Higgs_masses(struct parameters* param);
int excluded_SUSY_mass_calculator(char name[]);
int excluded_SUSY_masses(struct parameters* param);
int excluded_mass_calculator(char name[]);
int excluded_masses(struct parameters* param);
int charged_LSP_calculator(char name[]);
int charged_LSP(struct parameters* param);

/* gmuon.c */
double F1N(double x);
double F2N(double x);
double F1C(double x);
double F2C(double x);
double fPS(double x);
double fS(double x);
double fft(double x);
double muonI1(double a);
double muonI2(double a);
double muonI3(double a);
double muonf(double z);
double muong(double z);
double muon_gm2_calculator(char name[]);
double muon_gm2(struct parameters* param);

/* btaunu.c */
double Blnu(int gen, struct parameters* param);
double Blnu_calculator(int gen, char name[]);
double RBlnu(int gen, struct parameters* param);
double RBlnu_calculator(int gen, char name[]);
double Btaunu(struct parameters* param);
double Btaunu_calculator(char name[]);
double RBtaunu(struct parameters* param);
double RBtaunu_calculator(char name[]);

/* bdtaunu.c */
double dGammaBDlnu_dq2(int gen, int charge, double q2, double obs[][3], struct parameters* param);
double BRBDlnu(int gen, int charge, double q2min, double q2max, double obs[], struct parameters* param);
double BRBDlnu_calculator(int gen, int charge, double q2min, double q2max, double obs[], char name[]);
double BRBDlnu_full(int gen, int charge, double obs[], struct parameters* param);
double BRBDlnu_full_calculator(int gen, int charge, double obs[], char name[]);
double dGammaBDstarlnu_dq2(int gen, int charge, double q2, double obs[][3], struct parameters* param);
double BRBDstarlnu(int gen, int charge, double q2min, double q2max, double obs[], struct parameters* param);
double BRBDstarlnu_calculator(int gen, int charge, double q2min, double q2max, double obs[], char name[]);
double BRBDstarlnu_full(int gen, int charge, double obs[], struct parameters* param);
double BRBDstarlnu_full_calculator(int gen, int charge, double obs[], char name[]);
double BDtaunu(struct parameters* param);
double BDtaunu_calculator(char name[]);
double BDtaunu_BDenu(struct parameters* param);
double BDtaunu_BDenu_calculator(char name[]);
double BDstartaunu(struct parameters* param);
double BDstartaunu_calculator(char name[]);
double BDstartaunu_BDstarenu(struct parameters* param);
double BDstartaunu_BDstarenu_calculator(char name[]);

/* bsmumu.c */
double Bsll(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[],struct parameters* param, double mu_b);
double Bsmumu(double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[],struct parameters* param, double mu_b);
double Bsll_untag(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[],struct parameters* param, double mu_b);
double Bsmumu_untag(double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[],struct parameters* param, double mu_b);
double Bdll(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double Bdmumu(double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double Bsmumu_calculator(char name[]);
double Bsmumu_untag_calculator(char name[]);
double Bdmumu_calculator(char name[]);

/* kmunu.c */
double Kmunu_pimunu(struct parameters* param);
double Kmunu_pimunu_calculator(char name[]);
double Rmu23(struct parameters* param);
double Rmu23_calculator(char name[]);

/* dslnu.c */
double Dslnu(int gen, struct parameters* param);
double Dslnu_calculator(int gen, char name[]);
double Dstaunu(struct parameters* param);
double Dstaunu_calculator(char name[]);
double Dsmunu(struct parameters* param);
double Dsmunu_calculator(char name[]);

/* dmunu.c */
double Dlnu(int gen, struct parameters* param);
double Dlnu_calculator(int gen, char name[]);
double Dmunu(struct parameters* param);
double Dmunu_calculator(char name[]);

/* 2hdmc.c */
int thdmc_types(double l1, double l2, double l3, double l4, double l5, double l6, double l7, double m12_2, double tanb, int type, char name[]);

/* spheno.c */
int spheno_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int spheno_gmsb(double Lambda, double Mmess, double tanb, int N5, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int spheno_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);

/* suspect.c */
int suspect_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_gmsb(double Lambda, double Mmess, double tanb, int N5, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_mssm(double m1, double m2, double m3, double tanb, double mA, double at, double ab, double atau, double mu, double mer, double mel, double mstaul, double mstaur, double mql, double mq3l, double mqur, double mqtr, double mqdr, double mqbr, double Q, double mtop, double mbot, double alphas_mz, char name[]);

/* fleshouches.c */
void flha_generator(char name[], char name_output[]);
int FLHA_Reader(char name[], struct flhaparam* paramflha);

/* nmssmtools.c */
int nmssmtools_cnmssm(double m0, double m12, double tanb, double A0, double lambda, double AK, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int nmssmtools_nnuhm(double m0, double m12, double tanb, double A0, double MHDGUT, double MHUGUT, double lambda, double AK, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int nmssmtools_ngmsb(double Lambda, double Mmess, double tanb, int N5, double lambda, double AL, double Del_h, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int NMSSM_collider_excluded(char name[]);
int NMSSM_theory_excluded(char name[]);
int NMSSM_upsilon_excluded(char name[]);
int NMSSM_etab_excluded(char name[]);

/* bsll.c */
double complex g_bsll(double z, double s);
double Rcchad(double s, struct parameters* param);
double complex g_bsll_parametrized(double z, double s, struct parameters* param);
double glambda_bsll(double z);
double grho_bsll(double z);
double f_bsll(double x);
double h_bsll(double z);
double kappa_bsll(double x, double alfas);
double sigma_bsll(double s);
double sigma7_bsll(double s, double L);
double sigma9_bsll(double s);
double f7_bsll(double s);
double f9_bsll(double s);
double complex Gm1_bsll(double t);
double complex G0_bsll(double t);
double complex Di23_bsll(double s, double w, double z);
double complex Di27_bsll(double s, double w, double z);
double tau77_bsll(double s);
double tau99_bsll(double s);
double tau79_bsll(double s);
double tau710_bsll(double s);
double tau910_bsll(double s);
double tau22_bsll(double w, double s, double z);
double integ_tau22(double s, double z);
double tau78_bsll(double s);
double tau88_bsll(double s);
double tau89_bsll(double s);
double complex tau27_bsll(double w, double s, double z);
double complex integ_tau27(double s, double z);
double complex tau28_bsll(double w, double s, double z);
double complex integ_tau28(double s, double z);
double complex tau29_bsll(double w, double s, double z);
double complex integ_tau29(double s, double z);
double complex tau810_bsll(double s);
double complex dtau210_bsll(double w, double s, double z);
double complex tau210_bsll(double s, double z);
double complex F_bsll(double r);
double complex F87_bsll(double s, double L);
double F89_bsll(double s);
double w77em(double s, double L);
double w79em(double s, double L);
double w710em(double s, double L);
double w99em(double s, double L);
double w910em(double s, double L);
double w1010em(double s, double L);
double w22em(double s, double L, double mub);
double complex w27em(double s, double L, double mub);
double complex w29em(double s, double L, double mub);
double complex w210em(double s, double L, double a, double mub);
double dBR_BXsll_dshat(int gen, double shat,  double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dBR_BXsll_dshat_calculator(int gen, double shat, char name[]);
double BRBXsll(int gen, double smin, double smax, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBXsll_lowq2(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBXsll_highq2(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBXsmumu_lowq2_calculator(char name[]);
double BRBXsmumu_highq2_calculator(char name[]);
double BRBXstautau_highq2_calculator(char name[]);
double A_BXsll(int gen, double shat, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double A_BXsll_bin(int gen, double smin, double smax, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double A_BXsll_lowq2(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double A_BXsll_highq2(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double A_BXsmumu_lowq2_calculator(char name[]);
double A_BXsmumu_highq2_calculator(char name[]);
double A_BXsll_zero(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double A_BXsmumu_zero_calculator(char name[]);

/* bsll_extra.c */
double complex F17_bsll(double s, double z, double L);
double complex F27_bsll(double s, double z, double L);
double complex F19_bsll(double s, double z, double L);
double complex F29_bsll(double s, double z, double L);
double complex F19_bkll(double s, double z, double L);
double complex F29_bkll(double s, double z, double L);
double F_17re(double muh, double z, double sh, int maxpow);
double F_17re_cache(double muh, double z, double sh, int maxpow);
double F_17im(double muh, double z, double sh, int maxpow);
double F_17im_cache(double muh, double z, double sh, int maxpow);
double F_19re(double muh, double z, double sh, int maxpow);
double F_19re_cache(double muh, double z, double sh, int maxpow);
double F_19im(double muh, double z, double sh, int maxpow);
double F_19im_cache(double muh, double z, double sh, int maxpow);
double F_27re(double muh, double z, double sh, int maxpow);
double F_27re_cache(double muh, double z, double sh, int maxpow);
double F_27im(double muh, double z, double sh, int maxpow);
double F_27im_cache(double muh, double z, double sh, int maxpow);
double F_29re(double muh, double z, double sh, int maxpow);
double F_29re_cache(double muh, double z, double sh, int maxpow);
double F_29im(double muh, double z, double sh, int maxpow);
double F_29im_cache(double muh, double z, double sh, int maxpow);
double DeltaF_29re(double muh, double z, double sh, int maxpow);
double DeltaF_29re_cache(double muh, double z, double sh, int maxpow);
double DeltaF_29im(double muh, double z, double sh, int maxpow);
double DeltaF_29im_cache(double muh, double z, double sh, int maxpow);
double DeltaF_19re(double muh, double z, double sh, int maxpow);
double DeltaF_19re_cache(double muh, double z, double sh, int maxpow);
double DeltaF_19im(double muh, double z, double sh, int maxpow);
double DeltaF_19im_cache(double muh, double z, double sh, int maxpow);

/* bkstarll.c */
double complex h_bkll(double q2, double mq, double mu);
double phi_Kstar(double u, double a1, double a2);
double complex B0_bkll(double s, double mq);
double complex L1_bkll(double complex x);
double complex I1_bkll(double u, double mq, double q2, struct parameters* param);
double complex tperp_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param);
double complex tpar_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param);
double complex F27u(double shat, double l);
double complex F19u(double shat, double l);
double complex F29u(double shat, double l);
double dAI_BKstarmumu_dq2(double q2, double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu(double smin, double smax, double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_lowq2(double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_highq2(double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_zero(double complex C0b[], double complex C1b[], double complex C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_lowq2_calculator(char name[]);
double AI_BKstarmumu_highq2_calculator(char name[]);
double AI_BKstarmumu_zero_calculator(char name[]);
/* soft */
double dGamma_BKstarll_dq2_soft(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKstarll_dq2_soft_calculator(int gen, int charge, double q2, double obs[][3], char name[]);
double BRBKstarll_soft(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarll_lowq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarll_highq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_BKstarmumu_lowq2_soft_calculator(char name[],double obs[]);
double BRobs_BKstarmumu_highq2_soft_calculator(char name[],double obs[]);
/* full */
double dGamma_BKstarll_dq2_full(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKstarll_dq2_full_calculator(int gen, int charge, double q2, double obs[][3], char name[]);
double BRBKstarll_full(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarll_lowq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarll_highq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_BKstarmumu_lowq2_full_calculator(char name[], double obs[]);
double BRobs_BKstarmumu_highq2_full_calculator(char name[], double obs[]);
/* wrapper */
double dGamma_BKstarll_dq2(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKstarll_dq2_calculator(int gen, int charge, double q2, double obs[][3], char name[]);
double BRBKstarll(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarll_lowq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarll_highq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_BKstarmumu_lowq2_calculator(char name[], double obs[]);
double BRobs_BKstarmumu_highq2_calculator(char name[], double obs[]);

/* bkstargamma.c */
double BR_BKstargamma(int charge, double complex C0b[], double complex C1b[], double complex C2b[], double complex Cpb[], struct parameters* param, double mu_b);
double BR_BKstargamma_calculator(int charge, char name[]);

/* bkll.c */
/* soft */
double dGamma_BKll_dq2_soft(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKll_dq2_soft_calculator(int gen, int charge, double q2, double obs[][3], char name[]);
double BRBKll_soft(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double RK_BKll_soft(int charge, double smin, double smax, double complex C0be[], double complex C1be[], double complex C2be[], double complex CQ0be[], double complex CQ1be[], double complex Cpbe[], double complex CQpbe[], double complex C0bmu[], double complex C1bmu[], double complex C2bmu[], double complex CQ0bmu[], double complex CQ1bmu[], double complex Cpbmu[], double complex CQpbmu[], struct parameters* param, double mu_b);
double BRBKll_lowq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKll_highq2_soft(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_BKmumu_lowq2_soft_calculator(char name[], double obs[]);
double BRobs_BKmumu_highq2_soft_calculator(char name[], double obs[]);
double RK_BKll_soft_calculator(int charge, double smin, double smax,char name[]);
/* full */
double dGamma_BKll_dq2_full(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKll_dq2_full_calculator(int gen, int charge, double q2, double obs[][3], char name[]);
double BRBKll_full(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double RK_BKll_full(int charge, double smin, double smax, double complex C0be[], double complex C1be[], double complex C2be[], double complex CQ0be[], double complex CQ1be[], double complex Cpbe[], double complex CQpbe[], double complex C0bmu[], double complex C1bmu[], double complex C2bmu[], double complex CQ0bmu[], double complex CQ1bmu[], double complex Cpbmu[], double complex CQpbmu[], struct parameters* param, double mu_b);
double BRBKll_lowq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKll_highq2_full(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_BKmumu_lowq2_full_calculator(char name[], double obs[]);
double BRobs_BKmumu_highq2_full_calculator(char name[], double obs[]);
double RK_BKll_full_calculator(int charge, double smin, double smax,char name[]);
/* wrapper */
double dGamma_BKll_dq2(int gen, int charge, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKll_dq2_calculator(int gen, int charge, double q2, double obs[][3], char name[]);
double BRBKll(int gen, int charge, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double RK_BKll(int charge, double smin, double smax, double complex C0be[], double complex C1be[], double complex C2be[], double complex CQ0be[], double complex CQ1be[], double complex Cpbe[], double complex CQpbe[], double complex C0bmu[], double complex C1bmu[], double complex C2bmu[], double complex CQ0bmu[], double complex CQ1bmu[], double complex Cpbmu[], double complex CQpbmu[], struct parameters* param, double mu_b);
double BRBKll_lowq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKll_highq2(int gen, int charge, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_BKmumu_lowq2_calculator(char name[], double obs[]);
double BRobs_BKmumu_highq2_calculator(char name[], double obs[]);
double RK_BKll_calculator(int charge, double smin, double smax,char name[]);

/* bsphill.c */
/* soft */
double dGamma_Bsphill_dq2_soft(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_Bsphill_dq2_soft_calculator(int gen, double q2, double obs[][3], char name[]);
double BRBsphill_soft(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBsphill_lowq2_soft(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBsphill_highq2_soft(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_Bsphill_lowq2_soft_calculator(char name[],double obs[]);
double BRobs_Bsphill_highq2_soft_calculator(char name[],double obs[]);
/* full.c */
double dGamma_Bsphill_dq2_full(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_Bsphill_dq2_full_calculator(int gen, double q2, double obs[][3], char name[]);
double BRBsphill_full(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBsphill_lowq2_full(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBsphill_highq2_full(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_Bsphill_lowq2_full_calculator(char name[],double obs[]);
double BRobs_Bsphill_highq2_full_calculator(char name[],double obs[]);
/* wrapper */
double dGamma_Bsphill_dq2(int gen, double q2, double obs[][3], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_Bsphill_dq2_calculator(int gen, double q2, double obs[][3], char name[]);
double BRBsphill(int gen, double smin, double smax, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBsphill_lowq2(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBsphill_highq2(int gen, double obs[], double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRobs_Bsphill_lowq2_calculator(char name[],double obs[]);
double BRobs_Bsphill_highq2_calculator(char name[],double obs[]);

/* chi2.c */
int observables(int ke, double complex deltaC[], double complex deltaCp[], double complex deltaCQ[], double complex deltaCQp[], char* names[], double values[], int onoff[], double values_ref[], struct parameters* param);
void activate_observables(int nbobs0, char* names0[], int nbobsin, char** namesin, int onoff[]);
void correlation_matrices(int fullFF, double **correlation_matrix_FFKstar, double **correlation_matrix_FFphi, double **correlation_matrix_FFK);
void experimental_covariance(int likelihood, char* namesin[], int nbobsin, double *central_exp, double *errors_exp, double **correlations);
int get_covariance(double ***covariance_th, double ***covariance_exp, double **central_exp, char* namesin[], int nbobsin, double complex deltaC[], double complex deltaCp[], double complex deltaCQ[], double complex deltaCQp[], struct parameters* param, int *nbobs_out);
void get_covtot(double ***covariance_th, double ***covariance_exp, double ***covariance_tot, int nbobs);
int reduce_covariance(double ***covariance_in, char* namesin[], int nbobsin, double ***covariance_out, char* namesout[], int nbobsout);
int reduce_values(double **values_in, char* namesin[], int nbobsin, double **values_out, char* namesout[], int nbobsout);
int get_invcovtot(double ***covariance_tot, double ***inv_cov_tot, int nbobs);
void print_covariance(char filename[], double **covariance_in, double *central_in, char* namesin[], int nbobsin);
int read_nbobs(char filename[]);
int read_covariance(char filename[], double ***covariance, double **central, char* names[], int *nbobs);
int get_predictions(char* namesin[], int nbobsin, double** predictions, double complex deltaC[], double complex deltaCp[], double complex deltaCQ[], double complex deltaCQp[], struct parameters* param);
double get_chi2(double **inv_cov_tot, double *predictions, double *central_exp, int nbobs);
double chi2(char* namesin[], int nbobsin, struct parameters* param, int *nbobs);

/* zwidths.c */
double Zwidth_neu1neu1(struct parameters* param);
double Zwidth_invisible(struct parameters* param);
double Zwidth_b1b1(struct parameters* param);
double Zwidth_neu1neu1_calculator(char name[]);
double Zwidth_invisible_calculator(char name[]);
double Zwidth_b1b1_calculator(char name[]);

#endif
