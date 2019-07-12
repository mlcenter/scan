#include "include.h"

void Init_param(struct parameters* param)
/* initializes the "param" structure by initializing the parameters with either 0 or a value from the PDG */
{
	int ie,je;
	
	param->SM=0;
	param->widthcalc=0; /* 0 = automatic choice, 1 = Hdecay, 2 = FeynHiggs, 3 = FeynHiggs Tree */
	param->model=-3; /* this parameter is used to test whether the scan of the SLHA file succeeds. */
	param->generator=0;
	param->Q=0.;
	param->m0=0.;
	param->m12=0.;
	param->tan_beta=0.;
	param->sign_mu=0.;
	param->A0=0.;
	param->mass_W=0.;
	param->Lambda=0.;
	param->Mmess=0.;
	param->N5=0.;
	param->cgrav=0.;
	param->m32=0.;
	param->mass_Z=0.;
	param->mass_b=0.;
	param->mass_top_pole=0.;
	param->mass_tau_pole=0.;
	param->inv_alpha_em=0.;
	param->alphas_MZ=0.;
	param->alpha=0.;
	param->Gfermi=0.;
	param->GAUGE_Q=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++)
	{
		param->charg_Umix[ie][je]=0.;
		param->charg_Vmix[ie][je]=0.;
		param->stop_mix[ie][je]=0.;
		param->sbot_mix[ie][je]=0.;
		param->stau_mix[ie][je]=0.;
	}
	for(ie=1;ie<=5;ie++) for(je=1;je<=5;je++) param->neut_mix[ie][je]=0.;
	for(ie=1;ie<=5;ie++) param->mass_neut[ie]=0.;
	for(ie=1;ie<=3;ie++) param->yub[ie]=param->yut[ie]=param->yutau[ie]=0.;
		
	param->Min=0.;
	param->M1_Min=0.;
	param->M2_Min=0.;
	param->M3_Min=0.;
	param->At_Min=0.;
	param->Ab_Min=0.;
	param->Atau_Min=0.;
	param->M2H1_Min=0.;
	param->M2H2_Min=0.;
	param->mu_Min=0.;
	param->M2A_Min=0.;
	param->tb_Min=0.;
	param->mA_Min=0.;
	param->MeL_Min=0.;
	param->MmuL_Min=0.;
	param->MtauL_Min=0.;
	param->MeR_Min=0.;
	param->MmuR_Min=0.;
	param->MtauR_Min=0.;
	param->MqL1_Min=0.;
	param->MqL2_Min=0.;
	param->MqL3_Min=0.;
	param->MuR_Min=0.;
	param->McR_Min=0.;
	param->MtR_Min=0.;
	param->MdR_Min=0.;
	param->MsR_Min=0.;
	param->MbR_Min=0.;
	param->N51=0.;
	param->N52=0.;
	param->N53=0.;
	param->M2H1_Q=0.;
	param->M2H2_Q=0.;
	param->mass_h0=0.;
	param->mass_H0=0.;
	param->mass_A0=0.;
	param->mass_H=0.;
	param->mass_dnl=0.;
	param->mass_upl=0.;
	param->mass_stl=0.;
	param->mass_chl=0.;
	param->mass_b1=0.;
	param->mass_t1=0.;
	param->mass_el=0.;
	param->mass_nuel=0.;
	param->mass_mul=0.;
	param->mass_numl=0.;
	param->mass_tau1=0.;
	param->mass_nutl=0.;
	param->mass_gluino=0.;
	param->mass_cha1=0.;
	param->mass_cha2=0.;
	param->mass_dnr=0.;
	param->mass_upr=0.;
	param->mass_str=0.;
	param->mass_chr=0.;
	param->mass_b2=0.;
	param->mass_t2=0.;
	param->mass_er=0.;
	param->mass_mur=0.;
	param->mass_tau2=0.;
	param->gp=0.;
	param->g2=0.;
	param->gp_Q=0.;
	param->g2_Q=0.;
	param->g3_Q=0.;
	param->YU_Q=0.;
	param->YD_Q=0.;
	param->YE_Q=0.;
	param->HMIX_Q=0.;
	param->mu_Q=0.;
	param->tanb_GUT=0.;
	param->Higgs_VEV=0.;
	param->mA2_Q=0.;
	param->MSOFT_Q=0.;
	param->M1_Q=0.;
	param->M2_Q=0.;
	param->M3_Q=0.;
	param->MeL_Q=0.;
	param->MmuL_Q=0.;
	param->MtauL_Q=0.;
	param->MeR_Q=0.;
	param->MmuR_Q=0.;
	param->MtauR_Q=0.;
	param->MqL1_Q=0.;
	param->MqL2_Q=0.;
	param->MqL3_Q=0.;
	param->MuR_Q=0.;
	param->McR_Q=0.;
	param->MtR_Q=0.;
	param->MdR_Q=0.;
	param->MsR_Q=0.;
	param->MbR_Q=0.;
	param->AU_Q=0.;
	param->A_u=0.;
	param->A_c=0.;
	param->A_t=0.;
	param->AD_Q=0.;
	param->A_d=0.;
	param->A_s=0.;
	param->A_b=0.;
	param->AE_Q=0.;
	param->A_e=0.;
	param->A_mu=0.;
	param->A_tau=0.;
	param->mass_graviton=0.;
	param->mass_gravitino=0.;
	param->mass_nuer=0.;
	param->mass_numr=0.;
	param->mass_nutr=0.;
	param->mass_t=0.;
	param->mass_tau=0.;
	param->mass_gluon=0.;
	param->mass_nue=0.;
	param->mass_num=0.;
	param->mass_nut=0.;
	param->mass_photon=0.;
	param->mass_Z0=0.;

	/* SLHA2 */
	param->NMSSM=0;
	param->RV=0;
	param->CPV=0;
	param->FV=0;
	param->CKM_lambda=0.;
	param->CKM_A=0.;
	param->CKM_rhobar=0.;
	param->CKM_etabar=0.;
	param->PMNS_theta12=0.;
	param->PMNS_theta23=0.;
	param->PMNS_theta13=0.;
	param->PMNS_delta13=0.;
	param->PMNS_alpha1=0.;
	param->PMNS_alpha2=0.;
	param->lambdaNMSSM_Min=0.;
	param->kappaNMSSM_Min=0.;
	param->AlambdaNMSSM_Min=0.;
	param->AkappaNMSSM_Min=0.;
	param->lambdaSNMSSM_Min=0.;
	param->xiFNMSSM_Min=0.;
	param->xiSNMSSM_Min=0.;
	param->mupNMSSM_Min=0.;
	param->mSp2NMSSM_Min=0.;
	param->mS2NMSSM_Min=0.;
	param->mass_H03=0.;
	param->mass_A02=0.;
	param->NMSSMRUN_Q=0.;
	param->lambdaNMSSM=0.;
	param->kappaNMSSM=0.;
	param->AlambdaNMSSM=0.;
	param->AkappaNMSSM=0.;
	param->lambdaSNMSSM=0.;
	param->xiFNMSSM=0.;
	param->xiSNMSSM=0.;
	param->mupNMSSM=0.;
	param->mSp2NMSSM=0.;
	param->mS2NMSSM=0.;
	param->PMNSU_Q=0.;
	param->CKM_Q=0.;
	param->IMCKM_Q=0.;
	param->MSE2_Q=0.;
	param->MSU2_Q=0.;
	param->MSD2_Q=0.;
	param->MSL2_Q=0.;
	param->MSQ2_Q=0.;
	param->TU_Q=0.;
	param->TD_Q=0.;
	param->TE_Q=0.;
	param->NMSSMcoll=0;
	param->NMSSMtheory=0;
	param->NMSSMups1S=0;
	param->NMSSMetab1S=0;
	
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		param->sU_mix[ie][je]=0.;
		param->sD_mix[ie][je]=0.;
		param->sE_mix[ie][je]=0.;
	}
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		param->H0_mix[ie][je]=0.;
		param->A0_mix[ie][je]=0.;
		param->sNU_mix[ie][je]=0.;
		param->sCKM_msq2[ie][je]=0.;
		param->sCKM_msl2[ie][je]=0.;
		param->sCKM_msd2[ie][je]=0.;
		param->sCKM_msu2[ie][je]=0.;
		param->sCKM_mse2[ie][je]=0.;
		param->PMNS_U[ie][je]=0.;
		param->CKM[ie][je]=0.;
		param->IMCKM[ie][je]=0.;
		param->TU[ie][je]=0.;
		param->TD[ie][je]=0.;
		param->TE[ie][je]=0.;
	}
	
	/* non-SLHA*/
	param->mass_c_pole=0.;
	param->mass_b_1S=0.;
	param->mass_b_pole=0.;
	param->mtmt=0.;
	param->Lambda3=0.;
	param->Lambda4=0.;
	param->Lambda5=0.;
	param->Lambda6=0.;
	param->alphasMZ_Lambda3=0.;
	param->alphasMZ_Lambda4=0.;
	param->alphasMZ_Lambda5=0.;
	param->alphasMZ_Lambda6=0.;
	
	/* widths */
	param->width_h0=0.;
	param->width_H0=0.;
	param->width_A0=0.;
	param->width_H=0.;
	param->width_H03=0.;
	param->width_A02=0.;
	param->width_gluino=0.1;
	param->width_t1=0.1;
	param->width_t2=0.1;
	param->width_b1=0.1;
	param->width_b2=0.1;
	param->width_ul=0.1;
	param->width_ur=0.1;
	param->width_dl=0.1;
	param->width_dr=0.1;;
	param->width_cl=0.1;
	param->width_cr=0.1;
	param->width_sl=0.1;
	param->width_sr=0.1;
	param->width_el=0.1;
	param->width_er=0.1;
	param->width_ml=0.1;
	param->width_mr=0.1;
	param->width_tau1=0.1;
	param->width_tau2=0.1;
	param->width_nuel=0.1;
	param->width_numl=0.1;
	param->width_nutaul=0.1;
	param->width_c1=0.1;
	param->width_c2=0.1;
	param->width_o1=0.1;
	param->width_o2=0.1;
	param->width_o3=0.1;
	param->width_o4=0.1;
	param->width_o5=0.1;
		
	/* BR */
	param->BRtbW=param->BRtbH=param->BRtt1o1=param->BRtt1o2=param->BRtt1o3=param->BRtt1o4=param->BRtt2o1=param->BRtt2o2=param->BRtt2o3=param->BRtt2o4=param->BRgluinot1tbar=param->BRgluinot1bart=param->BRgluinodldbar=param->BRgluinodlbard=param->BRgluinodrdbar=param->BRgluinodrbard=param->BRgluinoulubar=param->BRgluinoulbaru=param->BRgluinourubar=param->BRgluinourbaru=param->BRgluinoslsbar=param->BRgluinoslbars=param->BRgluinosrsbar=param->BRgluinosrbars=param->BRgluinoclcbar=param->BRgluinoclbarc=param->BRgluinocrcbar=param->BRgluinocrbarc=param->BRgluinob1bbar=param->BRgluinob1barb=param->BRgluinob2bbar=param->BRgluinob2barb=param->BRgluinot2tbar=param->BRgluinot2bart=param->BRgluinoo1g=param->BRgluinoo2g=param->BRgluinoo3g=param->BRgluinoo4g=param->BRgluinoo1ddbar=param->BRgluinoo1uubar=param->BRgluinoo1ssbar=param->BRgluinoo1ccbar=param->BRgluinoo1bbbar=param->BRgluinoo1ttbar=param->BRgluinoo2ddbar=param->BRgluinoo2uubar=param->BRgluinoo2ssbar=param->BRgluinoo2ccbar=param->BRgluinoo2bbbar=param->BRgluinoo2ttbar=param->BRgluinoo3ddbar=param->BRgluinoo3uubar=param->BRgluinoo3ssbar=param->BRgluinoo3ccbar=param->BRgluinoo3bbbar=param->BRgluinoo3ttbar=param->BRgluinoo4ddbar=param->BRgluinoo4uubar=param->BRgluinoo4ssbar=param->BRgluinoo4ccbar=param->BRgluinoo4bbbar=param->BRgluinoo4ttbar=param->BRgluinoc1dubar=param->BRgluinoc1udbar=param->BRgluinoc1scbar=param->BRgluinoc1csbar=param->BRgluinoc1btbar=param->BRgluinoc1tbbar=param->BRgluinoc2dubar=param->BRgluinoc2udbar=param->BRgluinoc2scbar=param->BRgluinoc2csbar=param->BRgluinoc2btbar=param->BRgluinoc2tbbar=param->BRgluinot1barW=param->BRgluinot1W=param->BRgluinot1barH=param->BRgluinot1H=param->BRt1o1t=param->BRt1o2t=param->BRt1o3t=param->BRt1o4t=param->BRt1c1b=param->BRt1c2b=param->BRt1o1c=param->BRt1o1u=param->BRt1gluinoc=param->BRt2o1t=param->BRt2o2t=param->BRt2o3t=param->BRt2o4t=param->BRt2c1b=param->BRt2c2b=param->BRt2t1h=param->BRt2t1Z=param->BRt2b1W=param->BRt2o1c=param->BRt2o1u=param->BRt2gluinoc=param->BRb1o1b=param->BRb1o2b=param->BRb1o3b=param->BRb1o4b=param->BRb1c1t=param->BRb1c2t=param->BRb1gluinob=param->BRb1t1W=param->BRb2o1b=param->BRb2o2b=param->BRb2o3b=param->BRb2o4b=param->BRb2c1t=param->BRb2c2t=param->BRb2gluinob=param->BRb2b1h=param->BRb2b1Z=param->BRb2t1W=param->BRb2t2W=param->BRulo1u=param->BRulo2u=param->BRulo3u=param->BRulo4u=param->BRulc1d=param->BRulc2d=param->BRulgluinou=param->BRuro1u=param->BRuro2u=param->BRuro3u=param->BRuro4u=param->BRurc1d=param->BRurgluinou=param->BRurc2d=param->BRdlo1d=param->BRdlo2d=param->BRdlo3d=param->BRdlo4d=param->BRdlc1u=param->BRdlc2u=param->BRdlgluinod=param->BRdro1d=param->BRdro2d=param->BRdro3d=param->BRdro4d=param->BRdrgluinod=param->BRdrc1u=param->BRdrc2u=param->BRclo1c=param->BRclo2c=param->BRclo3c=param->BRclo4c=param->BRclc1s=param->BRclc2s=param->BRclgluinoc=param->BRcro1c=param->BRcro2c=param->BRcro3c=param->BRcro4c=param->BRcrc1s=param->BRcrgluinoc=param->BRcrc2s=param->BRslo1s=param->BRslo2s=param->BRslo3s=param->BRslo4s=param->BRslc1c=param->BRslc2c=param->BRslgluinos=param->BRsro1s=param->BRsro2s=param->BRsro3s=param->BRsro4s=param->BRsrgluinos=param->BRsrc1c=param->BRsrc2c=param->BRelo1e=param->BRelo2e=param->BRelo3e=param->BRelo4e=param->BRelc1nue=param->BRelc2nue=param->BRero1e=param->BRero2e=param->BRero3e=param->BRero4e=param->BRerc1nue=param->BRerc2nue=param->BRmlo1m=param->BRmlo2m=param->BRmlo3m=param->BRmlo4m=param->BRmlc1num=param->BRmlc2num=param->BRmro1m=param->BRmro2m=param->BRmro3m=param->BRmro4m=param->BRmrc1num=param->BRmrc2num=param->BRtau1o1tau=param->BRtau1o2tau=param->BRtau1o3tau=param->BRtau1o4tau=param->BRtau1c1nutau=param->BRtau1c2nutau=param->BRtau1nutaulH=param->BRtau1nutaulW=param->BRtau2o1tau=param->BRtau2o2tau=param->BRtau2o3tau=param->BRtau2o4tau=param->BRtau2c1nutau=param->BRtau2c2nutau=param->BRtau2tau1h=param->BRtau2tau1Z=param->BRtau2nutaulH=param->BRtau2nutaulW=param->BRtau2tau1H=param->BRtau2tau1A=param->BRnuelo1nue=param->BRnuelo2nue=param->BRnuelo3nue=param->BRnuelo4nue=param->BRnuelc1e=param->BRnuelc2e=param->BRnumlo1num=param->BRnumlo2num=param->BRnumlo3num=param->BRnumlo4num=param->BRnumlc1m=param->BRnumlc2m=param->BRnutaulo1nutau=param->BRnutaulo2nutau=param->BRnutaulo3nutau=param->BRnutaulo4nutau=param->BRnutaulc1tau=param->BRnutaulc2tau=param->BRnutaultau1W=param->BRnutaultau1H=param->BRnutaultau2H=param->BRnutaultau2W=param->BRc1o1W=param->BRc1tau1nutau=param->BRc1o1udbar=param->BRc1o1csbar=param->BRc1o1enue=param->BRc1o1mnum=param->BRc1o1taunutau=param->BRc1o2udbar=param->BRc1o2csbar=param->BRc1o2enue=param->BRc1o2mnum=param->BRc1o2taunutau=param->BRc1o3udbar=param->BRc1o3csbar=param->BRc1o3enue=param->BRc1o3mnum=param->BRc1o3taunutau=param->BRc1o4udbar=param->BRc1o4csbar=param->BRc1o4enue=param->BRc1o4mnum=param->BRc1o4taunutau=param->BRc1nuele=param->BRc1numlm=param->BRc1elnue=param->BRc1mlnum=param->BRc1tau2nutau=param->BRc1c1Z=param->BRc1c1h=param->BRc1nutaultau=param->BRc1o2W=param->BRc2c1Z=param->BRc2o1W=param->BRc2o2W=param->BRc2c1h=param->BRc2nuele=param->BRc2numlm=param->BRc2nutaultau=param->BRc2elnue=param->BRc2mlnum=param->BRc2tau1nutau=param->BRc2tau2nutau=param->BRo2o1Z=param->BRo2o1h=param->BRo2tau1taubar=param->BRo2tau1bartau=param->BRo2o1gamma=param->BRo2o1ubaru=param->BRo2o1dbard=param->BRo2o1cbarc=param->BRo2o1sbars=param->BRo2o1bbarb=param->BRo2o1tbart=param->BRo2o1ebare=param->BRo2o1mbarm=param->BRo2o1taubartau=param->BRo2o1nuebarnue=param->BRo2o1numbarnum=param->BRo2o1nutaubarnutau=param->BRo2c1ubard=param->BRo2c1dbaru=param->BRo2c1cbars=param->BRo2c1sbarc=param->BRo2c1tbarb=param->BRo2c1bbart=param->BRo2c1nuebare=param->BRo2c1nueebar=param->BRo2c1numbarm=param->BRo2c1nummbar=param->BRo2c1nutaubartau=param->BRo2c1nutautaubar=param->BRo2c2ubard=param->BRo2c2dbaru=param->BRo2c2cbars=param->BRo2c2sbarc=param->BRo2c2tbarb=param->BRo2c2bbart=param->BRo2c2nuebare=param->BRo2c2nueebar=param->BRo2c2numbarm=param->BRo2c2nummbar=param->BRo2c2nutaubartau=param->BRo2c2nutautaubar=param->BRo2elebar=param->BRo2elbare=param->BRo2erebar=param->BRo2erbare=param->BRo2mlmbar=param->BRo2mlbarm=param->BRo2mrmbar=param->BRo2mrbarm=param->BRo2tau2taubar=param->BRo2tau2bartau=param->BRo2nuelnuebar=param->BRo2nuelbarnue=param->BRo2numlnumbar=param->BRo2numlbarnum=param->BRo2nutaulnutaubar=param->BRo2nutaulbarnutau=param->BRo3o1Z=param->BRo3o2Z=param->BRo3c1W=param->BRo3c1barW=param->BRo3o1h=param->BRo3o2h=param->BRo3elebar=param->BRo3elbare=param->BRo3erebar=param->BRo3erbare=param->BRo3mlmbar=param->BRo3mlbarm=param->BRo3mrmbar=param->BRo3mrbarm=param->BRo3tau1taubar=param->BRo3tau1bartau=param->BRo3tau2taubar=param->BRo3tau2bartau=param->BRo3nuelnuebar=param->BRo3nuelbarnue=param->BRo3numlnumbar=param->BRo3numlbarnum=param->BRo3nutaulnutaubar=param->BRo3nutaulbarnutau=param->BRo3o1gamma=param->BRo3o2gamma=param->BRo4o1Z=param->BRo4o2Z=param->BRo4c1W=param->BRo4c1barW=param->BRo4o1h=param->BRo4o2h=param->BRo4elebar=param->BRo4elbare=param->BRo4erebar=param->BRo4erbare=param->BRo4mlmbar=param->BRo4mlbarm=param->BRo4mrmbar=param->BRo4mrbarm=param->BRo4tau1taubar=param->BRo4tau1bartau=param->BRo4tau2taubar=param->BRo4tau2bartau=param->BRo4nuelnuebar=param->BRo4nuelbarnue=param->BRo4numlnumbar=param->BRo4numlbarnum=param->BRo4nutaulnutaubar=param->BRo4nutaulbarnutau=param->BRo4o1gamma=param->BRo4o2gamma=param->BRo4o3gamma=param->BRo5o1Z=param->BRo5o2Z=param->BRo5c1W=param->BRo5c1barW=param->BRo5o1h=param->BRo5o2h=param->BRo5elebar=param->BRo5elbare=param->BRo5erebar=param->BRo5erbare=param->BRo5mlmbar=param->BRo5mlbarm=param->BRo5mrmbar=param->BRo5mrbarm=param->BRo5tau1taubar=param->BRo5tau1bartau=param->BRo5tau2taubar=param->BRo5tau2bartau=param->BRo5nuelnuebar=param->BRo5nuelbarnue=param->BRo5numlnumbar=param->BRo5numlbarnum=param->BRo5nutaulnutaubar=param->BRo5nutaulbarnutau=param->BRo5o1gamma=param->BRo5o2gamma=param->BRo5o3gamma=param->BRh0bb_SM=param->BRh0tautau_SM=param->BRh0WW_SM=param->BRh0gg_SM=param->BRh0gaga_SM=param->BRh0ZZ_SM=param->BRh0bb=param->BRh0tautau=param->BRh0WW=param->BRh0gg=param->BRh0gaga=param->BRh0ZZ=param->BRH0bb_SM=param->BRH0tautau_SM=param->BRH0WW_SM=param->BRH0gg_SM=param->BRH0gaga_SM=param->BRH0ZZ_SM=param->BRH0bb=param->BRH0tautau=param->BRH0WW=param->BRH0gg=param->BRH0gaga=param->BRH0ZZ=param->BRA0bb_SM=param->BRA0tautau_SM=param->BRA0WW_SM=param->BRA0gg_SM=param->BRA0gaga_SM=param->BRA0ZZ_SM=param->BRA0bb=param->BRA0tautau=param->BRA0WW=param->BRA0gg=param->BRA0gaga=param->BRA0ZZ=param->BRh0mumu=param->BRh0ss=param->BRh0cc=param->BRh0tt=param->BRh0gaZ=param->BRh0n1n2=param->BRh0n1n3=param->BRh0n1n4=param->BRh0n2n3=param->BRh0n2n4=param->BRh0c1c1=param->BRh0c1c2=param->BRh0n1n1=param->BRh0n2n2=param->BRH0mumu=param->BRH0ss=param->BRH0cc=param->BRH0tt=param->BRH0gaZ=param->BRH0hZ=param->BRH0n1n2=param->BRH0n1n3=param->BRH0n1n4=param->BRH0n2n3=param->BRH0n2n4=param->BRH0c1c1=param->BRH0c1c2=param->BRH0n1n1=param->BRH0n2n2=param->BRH0hh=param->BRA0mumu=param->BRA0ss=param->BRA0cc=param->BRA0tt=param->BRA0gaZ=param->BRA0hZ=param->BRA0n1n2=param->BRA0n1n3=param->BRA0n1n4=param->BRA0n2n3=param->BRA0n2n4=param->BRA0c1c1=param->BRA0c1c2=param->BRA0n1n1=param->BRA0n2n2=param->BRA0hh=param->BRHmunu=param->BRHtaunu=param->BRHub=param->BRHus=param->BRHcs=param->BRHcb=param->BRHtb=param->BRHWh=param->BRHWA=param->BRHc1n1=param->BRHc1n2=param->BRHc1n3=param->BRHc1n4=param->BRHc2n1=param->BRHc2n2=param->BRh0mumu_SM=param->BRh0ss_SM=param->BRh0cc_SM=param->BRh0tt_SM=param->BRh0gaZ_SM=param->BRh0stau1stau1=param->BRh0stau1stau2=param->BRh0stau2stau2=param->BRH0stau1stau1=param->BRH0stau1stau2=param->BRH0stau2stau2=param->BRA0stau1stau1=param->BRA0stau1stau2=param->BRA0stau2stau2=param->BRh0b1b1=param->BRh0b1b2=param->BRh0b2b2=param->BRH0b1b1=param->BRH0b1b2=param->BRH0b2b2=param->BRA0b1b1=param->BRA0b1b2=param->BRA0b2b2=-1.;
	param->BRH03bb=param->BRH03tautau=param->BRH03WW=param->BRH03gg=param->BRH03gaga=param->BRH03ZZ=param->BRA02bb=param->BRA02tautau=param->BRA02WW=param->BRA02gg=param->BRA02gaga=param->BRA02ZZ=param->BRH03mumu=param->BRH03ss=param->BRH03cc=param->BRH03tt=param->BRH03gaZ=param->BRH03hZ=param->BRH03n1n2=param->BRH03n1n3=param->BRH03n1n4=param->BRH03n2n3=param->BRH03n2n4=param->BRH03c1c1=param->BRH03c1c2=param->BRH03n1n1=param->BRH03n2n2=param->BRH03hh=param->BRA02mumu=param->BRA02ss=param->BRA02cc=param->BRA02tt=param->BRA02gaZ=param->BRA02hZ=param->BRA02n1n2=param->BRA02n1n3=param->BRA02n1n4=param->BRA02n2n3=param->BRA02n2n4=param->BRA02c1c1=param->BRA02c1c2=param->BRA02n1n1=param->BRA02n2n2=param->BRA02hh=param->BRH03stau1stau1=param->BRH03stau1stau2=param->BRH03stau2stau2=param->BRA02stau1stau1=param->BRA02stau1stau2=param->BRA02stau2stau2=param->BRH03b1b1=param->BRH03b1b2=param->BRH03b2b2=param->BRA02b1b1=param->BRA02b1b2=param->BRA02b2b2=-1.;
	/* SM */
	param->width_h0SM=param->mass_H0SM=param->width_H0SM=param->mass_A0SM=param->width_A0SM=0.;

	/* 2HDM */
	param->THDM_model=0;
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		param->lambda_u[ie][je]=0.;
		param->lambda_d[ie][je]=0.;
		param->lambda_l[ie][je]=0.;
	}
	
	/* NP contributions to Wilson coefficients */
    for(ie=1;ie<=20;ie++) param->deltaC[ie]=param->deltaCp[ie]=0.;
	for(ie=1;ie<=4;ie++) param->deltaCQ[ie]=param->deltaCQp[ie]=0.;
	
	/* Flavour physics */
	param->f_B=0.1905;
	param->f_Bs=0.2277;
	param->f_Ds=0.2486;
	param->f_D=0.2135;
	param->f_K=0.156;
	param->fK_fpi=1.193;
	param->f_Kstar_par=0.204;
	param->f_Kstar_perp=0.159;
	param->f_phi_par=0.233;
	param->f_phi_perp=0.191;
	param->m_B=5.27926;
	param->m_Bs=5.36677;
	param->m_Bd=5.27958;
	param->m_pi=0.13957;
	param->m_K=0.493677;
	param->m_K0=0.497614; 
	param->m_Kstar=0.89166;
	param->m_Kstar0=0.89581;
	param->m_D0=1.86484;
	param->m_D=1.86961;
	param->m_Dstar=2.01027;
	param->m_Dstar0=2.00697;
	param->m_Ds=1.9683;
	param->m_phi=1.019461;
	param->life_pi=2.6033e-8;
	param->life_K=1.2380e-8;  
	param->life_B=1.638e-12;
	param->life_Bs=1.511e-12;
	param->life_Bd=1.519e-12;
	param->life_D=1.040e-12;
	param->life_Ds=5.e-13;
	param->a1perp=0.04;
	param->a2perp=0.10;
	param->a1par=0.06;
	param->a2par=0.16;
	param->a1K=0.06;
	param->a2K=0.25;
	param->a1phi_perp=0.;
	param->a1phi_par=0.;
	param->a2phi_perp=0.14;
	param->a2phi_par=0.23;
	param->zeta3A=0.032;
	param->zeta3V=0.013;
	param->wA10=-2.1;
	param->deltatp=0.16;
	param->deltatm=-0.16;
	param->deltatp_phi=0.33;
	param->deltatm_phi=0.;
	param->lambda_Bp=0.46;
	param->lambda_Bsp=0.46;
	param->rho1=0.06;
	param->lambda2=0.12;
	param->BR_BXclnu_exp=0.1065;
	param->fullFF=1;
	param->likelihoodBKstarmumu=1;
	
	/* b -> s gamma parameters */
	param->mu_G2_bsg=0.336;
	param->rho_D3_bsg=0.153;
	param->rho_LS3_bsg=-0.145;
	param->mu_c_bsg=2.;

	/* B -> Kstar FF */
	param->T1_BKstar=0.312;/*0.282 LCSR only, 0.312 LCSR+Lattice, 0.35 Beneke et al. 0106067 */
	param->hadrerrBKstar=10.; /* in percent */
	param->a0V_BKstar=0.376313;
	param->a1V_BKstar=-1.16597;
	param->a2V_BKstar=2.42443;
	param->MV_BKstar=5.415;
	param->a0A1_BKstar=0.29725;
	param->a1A1_BKstar=0.392378;
	param->a2A1_BKstar= 1.18916;
	param->MA1_BKstar=5.829;
	param->a0A12_BKstar=0.265375;
	param->a1A12_BKstar=0.533638;
	param->a2A12_BKstar=0.483166;
	param->MA12_BKstar=5.829;
	param->a0A0_BKstar=0.369196;
	param->a1A0_BKstar=-1.36584;
	param->a2A0_BKstar=0.128191;
	param->MA0_BKstar=5.366;
	param->a0T1_BKstar=0.312055;
	param->a1T1_BKstar=-1.00893;
	param->a2T1_BKstar=1.5272;
	param->MT1_BKstar=5.415;
	param->a0T2_BKstar=0.312055;
	param->a1T2_BKstar=0.496846;
	param->a2T2_BKstar=1.61431;
	param->MT2_BKstar=5.829;
	param->a0T23_BKstar=0.667412;
	param->a1T23_BKstar=1.31812;
	param->a2T23_BKstar=3.82334;
	param->MT23_BKstar=5.829;	

	/* Bs -> phi FF */
	param->hadrerrBsphi=10.; /* in percent */
	param->a0V_Bsphi=0.364478;
    param->a1V_Bsphi=-1.22389;
    param->a2V_Bsphi=3.74061;
	param->MV_Bsphi=5.415;
	param->a0A1_Bsphi=0.288007;
    param->a1A1_Bsphi=0.350826;
    param->a2A1_Bsphi=1.69688;
	param->MA1_Bsphi=5.829;
	param->a0A12_Bsphi=0.267053;
    param->a1A12_Bsphi=0.954402;
    param->a2A12_Bsphi=2.15263;		
	param->MA12_Bsphi=5.829;
	param->a0A0_Bsphi=0.421328;
    param->a1A0_Bsphi=-0.976454;
    param->a2A0_Bsphi=3.2714;
	param->MA0_Bsphi=5.366;
	param->a0T1_Bsphi=0.299475;
    param->a1T1_Bsphi=-1.1013;
    param->a2T1_Bsphi=0.58459;
	param->MT1_Bsphi=5.415;
	param->a0T2_Bsphi=0.299475;
    param->a1T2_Bsphi=0.403564;
    param->a2T2_Bsphi=1.03987;
	param->MT2_Bsphi=5.829;
	param->a0T23_Bsphi=0.65233;
    param->a1T23_Bsphi=2.09622;
    param->a2T23_Bsphi=6.73572;
	param->MT23_Bsphi=5.829;
		
	/* B -> K FF */
	param->hadrerrBK=10.; /* in percent */
	param->a00_BK=0.54;
	param->a10_BK=-1.91;
	param->a20_BK=1.83;
	param->a30_BK=-0.02;
	param->a0p_BK=0.43;
	param->a1p_BK=-0.67;
	param->a2p_BK=-1.12;
	param->a0T_BK=0.4;
	param->a1T_BK=-0.53;
	param->a2T_BK=-0.29;
	param->DmBp_BK=0.04578;
	param->DmBT_BK=0.052;
	
	/* hadronic uncertainties */
	param->BtoKstarlow_ALperp_err_noq2=param->BtoKstarlow_ARperp_err_noq2=param->BtoKstarlow_ALpar_err_noq2=param->BtoKstarlow_ARpar_err_noq2=param->BtoKstarlow_AL0_err_noq2=param->BtoKstarlow_AR0_err_noq2=param->BtoKstarlow_At_err_noq2=param->BtoKstarlow_AS_err_noq2=param->BtoKstarlow_ALperp_err_q2=param->BtoKstarlow_ARperp_err_q2=param->BtoKstarlow_ALpar_err_q2=param->BtoKstarlow_ARpar_err_q2=param->BtoKstarlow_AL0_err_q2=param->BtoKstarlow_AR0_err_q2=param->BtoKstarlow_At_err_q2=param->BtoKstarlow_AS_err_q2=param->BtoKstarhigh_ALperp_err=param->BtoKstarhigh_ARperp_err=param->BtoKstarhigh_ALpar_err=param->BtoKstarhigh_ARpar_err=param->BtoKstarhigh_AL0_err=param->BtoKstarhigh_AR0_err=param->BtoKstarhigh_At_err=param->BtoKstarhigh_AS_err=0.;

	param->BtoKlow_FV_err_noq2=param->BtoKlow_FA_err_noq2=param->BtoKlow_FS_err_noq2=param->BtoKlow_FP_err_noq2=param->BtoKlow_FV_err_q2=param->BtoKlow_FA_err_q2=param->BtoKlow_FS_err_q2=param->BtoKlow_FP_err_q2=param->BtoKhigh_FV_err=param->BtoKhigh_FA_err=param->BtoKhigh_FS_err=param->BtoKhigh_FP_err=0.;	

	param->Bstophilow_ALperp_err_noq2=param->Bstophilow_ARperp_err_noq2=param->Bstophilow_ALpar_err_noq2=param->Bstophilow_ARpar_err_noq2=param->Bstophilow_AL0_err_noq2=param->Bstophilow_AR0_err_noq2=param->Bstophilow_At_err_noq2=param->Bstophilow_AS_err_noq2=param->Bstophilow_ALperp_err_q2=param->Bstophilow_ARperp_err_q2=param->Bstophilow_ALpar_err_q2=param->Bstophilow_ARpar_err_q2=param->Bstophilow_AL0_err_q2=param->Bstophilow_AR0_err_q2=param->Bstophilow_At_err_q2=param->Bstophilow_AS_err_q2=param->Bstophihigh_ALperp_err=param->Bstophihigh_ARperp_err=param->Bstophihigh_ALpar_err=param->Bstophihigh_ARpar_err=param->Bstophihigh_AL0_err=param->Bstophihigh_AR0_err=param->Bstophihigh_At_err=param->Bstophihigh_AS_err=0.;
	
	/* power corrections implementations */
	param->BKstar_implementation=1;
	
	param->hplus0=param->hminus0=param->hplus1=param->hminus1=param->hplus2=param->hminus2=param->hzero0=param->hzero1=param->hzero2=0.; /* hadronic parameters */
	
	param->real_alpha_perp0=-0.06e-4; /*err: 0.21e-4 */
	param->real_alpha_perp1=-6.77e-4; /*err: 0.27e-4 */
	param->real_alpha_perp2=18.96e-4; /*err: 0.59e-4 */
	param->real_alpha_par0=-0.35e-4; /*err: 0.62e-4 */
	param->real_alpha_par1=-3.13e-4; /*err: 0.41e-4 */
	param->real_alpha_par2=12.20e-4; /*err: 1.34e-4 */
	param->real_alpha_zero0= 0.05e-4; /*err: 1.52e-4 */
	param->real_alpha_zero1=17.26e-4; /*err: 1.64e-4 */
	param->imag_alpha_perp0=-0.21e-4; /*err: 2.25e-4 */
	param->imag_alpha_perp1=1.17e-4; /*err: 3.58e-4 */
	param->imag_alpha_perp2=-0.08e-4; /*err: 2.24e-4 */
	param->imag_alpha_par0=-0.04e-4; /*err: 3.67e-4 */
	param->imag_alpha_par1=-2.14e-4; /*err: 2.46e-4 */
	param->imag_alpha_par2=6.03e-4; /*err: 2.50e-4 */
	param->imag_alpha_zero0=-0.05e-4; /*err: 4.99e-4  */
	param->imag_alpha_zero1=4.29e-4; /*err: 3.14e-4 */
	
	param->DeltaC9_M1_q2bar=0.72; /*err: +0.57  -0.37 */
	param->r1_M1=0.10; /*err: +0.02  -0.00 */
	param->r2_M1=1.13; /*err: +0.00  -0.01 */
	param->DeltaC9_M2_q2bar=0.76; /*err: +0.70  -0.41 */
	param->r1_M2=0.09; /*err: +0.01  -0.00 */
	param->r2_M2=1.12; /*err: +0.00  -0.01 */
	param->DeltaC9_M3_q2bar=1.11; /*err: +1.14  -0.70 */
	param->r1_M3=0.06; /*err: +0.04  -0.10 */
	param->r2_M3=1.05; /*err: +0.05  -0.04 */
	
	/* B -> D parameters */
	param->Delta_BD=1.;
	param->rho_D2_BD=1.186;
	param->V1_1_BD=1.074;
			
	/* B -> D* parameters */
	param->Delta_BDstar=1.;
	param->rho_Dstar2_BDstar=1.214;
	param->R1_1_BDstar=1.403;
	param->R2_1_BDstar=0.864;
	param->R3_1_BDstar=0.97;
	param->V1_1_BDstar=1.074;
	param->hA1_1_BDstar=0.921;
	
	/* CKM matrix with Wolfenstein parametrisation */
	param->CKM_lambda=0.22506;
	param->CKM_A=0.811;
	param->CKM_rhobar=0.124;
	param->CKM_etabar=0.356;
	
	/* masses and couplings from PDG 2016 */
	param->mass_u = 2.2e-3;
	param->mass_d = 4.7e-3;
	param->mass_s = 0.096;
	param->mass_c = 1.27;
	param->scheme_c_mass = -1;
	param->mass_b = 4.18;
	param->mass_top_pole = 173.34;
	
	param->mass_e = 0.511e-3;
	param->mass_mu= 0.105658;
	param->mass_tau_pole=1.77686;
	param->mass_tau=param->mass_tau_pole;
	
	param->mass_Z=91.1876;
	param->alphas_MZ=0.1181;
	param->mass_W=80.385;

	param->mass_h0=param->mass_h0SM=125.09;

	param->gp=param->gp_Q=3.57458e-1;
	param->g2=param->g2_Q=6.51908e-1;	
	param->inv_alpha_em=1.27916e2;
	param->Gfermi=1.16637000e-5;

	param->width_Z=2.4952;
	param->width_W=2.085;
	param->width_top=2.;

	return;
}

/*--------------------------------------------------------------------*/

int Les_Houches_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[500],dummy2[500];
	double version,Qtemp;
	int ie,je;
	
	if(!test_file(name))
	{
		param->model=-4;
		return 0;
	}
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MODSEL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")))
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 0: fscanf(lecture,"%d",&param->THDM_model); break;
					case 1:	fscanf(lecture,"%d",&param->model); break;
					case 3:	fscanf(lecture,"%d",&param->NMSSM); break;
					case 4:	fscanf(lecture,"%d",&param->RV); break;
					case 5:	fscanf(lecture,"%d",&param->CPV); break;
					case 6:	fscanf(lecture,"%d",&param->FV); break;
					case 12: fscanf(lecture,"%lf",&param->Q); break;
				}	
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);		

	if(param->NMSSM != 0) param->model=param->NMSSM; 
	if(param->RV != 0) param->model=-2;
	if(param->CPV != 0) param->model=-2;
	if(param->THDM_model !=0) param->model=param->THDM_model;
	
	if(param->model<0) return 0;	
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"SPINFO"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: 	fscanf(lecture,"%s",dummy); 
							if(!strncasecmp(dummy,"ISA",3)) param->generator=1; 
							if(!strncasecmp(dummy,"SOFTSUSY",8)) param->generator=3; 
							if(!strncasecmp(dummy,"SPHENO",6)) param->generator=4; 
							if(!strncasecmp(dummy,"SUSPECT",7)) param->generator=5; 
							if(!strncasecmp(dummy,"NMSSMTOOLS",10)) param->generator=6; 
							if(!strncasecmp(dummy,"2HDMC",5)) param->generator=10; 
							break;
					case 2: if(param->generator==1) 
						{
							fscanf(lecture,"%lf",&version); 
							if(version>=7.80) param->generator=2;
						}
						break;
						
					case 3: if(param->generator==6)
						{
							if(EOF!=fscanf(lecture,"%s",dummy)) sprintf(dummy2,"%s",dummy);	
							while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) sprintf(dummy2,"%s%s",dummy2,dummy);
							if(!strcasecmp("# Chargino too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Neutralinos too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Charged Higgs too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strncasecmp("# Excluded by ee",dummy2,16)) {param->NMSSMcoll=1; break;}
							if(!strncasecmp("# Excluded by stop",dummy2,18)) {param->NMSSMcoll=1; break;}
							if(!strncasecmp("# Excluded by sbottom",dummy2,21)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Squark/gluino too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Selectron/smuon too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Stau too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Landau Pole below MGUT",dummy2)) {param->NMSSMtheory=1; break;}
							if(!strcasecmp("# Unphysical global minimum",dummy2)) {param->NMSSMtheory=1; break;}
							if(!strcasecmp("# Higgs soft masses >> Msusy",dummy2)) {param->NMSSMtheory=1; break;}
							if(!strncasecmp("# Excluded by Upsilon",dummy2,21)) {param->NMSSMups1S=1; break;}
							if(!strncasecmp("# Excluded etab(1S)",dummy2,19)) {param->NMSSMetab1S=1; break;}
						}
						break;
						
					case 4: if(param->generator==3)
						{
							if(EOF!=fscanf(lecture,"%s",dummy)) sprintf(dummy2,"%s",dummy);
							while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) sprintf(dummy2,"%s%c",dummy2,dummy[0]);
							if(strcasecmp("Point invalid: stau LSP",dummy2)) {param->model=-1; fclose(lecture); return 0;}
						}
						else
						{
							param->model=-1; fclose(lecture); return 0;
						}
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SMINPUTS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->inv_alpha_em); break;
					case 2: fscanf(lecture,"%lf",&param->Gfermi); break;
					case 3: fscanf(lecture,"%lf",&param->alphas_MZ); break;
					case 4: fscanf(lecture,"%lf",&param->mass_Z); break;
					case 5: fscanf(lecture,"%lf",&param->mass_b); break;
					case 6: fscanf(lecture,"%lf",&param->mass_top_pole); break;
					case 7: fscanf(lecture,"%lf",&param->mass_tau_pole); break;
					case 8: fscanf(lecture,"%lf",&param->mass_nut); break;
					case 11: fscanf(lecture,"%lf",&param->mass_e); break;
					case 12: fscanf(lecture,"%lf",&param->mass_nue); break;
					case 13: fscanf(lecture,"%lf",&param->mass_mu); break;
					case 14: fscanf(lecture,"%lf",&param->mass_num); break;
					case 21: fscanf(lecture,"%lf",&param->mass_d); break;
					case 22: fscanf(lecture,"%lf",&param->mass_u); break;
					case 23: fscanf(lecture,"%lf",&param->mass_s); break;
					case 24: fscanf(lecture,"%lf",&param->mass_c);param->scheme_c_mass=1; break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"VCKMIN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->CKM_lambda); break;
					case 2: fscanf(lecture,"%lf",&param->CKM_A); break;
					case 3: fscanf(lecture,"%lf",&param->CKM_rhobar); break;
					case 4: fscanf(lecture,"%lf",&param->CKM_etabar); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UPMNSIN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->PMNS_theta12); break;
					case 2: fscanf(lecture,"%lf",&param->PMNS_theta23); break;
					case 3: fscanf(lecture,"%lf",&param->PMNS_theta13); break;
					case 4: fscanf(lecture,"%lf",&param->PMNS_delta13); break;
					case 5: fscanf(lecture,"%lf",&param->PMNS_alpha1); break;
					case 6: fscanf(lecture,"%lf",&param->PMNS_alpha2); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MINPAR"))
		{
			switch(param->model)
			{
				case 1:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)){while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));}
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
							case 1: fscanf(lecture,"%lf",&param->m0); break;
							case 2: fscanf(lecture,"%lf",&param->m12); break;
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
							case 4: fscanf(lecture,"%lf",&param->sign_mu); break;
							case 5: fscanf(lecture,"%lf",&param->A0); break;
						}
				}
					break;
				}
	
				case 2:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
							case 1: fscanf(lecture,"%lf",&param->Lambda); break;
							case 2: fscanf(lecture,"%lf",&param->Mmess); break;
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
							case 4: fscanf(lecture,"%lf",&param->sign_mu); break;
							case 5: fscanf(lecture,"%lf",&param->N5); break;
							case 6: fscanf(lecture,"%lf",&param->cgrav); break;
						}
					}
					break;
				}
	
				case 3:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
							case 1: fscanf(lecture,"%lf",&param->m32); break;
							case 2: fscanf(lecture,"%lf",&param->m0); break;
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
							case 4: fscanf(lecture,"%lf",&param->sign_mu); break;
						}
					}
					break;
				}

				default:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
						}
					}
					break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"EXTPAR"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 0: fscanf(lecture,"%lf",&param->Min); break;
					case 1: fscanf(lecture,"%lf",&param->M1_Min); break;
					case 2: fscanf(lecture,"%lf",&param->M2_Min); break;
					case 3: fscanf(lecture,"%lf",&param->M3_Min); break;	
					case 11: fscanf(lecture,"%lf",&param->At_Min); break;
					case 12: fscanf(lecture,"%lf",&param->Ab_Min); break;
					case 13: fscanf(lecture,"%lf",&param->Atau_Min); break;
					case 21: fscanf(lecture,"%lf",&param->M2H1_Min); break;
					case 22: fscanf(lecture,"%lf",&param->M2H2_Min); break;
					case 23: fscanf(lecture,"%lf",&param->mu_Min); break;
					case 24: fscanf(lecture,"%lf",&param->M2A_Min); break;
					case 25: fscanf(lecture,"%lf",&param->tb_Min); break;
					case 26: fscanf(lecture,"%lf",&param->mA_Min); break;
					case 31: fscanf(lecture,"%lf",&param->MeL_Min); break;
					case 32: fscanf(lecture,"%lf",&param->MmuL_Min); break;
					case 33: fscanf(lecture,"%lf",&param->MtauL_Min); break;
					case 34: fscanf(lecture,"%lf",&param->MeR_Min); break;
					case 35: fscanf(lecture,"%lf",&param->MmuR_Min); break;
					case 36: fscanf(lecture,"%lf",&param->MtauR_Min); break;
					case 41: fscanf(lecture,"%lf",&param->MqL1_Min); break;
					case 42: fscanf(lecture,"%lf",&param->MqL2_Min); break;
					case 43: fscanf(lecture,"%lf",&param->MqL3_Min); break;
					case 44: fscanf(lecture,"%lf",&param->MuR_Min); break;
					case 45: fscanf(lecture,"%lf",&param->McR_Min); break;
					case 46: fscanf(lecture,"%lf",&param->MtR_Min); break;
					case 47: fscanf(lecture,"%lf",&param->MdR_Min); break;
					case 48: fscanf(lecture,"%lf",&param->MsR_Min); break;
					case 49: fscanf(lecture,"%lf",&param->MbR_Min); break;
					case 51: fscanf(lecture,"%lf",&param->N51); break;
					case 52: fscanf(lecture,"%lf",&param->N52); break;
					case 53: fscanf(lecture,"%lf",&param->N53); break;
					case 61: fscanf(lecture,"%lf",&param->lambdaNMSSM_Min); break;
					case 62: fscanf(lecture,"%lf",&param->kappaNMSSM_Min); break;
					case 63: fscanf(lecture,"%lf",&param->AlambdaNMSSM_Min); break;
					case 64: fscanf(lecture,"%lf",&param->AkappaNMSSM_Min); break;
					case 65: fscanf(lecture,"%lf",&param->lambdaSNMSSM_Min); break;
					case 66: fscanf(lecture,"%lf",&param->xiFNMSSM_Min); break;
					case 67: fscanf(lecture,"%lf",&param->xiSNMSSM_Min); break;
					case 68: fscanf(lecture,"%lf",&param->mupNMSSM_Min); break;
					case 69: fscanf(lecture,"%lf",&param->mSp2NMSSM_Min); break;
					case 70: fscanf(lecture,"%lf",&param->mS2NMSSM_Min); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->mass_d); if(isnan(param->mass_d)) {param->model=-3; return 0;} break;
					case 2: fscanf(lecture,"%lf",&param->mass_u); if(isnan(param->mass_u)) {param->model=-3; return 0;} break;
					case 3: fscanf(lecture,"%lf",&param->mass_s); if(isnan(param->mass_s)) {param->model=-3; return 0;} break;
					case 4: fscanf(lecture,"%lf",&param->mass_c_pole); if(isnan(param->mass_c_pole)) {param->model=-3; return 0;} break;
					case 5: fscanf(lecture,"%lf",&param->mass_b_pole); if(isnan(param->mass_b_pole)) {param->model=-3; return 0;} break;
					case 6: fscanf(lecture,"%lf",&param->mass_top_pole); if(isnan(param->mass_top_pole)) {param->model=-3; return 0;} break;
					case 11: fscanf(lecture,"%lf",&param->mass_e); if(isnan(param->mass_e)) {param->model=-3; return 0;} break;
					case 12: fscanf(lecture,"%lf",&param->mass_nue); if(isnan(param->mass_nue)) {param->model=-3; return 0;} break;
					case 13: fscanf(lecture,"%lf",&param->mass_mu); if(isnan(param->mass_mu)) {param->model=-3; return 0;} break;
					case 14: fscanf(lecture,"%lf",&param->mass_num); if(isnan(param->mass_num)) {param->model=-3; return 0;} break;
					case 15: fscanf(lecture,"%lf",&param->mass_tau); if(isnan(param->mass_tau)) {param->model=-3; return 0;} break;
					case 16: fscanf(lecture,"%lf",&param->mass_nut); if(isnan(param->mass_nut)) {param->model=-3; return 0;} break;
					case 21: fscanf(lecture,"%lf",&param->mass_gluon); break;
					case 22: fscanf(lecture,"%lf",&param->mass_photon); break;
					case 23: fscanf(lecture,"%lf",&param->mass_Z0); if(isnan(param->mass_Z0)) {param->model=-3; return 0;} break;
					case 24: fscanf(lecture,"%lf",&param->mass_W); if(isnan(param->mass_W)) {param->model=-3; return 0;} break;
					case 25: fscanf(lecture,"%lf",&param->mass_h0); if(isnan(param->mass_h0)) {param->model=-3; return 0;} break;
					case 35: fscanf(lecture,"%lf",&param->mass_H0); if(isnan(param->mass_H0)) {param->model=-3; return 0;} break;
					case 36: fscanf(lecture,"%lf",&param->mass_A0); if(isnan(param->mass_A0)) {param->model=-3; return 0;} break;
					case 37: fscanf(lecture,"%lf",&param->mass_H); if(isnan(param->mass_H)) {param->model=-3; return 0;} break;
					case 39: fscanf(lecture,"%lf",&param->mass_graviton); break;
					case 45: fscanf(lecture,"%lf",&param->mass_H03); if(isnan(param->mass_H03)) {param->model=-3; return 0;} break;
					case 46: fscanf(lecture,"%lf",&param->mass_A02); if(isnan(param->mass_A02)) {param->model=-3; return 0;} break;
					case 1000001: fscanf(lecture,"%lf",&param->mass_dnl); if(isnan(param->mass_dnl)) {param->model=-3; return 0;} break;
					case 1000002: fscanf(lecture,"%lf",&param->mass_upl); if(isnan(param->mass_upl)) {param->model=-3; return 0;} break;
					case 1000003: fscanf(lecture,"%lf",&param->mass_stl); if(isnan(param->mass_stl)) {param->model=-3; return 0;} break;
					case 1000004: fscanf(lecture,"%lf",&param->mass_chl); if(isnan(param->mass_chl)) {param->model=-3; return 0;} break;
					case 1000005: fscanf(lecture,"%lf",&param->mass_b1); if(isnan(param->mass_b1)) {param->model=-3; return 0;} break;
					case 1000006: fscanf(lecture,"%lf",&param->mass_t1); if(isnan(param->mass_t1)) {param->model=-3; return 0;} break;
					case 1000011: fscanf(lecture,"%lf",&param->mass_el); if(isnan(param->mass_el)) {param->model=-3; return 0;} break;
					case 1000012: fscanf(lecture,"%lf",&param->mass_nuel); if(isnan(param->mass_nuel)) {param->model=-3; return 0;} break;
					case 1000013: fscanf(lecture,"%lf",&param->mass_mul); if(isnan(param->mass_mul)) {param->model=-3; return 0;} break;
					case 1000014: fscanf(lecture,"%lf",&param->mass_numl); if(isnan(param->mass_numl)) {param->model=-3; return 0;} break;
					case 1000015: fscanf(lecture,"%lf",&param->mass_tau1); if(isnan(param->mass_tau1)) {param->model=-3; return 0;} break;
					case 1000016: fscanf(lecture,"%lf",&param->mass_nutl); if(isnan(param->mass_nutl)) {param->model=-3; return 0;} break;
					case 1000021: fscanf(lecture,"%lf",&param->mass_gluino); if(isnan(param->mass_gluino)) {param->model=-3; return 0;} break;
					case 1000022: fscanf(lecture,"%lf",&param->mass_neut[1]); if(isnan(param->mass_neut[1])) {param->model=-3; return 0;} break;
					case 1000023: fscanf(lecture,"%lf",&param->mass_neut[2]); if(isnan(param->mass_neut[2])) {param->model=-3; return 0;} break;
					case 1000024: fscanf(lecture,"%lf",&param->mass_cha1); if(isnan(param->mass_cha1)) {param->model=-3; return 0;} break;
					case 1000025: fscanf(lecture,"%lf",&param->mass_neut[3]); if(isnan(param->mass_neut[3])) {param->model=-3; return 0;} break;
					case 1000035: fscanf(lecture,"%lf",&param->mass_neut[4]); if(isnan(param->mass_neut[4])) {param->model=-3; return 0;} break;
					case 1000037: fscanf(lecture,"%lf",&param->mass_cha2); if(isnan(param->mass_cha2)) {param->model=-3; return 0;} break;
					case 1000039: fscanf(lecture,"%lf",&param->mass_gravitino); break;
					case 1000045: fscanf(lecture,"%lf",&param->mass_neut[5]); if(isnan(param->mass_neut[5])) {param->model=-3; return 0;} break;
					case 2000001: fscanf(lecture,"%lf",&param->mass_dnr); if(isnan(param->mass_dnr)) {param->model=-3; return 0;} break;
					case 2000002: fscanf(lecture,"%lf",&param->mass_upr); if(isnan(param->mass_upr)) {param->model=-3; return 0;} break;
					case 2000003: fscanf(lecture,"%lf",&param->mass_str); if(isnan(param->mass_str)) {param->model=-3; return 0;} break;
					case 2000004: fscanf(lecture,"%lf",&param->mass_chr); if(isnan(param->mass_chr)) {param->model=-3; return 0;} break;
					case 2000005: fscanf(lecture,"%lf",&param->mass_b2); if(isnan(param->mass_b2)) {param->model=-3; return 0;} break;
					case 2000006: fscanf(lecture,"%lf",&param->mass_t2); if(isnan(param->mass_t2)) {param->model=-3; return 0;} break;
					case 2000011: fscanf(lecture,"%lf",&param->mass_er); if(isnan(param->mass_er)) {param->model=-3; return 0;} break;
					case 2000012: fscanf(lecture,"%lf",&param->mass_nuer); break;
					case 2000013: fscanf(lecture,"%lf",&param->mass_mur); if(isnan(param->mass_mur)) {param->model=-3; return 0;} break;
					case 2000014: fscanf(lecture,"%lf",&param->mass_numr); break;
					case 2000015: fscanf(lecture,"%lf",&param->mass_tau2); if(isnan(param->mass_tau2)) {param->model=-3; return 0;} break;
					case 2000016: fscanf(lecture,"%lf",&param->mass_nutr); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
		}
		else if(!strcasecmp(dummy,"STOPMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->stop_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SBOTMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sbot_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"STAUMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->stau_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->neut_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMNMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->neut_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));

				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->charg_Umix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"VMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->charg_Vmix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"GAUGE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->GAUGE_Q!=0.)&&(Qtemp>param->GAUGE_Q)) break;
					param->GAUGE_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->gp_Q); break;
					case 2: fscanf(lecture,"%lf",&param->g2_Q); break;
					case 3: fscanf(lecture,"%lf",&param->g3_Q); break;	
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"YU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->YU_Q!=0.)&&(Qtemp>param->YU_Q)) break;
					param->YU_Q=Qtemp;
				}
				else 				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%lf",&param->yut[ie]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"YD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->YD_Q!=0.)&&(Qtemp>param->YD_Q)) break;
					param->YD_Q=Qtemp;
				}
				else
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%lf",&param->yub[ie]);
					}
				}			
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"YE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->YE_Q!=0.)&&(Qtemp>param->YE_Q)) break;
					param->YE_Q=Qtemp;
				}
				else
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%lf",&param->yutau[ie]);
					}
				}			
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"HMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->HMIX_Q!=0.)&&(Qtemp>param->HMIX_Q)) break;
					param->HMIX_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->mu_Q); break;
					case 2: fscanf(lecture,"%lf",&param->tanb_GUT); break;
					case 3: fscanf(lecture,"%lf",&param->Higgs_VEV); break;	
					case 4: fscanf(lecture,"%lf",&param->mA2_Q); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMHMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->H0_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMAMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->A0_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSOFT"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSOFT_Q!=0.)&&(Qtemp>param->MSOFT_Q)) break;
					param->MSOFT_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->M1_Q); break;
					case 2: fscanf(lecture,"%lf",&param->M2_Q); break;
					case 3: fscanf(lecture,"%lf",&param->M3_Q); break;	
					case 21: fscanf(lecture,"%lf",&param->M2H1_Q); break;
					case 22: fscanf(lecture,"%lf",&param->M2H2_Q); break;
					case 31: fscanf(lecture,"%lf",&param->MeL_Q); break;
					case 32: fscanf(lecture,"%lf",&param->MmuL_Q); break;
					case 33: fscanf(lecture,"%lf",&param->MtauL_Q); break;
					case 34: fscanf(lecture,"%lf",&param->MeR_Q); break;
					case 35: fscanf(lecture,"%lf",&param->MmuR_Q); break;
					case 36: fscanf(lecture,"%lf",&param->MtauR_Q); break;
					case 41: fscanf(lecture,"%lf",&param->MqL1_Q); break;
					case 42: fscanf(lecture,"%lf",&param->MqL2_Q); break;
					case 43: fscanf(lecture,"%lf",&param->MqL3_Q); break;
					case 44: fscanf(lecture,"%lf",&param->MuR_Q); break;
					case 45: fscanf(lecture,"%lf",&param->McR_Q); break;
					case 46: fscanf(lecture,"%lf",&param->MtR_Q); break;
					case 47: fscanf(lecture,"%lf",&param->MdR_Q); break;
					case 48: fscanf(lecture,"%lf",&param->MsR_Q); break;
					case 49: fscanf(lecture,"%lf",&param->MbR_Q); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"AU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->AU_Q!=0.)&&(Qtemp>param->AU_Q)) break;
					param->AU_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->A_u); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 2: fscanf(lecture,"%lf",&param->A_c); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->A_t); break;	
						}
						break;
					}
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"AD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->AD_Q!=0.)&&(Qtemp>param->AD_Q)) break;
					param->AD_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->A_d); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 2: fscanf(lecture,"%lf",&param->A_s); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->A_b); break;	
						}
						break;
					}
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"AE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->AE_Q!=0.)&&(Qtemp>param->AE_Q)) break;
					param->AE_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->A_e); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 2: fscanf(lecture,"%lf",&param->A_mu); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->A_tau); break;	
						}
						break;
					}
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}		
		else if(!strcasecmp(dummy,"NMSSMRUN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->NMSSMRUN_Q!=0.)&&(Qtemp>param->NMSSMRUN_Q)) break;
					param->NMSSMRUN_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					default: while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));break;
					case 1: fscanf(lecture,"%lf",&param->lambdaNMSSM); break;
					case 2: fscanf(lecture,"%lf",&param->kappaNMSSM); break;
					case 3: fscanf(lecture,"%lf",&param->AlambdaNMSSM); break;
					case 4: fscanf(lecture,"%lf",&param->AkappaNMSSM); break;
					case 5: fscanf(lecture,"%lf",&param->lambdaSNMSSM); break;
					case 6: fscanf(lecture,"%lf",&param->xiFNMSSM); break;
					case 7: fscanf(lecture,"%lf",&param->xiSNMSSM); break;
					case 8: fscanf(lecture,"%lf",&param->mupNMSSM); break;
					case 9: fscanf(lecture,"%lf",&param->mSp2NMSSM); break;
					case 10: fscanf(lecture,"%lf",&param->mS2NMSSM); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"USQMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sU_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"DSQMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sD_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SELMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sE_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SNUMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sNU_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSQ2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSQ2_Q!=0.)&&(Qtemp>param->MSQ2_Q)) break;
					param->MSQ2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msq2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSL2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSL2_Q!=0.)&&(Qtemp>param->MSL2_Q)) break;
					param->MSL2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msl2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSD2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSD2_Q!=0.)&&(Qtemp>param->MSD2_Q)) break;
					param->MSD2_Q=Qtemp;
				}
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%lf",&param->MSD2_Q);
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msd2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSU2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSU2_Q!=0.)&&(Qtemp>param->MSU2_Q)) break;
					param->MSU2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msu2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSE2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSE2_Q!=0.)&&(Qtemp>param->MSE2_Q)) break;
					param->MSE2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_mse2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"VCKM"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->CKM_Q!=0.)&&(Qtemp>param->CKM_Q)) break;
					param->CKM_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->CKM[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"IMVCKM"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->IMCKM_Q!=0.)&&(Qtemp>param->IMCKM_Q)) break;
					param->IMCKM_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->IMCKM[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UPMNS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->PMNSU_Q!=0.)&&(Qtemp>param->PMNSU_Q)) break;
					param->PMNSU_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->PMNS_U[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"TU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->TU_Q!=0.)&&(Qtemp>param->TU_Q)) break;
					param->TU_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->TU[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"TD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->TD_Q!=0.)&&(Qtemp>param->TD_Q)) break;
					param->TD_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->TD[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"TE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->TE_Q!=0.)&&(Qtemp>param->TE_Q)) break;
					param->TE_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->TE[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UCOUPL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->lambda_u[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"DCOUPL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->lambda_d[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"LCOUPL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->lambda_l[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}		
	}
	fclose(lecture);
	
	if(param->model<0) return 0;	

	if(alphas_running(param->mass_b/2.,param->mass_top_pole, param->mass_b,param)<0.) 
	{	
		param->model=-3;
		return 0;
	}
	
	if(alphas_running(2.*param->mass_top_pole,param->mass_top_pole, param->mass_b,param)<0.)
	{	
		param->model=-3;
		return 0;
	}

	if(param->mass_c_pole>0.&&param->scheme_c_mass<0)
	{
			if(param->mass_c_pole<1.5) param->mass_c=mcmc_from_pole(param->mass_c_pole,1,param);
			else if(param->mass_c_pole<1.75) param->mass_c=mcmc_from_pole(param->mass_c_pole,2,param);
			else param->mass_c=mcmc_from_pole(param->mass_c_pole,3,param);
	}
	
	SM_Decays_Reader(name,param);
	Higgs_Decays_Reader(name,param);
	SUSY_Decays_Reader(name,param);
	
	slha_adjust(param);
	
	return 1;	
}

/*--------------------------------------------------------------------*/

void slha_adjust(struct parameters* param)
{
	double dum;
	double mass[7];
	int ie,je;
	int iemax=0;
	int iemax1=0;
	int iemax2=0;

	if(param->mass_Z0==0.) param->mass_Z0=param->mass_Z;

	if(param->MSOFT_Q==0.) param->MSOFT_Q=max(param->TU_Q,max(param->TD_Q,max(param->TE_Q,max(param->PMNSU_Q,max(param->CKM_Q,max(param->MSE2_Q,max(param->MSU2_Q,max(param->MSD2_Q,max(param->MSL2_Q,max(param->MSQ2_Q,max(param->NMSSMRUN_Q,max(param->YU_Q,max(param->YD_Q,max(param->YE_Q,max(param->HMIX_Q,max(param->GAUGE_Q,max(param->AU_Q,max(param->AD_Q,param->AE_Q))))))))))))))))));
	
	if(param->Q==0.) param->Q=sqrt(param->mass_t1*param->mass_t2);
	
	if(param->MSOFT_Q==0.) param->TU_Q=param->TD_Q=param->TE_Q=param->PMNSU_Q=param->CKM_Q=param->MSE2_Q=param->MSU2_Q=param->MSD2_Q=param->MSL2_Q=param->MSQ2_Q=param->MSOFT_Q=param->NMSSMRUN_Q=param->YU_Q=param->YD_Q=param->YE_Q=param->HMIX_Q=param->GAUGE_Q=param->AU_Q=param->AD_Q=param->AE_Q=param->Q;			
					
	if(param->tan_beta==0.) param->tan_beta=param->tb_Min;
	if(param->tan_beta==0.) param->tan_beta=param->tanb_GUT;

	if(param->mu_Q==0.) param->mu_Q=param->mu_Min;

	if(param->mass_gluino==0.) param->mass_gluino=fabs(param->M3_Q);
	if(param->mass_gluino==0.) param->mass_gluino=fabs(param->M3_Min);

	if(param->mass_el==0.) param->mass_el=fabs(param->MeL_Q);
	if(param->mass_el==0.) param->mass_el=fabs(param->MeL_Min);
	
	if(param->mass_er==0.) param->mass_er=fabs(param->MeR_Q);
	if(param->mass_er==0.) param->mass_er=fabs(param->MeR_Min);
	
	if(param->mass_mul==0.) param->mass_mul=fabs(param->MmuL_Q);
	if(param->mass_mul==0.) param->mass_mul=fabs(param->MmuL_Min);
	
	if(param->mass_mur==0.) param->mass_mur=fabs(param->MmuR_Q);
	if(param->mass_mur==0.) param->mass_mur=fabs(param->MmuR_Min);
	
	if(param->mass_dnl==0.) param->mass_dnl=fabs(param->MqL1_Q);
	if(param->mass_dnl==0.) param->mass_dnl=fabs(param->MqL1_Min);
	
	if(param->mass_upl==0.) param->mass_upl=fabs(param->MqL1_Q);
	if(param->mass_upl==0.) param->mass_upl=fabs(param->MqL1_Min);
	
	if(param->mass_stl==0.) param->mass_stl=fabs(param->MqL2_Q);
	if(param->mass_stl==0.) param->mass_stl=fabs(param->MqL2_Min);
	
	if(param->mass_chl==0.) param->mass_chl=fabs(param->MqL2_Q);
	if(param->mass_chl==0.) param->mass_chl=fabs(param->MqL2_Min);
	
	if(param->mass_dnr==0.) param->mass_dnr=fabs(param->MdR_Q);
	if(param->mass_dnr==0.) param->mass_dnr=fabs(param->MdR_Min);
	
	if(param->mass_upr==0.) param->mass_upr=fabs(param->MuR_Q);
	if(param->mass_upr==0.) param->mass_upr=fabs(param->MuR_Min);
	
	if(param->mass_str==0.) param->mass_str=fabs(param->MsR_Q);
	if(param->mass_str==0.) param->mass_str=fabs(param->MsR_Min);
	
	if(param->mass_chr==0.) param->mass_chr=fabs(param->McR_Q);
	if(param->mass_chr==0.) param->mass_chr=fabs(param->McR_Min);

	if(param->mass_nuel==0.) param->mass_nuel=fabs(param->MeL_Q);
	if(param->mass_nuel==0.) param->mass_nuel=fabs(param->MeL_Min);
	
	if(param->mass_numl==0.) param->mass_numl=fabs(param->MmuL_Q);
	if(param->mass_numl==0.) param->mass_numl=fabs(param->MmuL_Min);
	
	if((param->tan_beta*param->MSOFT_Q)==0.) param->model=-3;
	
	if(param->MeL_Q==0.) param->MeL_Q=param->mass_el;
	if(param->MmuL_Q==0.) param->MmuL_Q=param->mass_mul;
	if(param->MtauL_Q==0.) param->MtauL_Q=param->mass_tau1;
	if(param->MeR_Q==0.) param->MeR_Q=param->mass_er;
	if(param->MmuR_Q==0.) param->MmuR_Q=param->mass_mur;
	if(param->MtauR_Q==0.) param->MtauR_Q=param->mass_tau2;
	if(param->MqL1_Q==0.) param->MqL1_Q=param->mass_dnl;
	if(param->MqL2_Q==0.) param->MqL2_Q=param->mass_stl;
	if(param->MqL3_Q==0.) param->MqL3_Q=param->mass_b1;
	if(param->MuR_Q==0.) param->MuR_Q=param->mass_upr;
	if(param->McR_Q==0.) param->McR_Q=param->mass_chr;
	if(param->MtR_Q==0.) param->MtR_Q=param->mass_t2;
	if(param->MdR_Q==0.) param->MdR_Q=param->mass_dnr;
	if(param->MsR_Q==0.) param->MsR_Q=param->mass_str;
	if(param->MbR_Q==0.) param->MbR_Q=param->mass_b2;
		
	if(param->A_tau==0.) param->A_tau=param->TE[3][3];
	if(param->A_b==0.) param->A_b=param->TD[3][3];
	if(param->A_t==0.) param->A_t=param->TU[3][3];

	
 	if(param->stop_mix[1][1]==0.) 
	{
		mass[1]=param->mass_upl;
		mass[2]=param->mass_chl;
		mass[3]=param->mass_t1;
		mass[4]=param->mass_upr;
		mass[5]=param->mass_chr;
		mass[6]=param->mass_t2;

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sU_mix[ie][1])>dum)
		{
			iemax=ie;
			dum=fabs(param->sU_mix[ie][1]);
		}
		param->mass_upl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sU_mix[ie][2])>dum)
		{
			iemax=ie;
			dum=fabs(param->sU_mix[ie][2]);
		}
		param->mass_chl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sU_mix[ie][4])>dum)
		{
			iemax=ie;
			dum=fabs(param->sU_mix[ie][4]);
		}
		param->mass_upr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sU_mix[ie][5])>dum)
		{
			iemax=ie;
			dum=fabs(param->sU_mix[ie][5]);
		}
		param->mass_chr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sU_mix[ie][3])>dum)
		{
			iemax1=ie;
			dum=fabs(param->sU_mix[ie][3]);
		}
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sU_mix[ie][6])>dum)
		{
			iemax2=ie;
			dum=fabs(param->sU_mix[ie][6]);
		}
		
		if(iemax1<iemax2)
		{
			param->mass_t1=mass[iemax1];
			param->mass_t2=mass[iemax2];		
			param->stop_mix[1][1]=param->sU_mix[iemax1][3];
			param->stop_mix[1][2]=param->sU_mix[iemax1][6];
		}
		else
		{
			param->mass_t1=mass[iemax2];
			param->mass_t2=mass[iemax1];		
			param->stop_mix[1][1]=param->sU_mix[iemax2][3];
			param->stop_mix[1][2]=param->sU_mix[iemax2][6];
		}
	}	
			
	dum=atan(param->stop_mix[1][2]/param->stop_mix[1][1]);
	
	if(param->generator==1) dum=atan(param->stop_mix[2][1]/param->stop_mix[1][1]);
	
	param->stop_mix[1][1]=cos(dum);
	param->stop_mix[2][1]=-sin(dum);
	param->stop_mix[1][2]=sin(dum);
	param->stop_mix[2][2]=cos(dum);

		
	if(param->sbot_mix[1][1]==0.)
	{
		mass[1]=param->mass_dnl;
		mass[2]=param->mass_stl;
		mass[3]=param->mass_b1;
		mass[4]=param->mass_dnr;
		mass[5]=param->mass_str;
		mass[6]=param->mass_b2;

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sD_mix[ie][1])>dum)
		{
			iemax=ie;
			dum=fabs(param->sD_mix[ie][1]);
		}
		param->mass_dnl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sD_mix[ie][2])>dum)
		{
			iemax=ie;
			dum=fabs(param->sD_mix[ie][2]);
		}
		param->mass_stl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sD_mix[ie][4])>dum)
		{
			iemax=ie;
			dum=fabs(param->sD_mix[ie][4]);
		}
		param->mass_dnr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sD_mix[ie][5])>dum)
		{
			iemax=ie;
			dum=fabs(param->sD_mix[ie][5]);
		}
		param->mass_str=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sD_mix[ie][3])>dum)
		{
			iemax1=ie;
			dum=fabs(param->sD_mix[ie][3]);
		}
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sD_mix[ie][6])>dum)
		{
			iemax2=ie;
			dum=fabs(param->sD_mix[ie][6]);
		}
		
		if(iemax1<iemax2)
		{
			param->mass_b1=mass[iemax1];
			param->mass_b2=mass[iemax2];		
			param->sbot_mix[1][1]=param->sD_mix[iemax1][3];
			param->sbot_mix[1][2]=param->sD_mix[iemax1][6];
		}
		else
		{
			param->mass_b1=mass[iemax2];
			param->mass_b2=mass[iemax1];		
			param->sbot_mix[1][1]=param->sD_mix[iemax2][3];
			param->sbot_mix[1][2]=param->sD_mix[iemax2][6];
		}	
	}	
			
	dum=atan(param->sbot_mix[1][2]/param->sbot_mix[1][1]);
	
	if(param->generator==1) dum=atan(param->sbot_mix[2][1]/param->sbot_mix[1][1]);
	
	param->sbot_mix[1][1]=cos(dum);
	param->sbot_mix[2][1]=-sin(dum);
	param->sbot_mix[1][2]=sin(dum);
	param->sbot_mix[2][2]=cos(dum);
	
			
	if(param->stau_mix[1][1]==0.)
	{
		mass[1]=param->mass_el;
		mass[2]=param->mass_mul;
		mass[3]=param->mass_tau1;
		mass[4]=param->mass_er;
		mass[5]=param->mass_mur;
		mass[6]=param->mass_tau2;

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sE_mix[ie][1])>dum)
		{
			iemax=ie;
			dum=fabs(param->sE_mix[ie][1]);
		}
		param->mass_el=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sE_mix[ie][2])>dum)
		{
			iemax=ie;
			dum=fabs(param->sE_mix[ie][2]);
		}
		param->mass_mul=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sE_mix[ie][4])>dum)
		{
			iemax=ie;
			dum=fabs(param->sE_mix[ie][4]);
		}
		param->mass_er=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sE_mix[ie][5])>dum)
		{
			iemax=ie;
			dum=fabs(param->sE_mix[ie][5]);
		}
		param->mass_mur=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sE_mix[ie][3])>dum)
		{
			iemax1=ie;
			dum=fabs(param->sE_mix[ie][3]);
		}
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(param->sE_mix[ie][6])>dum)
		{
			iemax2=ie;
			dum=fabs(param->sE_mix[ie][6]);
		}
				
		if(iemax1<iemax2)
		{
			param->mass_tau1=mass[iemax1];
			param->mass_tau2=mass[iemax2];		
			param->stau_mix[1][1]=param->sE_mix[iemax1][3];
			param->stau_mix[1][2]=param->sE_mix[iemax1][6];
		}
		else
		{
			param->mass_tau1=mass[iemax2];
			param->mass_tau2=mass[iemax1];		
			param->stau_mix[1][1]=param->sE_mix[iemax2][3];
			param->stau_mix[1][2]=param->sE_mix[iemax2][6];
		}
	}
			
	dum=atan(param->stau_mix[1][2]/param->stau_mix[1][1]);
	
	if(param->generator==1) dum=atan(param->stau_mix[2][1]/param->stau_mix[1][1]);
	
	param->stau_mix[1][1]=cos(dum);
	param->stau_mix[2][1]=-sin(dum);
	param->stau_mix[1][2]=sin(dum);
	param->stau_mix[2][2]=cos(dum);
	
	
	if((param->sNU_mix[1][1]!=0.)||(param->sNU_mix[1][2]!=0.)||(param->sNU_mix[1][3]!=0.))
	{
		mass[1]=param->mass_nuel;
		mass[2]=param->mass_numl;
		mass[3]=param->mass_nutl;

		dum=0.;
		for(ie=1;ie<=3;ie++) if(fabs(param->sNU_mix[ie][1])>dum)
		{
			iemax=ie;
			dum=fabs(param->sNU_mix[ie][1]);
		}
		param->mass_nuel=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=3;ie++) if(fabs(param->sNU_mix[ie][2])>dum)
		{
			iemax=ie;
			dum=fabs(param->sNU_mix[ie][2]);
		}
		param->mass_numl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=3;ie++) if(fabs(param->sNU_mix[ie][3])>dum)
		{
			iemax=ie;
			dum=fabs(param->sNU_mix[ie][3]);
		}
		param->mass_nutl=mass[iemax];		
	}	

	
	if(param->neut_mix[1][1]>0.) for(ie=1;ie<=5;ie++) for(je=1;je<=5;je++) param->neut_mix[ie][je]=-param->neut_mix[ie][je];
	
	param->mass_b_pole=mb_pole(param);
 	param->mass_b_1S=mb_1S(param);
	
	param->mass_t=param->mtmt=mt_mt(param);
	
	if(param->CKM_lambda*param->CKM_A*param->CKM_rhobar!=0.)
	{
		double s12,s13,s23,c12,c13,c23;
		double complex expid;
	
		s12=param->CKM_lambda;
		s23=param->CKM_A*param->CKM_lambda*param->CKM_lambda;
		s13=cabs(param->CKM_A*pow(param->CKM_lambda,3.)*(param->CKM_rhobar+I*param->CKM_etabar)*sqrt(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.))/sqrt(1.-param->CKM_lambda*param->CKM_lambda)/(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.)*(param->CKM_rhobar+I*param->CKM_etabar)));
		expid=(param->CKM_A*pow(param->CKM_lambda,3.)*(param->CKM_rhobar+I*param->CKM_etabar)*sqrt(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.))/sqrt(1.-param->CKM_lambda*param->CKM_lambda)/(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.)*(param->CKM_rhobar+I*param->CKM_etabar)))/s13;
	
		c12=sqrt(1.-s12*s12);
		c13=sqrt(1.-s13*s13);
		c23=sqrt(1.-s23*s23);

		param->CKM[1][1]=c12*c13;
		param->CKM[1][2]=s12*c13;
		param->CKM[1][3]=creal(s13/expid);
		param->IMCKM[1][3]=cimag(s13/expid);
		param->CKM[2][1]=creal(-s12*c23-c12*s23*s13*expid);
		param->IMCKM[2][1]=cimag(-s12*c23-c12*s23*s13*expid);
		param->CKM[2][2]=creal(c12*c23-s12*s23*s13*expid);
		param->IMCKM[2][2]=cimag(c12*c23-s12*s23*s13*expid);
		param->CKM[2][3]=s23*c13;
		param->CKM[3][1]=creal(s12*s23-c12*c23*s13*expid);
		param->IMCKM[3][1]=cimag(s12*s23-c12*c23*s13*expid);
		param->CKM[3][2]=creal(-c12*s23-s12*c23*s13*expid);
		param->IMCKM[3][2]=cimag(-c12*s23-s12*c23*s13*expid);
		param->CKM[3][3]=c23*c13;
	}

	param->Vud=param->CKM[1][1]+I*param->IMCKM[1][1];
	param->Vus=param->CKM[1][2]+I*param->IMCKM[1][2];
	param->Vub=param->CKM[1][3]+I*param->IMCKM[1][3];
	param->Vcd=param->CKM[2][1]+I*param->IMCKM[2][1];
	param->Vcs=param->CKM[2][2]+I*param->IMCKM[2][2];
	param->Vcb=param->CKM[2][3]+I*param->IMCKM[2][3];
	param->Vtd=param->CKM[3][1]+I*param->IMCKM[3][1];
	param->Vts=param->CKM[3][2]+I*param->IMCKM[3][2];
	param->Vtb=param->CKM[3][3]+I*param->IMCKM[3][3];

	if(param->FV!=0) param->widthcalc=2;
	
	return;
}

/*--------------------------------------------------------------------*/

int test_slha(char name[])
/* "container" function scanning the SLHA file "name" and checking if it is valid. A negative value indicates a problem with the SLHA file. */
{
	struct parameters param;
		
	Init_param(&param);
	if(Les_Houches_Reader(name,&param))
	{
		if(param.FV!=0) return 2;
		if(param.NMSSM!=0) return 3;
		if(param.THDM_model!=0) return 10;
		return 1;
	}

	return param.model;
}

/*--------------------------------------------------------------------*/

int SM_Decays_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500],id3[500];
	
	if(!test_file(name)) return 0;	

	param->width_Z=2.4952;
	param->width_W=2.085;
	param->width_top=2.;
	param->BRtbW=param->BRtbH=param->BRtt1o1=param->BRtt1o2=param->BRtt1o3=param->BRtt1o4=param->BRtt2o1=param->BRtt2o2=param->BRtt2o3=param->BRtt2o4=0.;
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
			{
				case 6:
				{
					fscanf(lecture,"%lf",&param->width_top);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{	
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");

							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"24"))) param->BRtbW=atof(dummy);
								else if((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"37"))) param->BRtbH=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"1000022"))) param->BRtt1o1=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"1000023"))) param->BRtt1o2=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"1000025"))) param->BRtt1o3=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"1000035"))) param->BRtt1o4=atof(dummy);
								else if((!strcasecmp(id1,"2000006"))&&(!strcasecmp(id2,"1000022"))) param->BRtt2o1=atof(dummy);
								else if((!strcasecmp(id1,"2000006"))&&(!strcasecmp(id2,"1000023"))) param->BRtt2o2=atof(dummy);
								else if((!strcasecmp(id1,"2000006"))&&(!strcasecmp(id2,"1000025"))) param->BRtt2o3=atof(dummy);
								else if((!strcasecmp(id1,"2000006"))&&(!strcasecmp(id2,"1000035"))) param->BRtt2o4=atof(dummy);
							}
						}
					}	
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}				
				case 23: fscanf(lecture,"%lf",&param->width_Z); break;
				case 24: fscanf(lecture,"%lf",&param->width_W); break;
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);
	
	return 1;	
}

/*--------------------------------------------------------------------*/

int Higgs_Decays_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500];
	
	if(!test_file(name)) return 0;	

	param->width_h0=0.;
	param->width_H0=0.;
	param->width_A0=0.;
	param->width_H=0.;
	param->width_H03=0.;
	param->width_A02=0.;
	
	param->BRh0bb_SM=param->BRh0tautau_SM=param->BRh0WW_SM=param->BRh0gg_SM=param->BRh0gaga_SM=param->BRh0ZZ_SM=param->BRh0bb=param->BRh0tautau=param->BRh0WW=param->BRh0gg=param->BRh0gaga=param->BRh0ZZ=param->BRH0bb_SM=param->BRH0tautau_SM=param->BRH0WW_SM=param->BRH0gg_SM=param->BRH0gaga_SM=param->BRH0ZZ_SM=param->BRH0bb=param->BRH0tautau=param->BRH0WW=param->BRH0gg=param->BRH0gaga=param->BRH0ZZ=param->BRA0bb_SM=param->BRA0tautau_SM=param->BRA0WW_SM=param->BRA0gg_SM=param->BRA0gaga_SM=param->BRA0ZZ_SM=param->BRA0bb=param->BRA0tautau=param->BRA0WW=param->BRA0gg=param->BRA0gaga=param->BRA0ZZ=param->BRh0mumu=param->BRh0ss=param->BRh0cc=param->BRh0tt=param->BRh0gaZ=param->BRh0n1n2=param->BRh0n1n3=param->BRh0n1n4=param->BRh0n2n3=param->BRh0n2n4=param->BRh0c1c1=param->BRh0c1c2=param->BRh0n1n1=param->BRh0n2n2=param->BRH0mumu=param->BRH0ss=param->BRH0cc=param->BRH0tt=param->BRH0gaZ=param->BRH0hZ=param->BRH0n1n2=param->BRH0n1n3=param->BRH0n1n4=param->BRH0n2n3=param->BRH0n2n4=param->BRH0c1c1=param->BRH0c1c2=param->BRH0n1n1=param->BRH0n2n2=param->BRH0hh=param->BRA0mumu=param->BRA0ss=param->BRA0cc=param->BRA0tt=param->BRA0gaZ=param->BRA0hZ=param->BRA0n1n2=param->BRA0n1n3=param->BRA0n1n4=param->BRA0n2n3=param->BRA0n2n4=param->BRA0c1c1=param->BRA0c1c2=param->BRA0n1n1=param->BRA0n2n2=param->BRA0hh=param->BRHmunu=param->BRHtaunu=param->BRHub=param->BRHus=param->BRHcs=param->BRHcb=param->BRHtb=param->BRHWh=param->BRHWA=param->BRHc1n1=param->BRHc1n2=param->BRHc1n3=param->BRHc1n4=param->BRHc2n1=param->BRHc2n2=param->BRh0mumu_SM=param->BRh0ss_SM=param->BRh0cc_SM=param->BRh0tt_SM=param->BRh0gaZ_SM=param->BRh0stau1stau1=param->BRh0stau1stau2=param->BRh0stau2stau2=param->BRH0stau1stau1=param->BRH0stau1stau2=param->BRH0stau2stau2=param->BRA0stau1stau1=param->BRA0stau1stau2=param->BRA0stau2stau2=param->BRh0b1b1=param->BRh0b1b2=param->BRh0b2b2=param->BRH0b1b1=param->BRH0b1b2=param->BRH0b2b2=param->BRA0b1b1=param->BRA0b1b2=param->BRA0b2b2=param->BRH03bb=param->BRH03tautau=param->BRH03WW=param->BRH03gg=param->BRH03gaga=param->BRH03ZZ=param->BRA02bb=param->BRA02tautau=param->BRA02WW=param->BRA02gg=param->BRA02gaga=param->BRA02ZZ=param->BRH03mumu=param->BRH03ss=param->BRH03cc=param->BRH03tt=param->BRH03gaZ=param->BRH03hZ=param->BRH03n1n2=param->BRH03n1n3=param->BRH03n1n4=param->BRH03n2n3=param->BRH03n2n4=param->BRH03c1c1=param->BRH03c1c2=param->BRH03n1n1=param->BRH03n2n2=param->BRH03hh=param->BRA02mumu=param->BRA02ss=param->BRA02cc=param->BRA02tt=param->BRA02gaZ=param->BRA02hZ=param->BRA02n1n2=param->BRA02n1n3=param->BRA02n1n4=param->BRA02n2n3=param->BRA02n2n4=param->BRA02c1c1=param->BRA02c1c2=param->BRA02n1n1=param->BRA02n2n2=param->BRA02hh=param->BRH03stau1stau1=param->BRH03stau1stau2=param->BRH03stau2stau2=param->BRA02stau1stau1=param->BRA02stau1stau2=param->BRA02stau2stau2=param->BRH03b1b1=param->BRH03b1b2=param->BRH03b2b2=param->BRA02b1b1=param->BRA02b1b2=param->BRA02b2b2=0.;
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
			{
				case 25:
				{
					fscanf(lecture,"%lf",&param->width_h0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{	
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRh0bb=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRh0tautau=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRh0mumu=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRh0ss=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRh0cc=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRh0tt=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRh0WW=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRh0gaZ=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000022")))) param->BRh0n1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000022")))) param->BRh0n1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000022")))) param->BRh0n1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000023")))) param->BRh0n2n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000023")))) param->BRh0n2n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000024")))||((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000024")))) param->BRh0c1c1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000037")))||((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1000024")))) param->BRh0c1c2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000037")))||((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-1000024")))) param->BRh0c1c2=atof(dummy);
							else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000022"))) param->BRh0n1n1=atof(dummy);
							else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000023"))) param->BRh0n2n2=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRh0gaga=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRh0gg=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRh0ZZ=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"1000015")))||((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRh0stau1stau1=atof(dummy);
							else if(((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-2000015")))||((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"1000015")))) param->BRh0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRh0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-2000015")))) param->BRh0stau2stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"1000005")))||((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRh0b1b1=atof(dummy);
							else if(((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-2000005")))||((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"1000005")))) param->BRh0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRh0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-2000005")))) param->BRh0b2b2=atof(dummy);
						}
					}	
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 26:
				{
					fscanf(lecture,"%lf",&param->width_h0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRh0bb=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRh0tautau=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRh0mumu=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRh0ss=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRh0cc=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRh0tt=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRh0WW=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRh0gaZ=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000022")))) param->BRh0n1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000022")))) param->BRh0n1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000022")))) param->BRh0n1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000023")))) param->BRh0n2n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000023")))) param->BRh0n2n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000024")))||((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000024")))) param->BRh0c1c1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000037")))||((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1000024")))) param->BRh0c1c2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000037")))||((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-1000024")))) param->BRh0c1c2=atof(dummy);
							else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000022"))) param->BRh0n1n1=atof(dummy);
							else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000023"))) param->BRh0n2n2=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRh0gaga=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRh0gg=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRh0ZZ=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"1000015")))||((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRh0stau1stau1=atof(dummy);
							else if(((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-2000015")))||((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"1000015")))) param->BRh0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRh0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-2000015")))) param->BRh0stau2stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"1000005")))||((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRh0b1b1=atof(dummy);
							else if(((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-2000005")))||((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"1000005")))) param->BRh0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRh0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-2000005")))) param->BRh0b2b2=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 35: 
				{
					fscanf(lecture,"%lf",&param->width_H0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRH0bb=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRH0tautau=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRH0mumu=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRH0ss=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRH0cc=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRH0tt=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRH0WW=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRH0gaZ=atof(dummy);
							else if(((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))) param->BRH0hZ=atof(dummy);
							else if(((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))) param->BRH0hZ=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000022")))) param->BRH0n1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000022")))) param->BRH0n1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000022")))) param->BRH0n1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000023")))) param->BRH0n2n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000023")))) param->BRH0n2n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000024")))||((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000024")))) param->BRH0c1c1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000037")))||((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1000024")))) param->BRH0c1c2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000037")))||((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-1000024")))) param->BRH0c1c2=atof(dummy);
							else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000022"))) param->BRH0n1n1=atof(dummy);
							else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000023"))) param->BRH0n2n2=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRH0gaga=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRH0gg=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRH0ZZ=atof(dummy);
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRH0hh=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRh0ZZ=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"1000015")))||((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRH0stau1stau1=atof(dummy);
							else if(((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-2000015")))||((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"1000015")))) param->BRH0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRH0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-2000015")))) param->BRH0stau2stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"1000005")))||((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRH0b1b1=atof(dummy);
							else if(((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-2000005")))||((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"1000005")))) param->BRH0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRH0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-2000005")))) param->BRH0b2b2=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 36:
				{
					fscanf(lecture,"%lf",&param->width_A0);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRA0bb=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRA0tautau=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRA0mumu=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRA0ss=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRA0cc=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRA0tt=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRA0WW=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRA0gaZ=atof(dummy);
							else if(((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))) param->BRA0hZ=atof(dummy);
							else if(((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))) param->BRA0hZ=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000022")))) param->BRA0n1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000022")))) param->BRA0n1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000022")))) param->BRA0n1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000023")))) param->BRA0n2n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000023")))) param->BRA0n2n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000024")))||((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000024")))) param->BRA0c1c1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000037")))||((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1000024")))) param->BRA0c1c2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000037")))||((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-1000024")))) param->BRA0c1c2=atof(dummy);
							else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000022"))) param->BRA0n1n1=atof(dummy);
							else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000023"))) param->BRA0n2n2=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRA0gaga=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRA0gg=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRA0ZZ=atof(dummy);
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRA0hh=atof(dummy);
							else if((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"26"))) param->BRA0hh=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"1000015")))||((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRA0stau1stau1=atof(dummy);
							else if(((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-2000015")))||((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"1000015")))) param->BRA0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRA0stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-2000015")))) param->BRA0stau2stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"1000005")))||((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRA0b1b1=atof(dummy);
							else if(((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-2000005")))||((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"1000005")))) param->BRA0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRA0b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-2000005")))) param->BRA0b2b2=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 37: 
				{
					fscanf(lecture,"%lf",&param->width_H);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"16"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"16")))) param->BRHtaunu=atof(dummy);
							else if(((!strcasecmp(id1,"14"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"14")))) param->BRHmunu=atof(dummy);
							else if(((!strcasecmp(id1,"2"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"2")))) param->BRHub=atof(dummy);
							else if(((!strcasecmp(id1,"2"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"2")))) param->BRHus=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"4")))) param->BRHcs=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"4")))) param->BRHcb=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"6")))) param->BRHtb=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"25")))||((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"24")))) param->BRHWh=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"26")))||((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"24")))) param->BRHWh=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"36")))||((!strcasecmp(id1,"36"))&&(!strcasecmp(id2,"24")))) param->BRHWA=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"1000022")))||((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000024")))) param->BRHc1n1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000024")))) param->BRHc1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000024")))) param->BRHc1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000024")))) param->BRHc1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"1000022")))||((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000037")))) param->BRHc2n1=atof(dummy);
							else if(((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000037")))) param->BRHc2n2=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 45: 
				{
					fscanf(lecture,"%lf",&param->width_H03);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRH03bb=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRH03tautau=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRH03mumu=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRH03ss=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRH03cc=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRH03tt=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRH03WW=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRH03gaZ=atof(dummy);
							else if(((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))) param->BRH03hZ=atof(dummy);
							else if(((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))) param->BRH03hZ=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000022")))) param->BRH03n1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000022")))) param->BRH03n1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000022")))) param->BRH03n1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000023")))) param->BRH03n2n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000023")))) param->BRH03n2n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000024")))||((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000024")))) param->BRH03c1c1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000037")))||((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1000024")))) param->BRH03c1c2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000037")))||((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-1000024")))) param->BRH03c1c2=atof(dummy);
							else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000022"))) param->BRH03n1n1=atof(dummy);
							else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000023"))) param->BRH03n2n2=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRH03gaga=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRH03gg=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRH03ZZ=atof(dummy);
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRH03hh=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"1000015")))||((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRH03stau1stau1=atof(dummy);
							else if(((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-2000015")))||((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"1000015")))) param->BRH03stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRH03stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-2000015")))) param->BRH03stau2stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"1000005")))||((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRH03b1b1=atof(dummy);
							else if(((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-2000005")))||((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"1000005")))) param->BRH03b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRH03b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-2000005")))) param->BRH03b2b2=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 46:
				{
					fscanf(lecture,"%lf",&param->width_A02);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRA02bb=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRA02tautau=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRA02mumu=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRA02ss=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRA02cc=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRA02tt=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRA02WW=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRA02gaZ=atof(dummy);
							else if(((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"25")))) param->BRA02hZ=atof(dummy);
							else if(((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"26")))) param->BRA02hZ=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000023")))||((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000022")))) param->BRA02n1n2=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000022")))) param->BRA02n1n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000022")))) param->BRA02n1n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000025")))||((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1000023")))) param->BRA02n2n3=atof(dummy);
							else if(((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000035")))||((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1000023")))) param->BRA02n2n4=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000024")))||((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000024")))) param->BRA02c1c1=atof(dummy);
							else if(((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-1000037")))||((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1000024")))) param->BRA02c1c2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1000037")))||((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-1000024")))) param->BRA02c1c2=atof(dummy);
							else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1000022"))) param->BRA02n1n1=atof(dummy);
							else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1000023"))) param->BRA02n2n2=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRA02gaga=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRA02gg=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRA02ZZ=atof(dummy);
							else if((!strcasecmp(id1,"25"))&&(!strcasecmp(id2,"25"))) param->BRA02hh=atof(dummy);
							else if((!strcasecmp(id1,"26"))&&(!strcasecmp(id2,"26"))) param->BRA02hh=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"1000015")))||((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRA02stau1stau1=atof(dummy);
							else if(((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-2000015")))||((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"1000015")))) param->BRA02stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-1000015")))) param->BRA02stau1stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"2000015")))||((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-2000015")))) param->BRA02stau2stau2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"1000005")))||((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRA02b1b1=atof(dummy);
							else if(((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-2000005")))||((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"1000005")))) param->BRA02b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-1000005")))) param->BRA02b1b2=atof(dummy);
							else if(((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"2000005")))||((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-2000005")))) param->BRA02b2b2=atof(dummy);
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);
	
	return 1;	
}

/*--------------------------------------------------------------------*/

int SUSY_Decays_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500],id3[500];
	
	if(!test_file(name)) return 0;	

	param->width_gluino=0.1;
	param->width_t1=0.1;
	param->width_t2=0.1;
	param->width_b1=0.1;
	param->width_b2=0.1;
	param->width_ul=0.1;
	param->width_ur=0.1;
	param->width_dl=0.1;
	param->width_dr=0.1;
	param->width_cl=0.1;
	param->width_cr=0.1;
	param->width_sl=0.1;
	param->width_sr=0.1;
	param->width_el=0.1;
	param->width_er=0.1;
	param->width_ml=0.1;
	param->width_mr=0.1;
	param->width_tau1=0.1;
	param->width_tau2=0.1;
	param->width_nuel=0.1;
	param->width_numl=0.1;
	param->width_nutaul=0.1;
	param->width_c1=0.1;
	param->width_c2=0.1;
	param->width_o1=0.1;
	param->width_o2=0.1;
	param->width_o3=0.1;
	param->width_o4=0.1;
	param->width_o5=0.1;	
	
	param->BRgluinot1tbar=param->BRgluinot1bart=param->BRgluinodldbar=param->BRgluinodlbard=param->BRgluinodrdbar=param->BRgluinodrbard=param->BRgluinoulubar=param->BRgluinoulbaru=param->BRgluinourubar=param->BRgluinourbaru=param->BRgluinoslsbar=param->BRgluinoslbars=param->BRgluinosrsbar=param->BRgluinosrbars=param->BRgluinoclcbar=param->BRgluinoclbarc=param->BRgluinocrcbar=param->BRgluinocrbarc=param->BRgluinob1bbar=param->BRgluinob1barb=param->BRgluinob2bbar=param->BRgluinob2barb=param->BRgluinot2tbar=param->BRgluinot2bart=param->BRgluinoo1g=param->BRgluinoo2g=param->BRgluinoo3g=param->BRgluinoo4g=param->BRgluinoo1ddbar=param->BRgluinoo1uubar=param->BRgluinoo1ssbar=param->BRgluinoo1ccbar=param->BRgluinoo1bbbar=param->BRgluinoo1ttbar=param->BRgluinoo2ddbar=param->BRgluinoo2uubar=param->BRgluinoo2ssbar=param->BRgluinoo2ccbar=param->BRgluinoo2bbbar=param->BRgluinoo2ttbar=param->BRgluinoo3ddbar=param->BRgluinoo3uubar=param->BRgluinoo3ssbar=param->BRgluinoo3ccbar=param->BRgluinoo3bbbar=param->BRgluinoo3ttbar=param->BRgluinoo4ddbar=param->BRgluinoo4uubar=param->BRgluinoo4ssbar=param->BRgluinoo4ccbar=param->BRgluinoo4bbbar=param->BRgluinoo4ttbar=param->BRgluinoc1dubar=param->BRgluinoc1udbar=param->BRgluinoc1scbar=param->BRgluinoc1csbar=param->BRgluinoc1btbar=param->BRgluinoc1tbbar=param->BRgluinoc2dubar=param->BRgluinoc2udbar=param->BRgluinoc2scbar=param->BRgluinoc2csbar=param->BRgluinoc2btbar=param->BRgluinoc2tbbar=param->BRgluinot1barW=param->BRgluinot1W=param->BRgluinot1barH=param->BRgluinot1H=param->BRt1o1t=param->BRt1o2t=param->BRt1o3t=param->BRt1o4t=param->BRt1c1b=param->BRt1c2b=param->BRt1o1c=param->BRt1o1u=param->BRt1gluinoc=param->BRt2o1t=param->BRt2o2t=param->BRt2o3t=param->BRt2o4t=param->BRt2c1b=param->BRt2c2b=param->BRt2t1h=param->BRt2t1Z=param->BRt2b1W=param->BRt2o1c=param->BRt2o1u=param->BRt2gluinoc=param->BRb1o1b=param->BRb1o2b=param->BRb1o3b=param->BRb1o4b=param->BRb1c1t=param->BRb1c2t=param->BRb1gluinob=param->BRb1t1W=param->BRb2o1b=param->BRb2o2b=param->BRb2o3b=param->BRb2o4b=param->BRb2c1t=param->BRb2c2t=param->BRb2gluinob=param->BRb2b1h=param->BRb2b1Z=param->BRb2t1W=param->BRb2t2W=param->BRulo1u=param->BRulo2u=param->BRulo3u=param->BRulo4u=param->BRulc1d=param->BRulc2d=param->BRulgluinou=param->BRuro1u=param->BRuro2u=param->BRuro3u=param->BRuro4u=param->BRurc1d=param->BRurgluinou=param->BRurc2d=param->BRdlo1d=param->BRdlo2d=param->BRdlo3d=param->BRdlo4d=param->BRdlc1u=param->BRdlc2u=param->BRdlgluinod=param->BRdro1d=param->BRdro2d=param->BRdro3d=param->BRdro4d=param->BRdrgluinod=param->BRdrc1u=param->BRdrc2u=param->BRclo1c=param->BRclo2c=param->BRclo3c=param->BRclo4c=param->BRclc1s=param->BRclc2s=param->BRclgluinoc=param->BRcro1c=param->BRcro2c=param->BRcro3c=param->BRcro4c=param->BRcrc1s=param->BRcrgluinoc=param->BRcrc2s=param->BRslo1s=param->BRslo2s=param->BRslo3s=param->BRslo4s=param->BRslc1c=param->BRslc2c=param->BRslgluinos=param->BRsro1s=param->BRsro2s=param->BRsro3s=param->BRsro4s=param->BRsrgluinos=param->BRsrc1c=param->BRsrc2c=param->BRelo1e=param->BRelo2e=param->BRelo3e=param->BRelo4e=param->BRelc1nue=param->BRelc2nue=param->BRero1e=param->BRero2e=param->BRero3e=param->BRero4e=param->BRerc1nue=param->BRerc2nue=param->BRmlo1m=param->BRmlo2m=param->BRmlo3m=param->BRmlo4m=param->BRmlc1num=param->BRmlc2num=param->BRmro1m=param->BRmro2m=param->BRmro3m=param->BRmro4m=param->BRmrc1num=param->BRmrc2num=param->BRtau1o1tau=param->BRtau1o2tau=param->BRtau1o3tau=param->BRtau1o4tau=param->BRtau1c1nutau=param->BRtau1c2nutau=param->BRtau1nutaulH=param->BRtau1nutaulW=param->BRtau2o1tau=param->BRtau2o2tau=param->BRtau2o3tau=param->BRtau2o4tau=param->BRtau2c1nutau=param->BRtau2c2nutau=param->BRtau2tau1h=param->BRtau2tau1Z=param->BRtau2nutaulH=param->BRtau2nutaulW=param->BRtau2tau1H=param->BRtau2tau1A=param->BRnuelo1nue=param->BRnuelo2nue=param->BRnuelo3nue=param->BRnuelo4nue=param->BRnuelc1e=param->BRnuelc2e=param->BRnumlo1num=param->BRnumlo2num=param->BRnumlo3num=param->BRnumlo4num=param->BRnumlc1m=param->BRnumlc2m=param->BRnutaulo1nutau=param->BRnutaulo2nutau=param->BRnutaulo3nutau=param->BRnutaulo4nutau=param->BRnutaulc1tau=param->BRnutaulc2tau=param->BRnutaultau1W=param->BRnutaultau1H=param->BRnutaultau2H=param->BRnutaultau2W=param->BRc1o1W=param->BRc1tau1nutau=param->BRc1o1udbar=param->BRc1o1csbar=param->BRc1o1enue=param->BRc1o1mnum=param->BRc1o1taunutau=param->BRc1o2udbar=param->BRc1o2csbar=param->BRc1o2enue=param->BRc1o2mnum=param->BRc1o2taunutau=param->BRc1o3udbar=param->BRc1o3csbar=param->BRc1o3enue=param->BRc1o3mnum=param->BRc1o3taunutau=param->BRc1o4udbar=param->BRc1o4csbar=param->BRc1o4enue=param->BRc1o4mnum=param->BRc1o4taunutau=param->BRc1nuele=param->BRc1numlm=param->BRc1elnue=param->BRc1mlnum=param->BRc1tau2nutau=param->BRc1c1Z=param->BRc1c1h=param->BRc1nutaultau=param->BRc1o2W=param->BRc2c1Z=param->BRc2o1W=param->BRc2o2W=param->BRc2c1h=param->BRc2nuele=param->BRc2numlm=param->BRc2nutaultau=param->BRc2elnue=param->BRc2mlnum=param->BRc2tau1nutau=param->BRc2tau2nutau=param->BRo2o1Z=param->BRo2o1h=param->BRo2tau1taubar=param->BRo2tau1bartau=param->BRo2o1gamma=param->BRo2o1ubaru=param->BRo2o1dbard=param->BRo2o1cbarc=param->BRo2o1sbars=param->BRo2o1bbarb=param->BRo2o1tbart=param->BRo2o1ebare=param->BRo2o1mbarm=param->BRo2o1taubartau=param->BRo2o1nuebarnue=param->BRo2o1numbarnum=param->BRo2o1nutaubarnutau=param->BRo2c1ubard=param->BRo2c1dbaru=param->BRo2c1cbars=param->BRo2c1sbarc=param->BRo2c1tbarb=param->BRo2c1bbart=param->BRo2c1nuebare=param->BRo2c1nueebar=param->BRo2c1numbarm=param->BRo2c1nummbar=param->BRo2c1nutaubartau=param->BRo2c1nutautaubar=param->BRo2c2ubard=param->BRo2c2dbaru=param->BRo2c2cbars=param->BRo2c2sbarc=param->BRo2c2tbarb=param->BRo2c2bbart=param->BRo2c2nuebare=param->BRo2c2nueebar=param->BRo2c2numbarm=param->BRo2c2nummbar=param->BRo2c2nutaubartau=param->BRo2c2nutautaubar=param->BRo2elebar=param->BRo2elbare=param->BRo2erebar=param->BRo2erbare=param->BRo2mlmbar=param->BRo2mlbarm=param->BRo2mrmbar=param->BRo2mrbarm=param->BRo2tau2taubar=param->BRo2tau2bartau=param->BRo2nuelnuebar=param->BRo2nuelbarnue=param->BRo2numlnumbar=param->BRo2numlbarnum=param->BRo2nutaulnutaubar=param->BRo2nutaulbarnutau=param->BRo3o1Z=param->BRo3o2Z=param->BRo3c1W=param->BRo3c1barW=param->BRo3o1h=param->BRo3o2h=param->BRo3elebar=param->BRo3elbare=param->BRo3erebar=param->BRo3erbare=param->BRo3mlmbar=param->BRo3mlbarm=param->BRo3mrmbar=param->BRo3mrbarm=param->BRo3tau1taubar=param->BRo3tau1bartau=param->BRo3tau2taubar=param->BRo3tau2bartau=param->BRo3nuelnuebar=param->BRo3nuelbarnue=param->BRo3numlnumbar=param->BRo3numlbarnum=param->BRo3nutaulnutaubar=param->BRo3nutaulbarnutau=param->BRo3o1gamma=param->BRo3o2gamma=param->BRo4o1Z=param->BRo4o2Z=param->BRo4c1W=param->BRo4c1barW=param->BRo4o1h=param->BRo4o2h=param->BRo4elebar=param->BRo4elbare=param->BRo4erebar=param->BRo4erbare=param->BRo4mlmbar=param->BRo4mlbarm=param->BRo4mrmbar=param->BRo4mrbarm=param->BRo4tau1taubar=param->BRo4tau1bartau=param->BRo4tau2taubar=param->BRo4tau2bartau=param->BRo4nuelnuebar=param->BRo4nuelbarnue=param->BRo4numlnumbar=param->BRo4numlbarnum=param->BRo4nutaulnutaubar=param->BRo4nutaulbarnutau=param->BRo4o1gamma=param->BRo4o2gamma=param->BRo4o3gamma=param->BRo5o1Z=param->BRo5o2Z=param->BRo5c1W=param->BRo5c1barW=param->BRo5o1h=param->BRo5o2h=param->BRo5elebar=param->BRo5elbare=param->BRo5erebar=param->BRo5erbare=param->BRo5mlmbar=param->BRo5mlbarm=param->BRo5mrmbar=param->BRo5mrbarm=param->BRo5tau1taubar=param->BRo5tau1bartau=param->BRo5tau2taubar=param->BRo5tau2bartau=param->BRo5nuelnuebar=param->BRo5nuelbarnue=param->BRo5numlnumbar=param->BRo5numlbarnum=param->BRo5nutaulnutaubar=param->BRo5nutaulbarnutau=param->BRo5o1gamma=param->BRo5o2gamma=param->BRo5o3gamma=0.;
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
			{
				case 1000021:
				{
					fscanf(lecture,"%lf",&param->width_gluino);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000001"))&&(!strcasecmp(id2,"-1"))) param->BRgluinodldbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000001"))&&(!strcasecmp(id2,"1"))) param->BRgluinodlbard=atof(dummy);
								else if((!strcasecmp(id1,"2000001"))&&(!strcasecmp(id2,"-1"))) param->BRgluinodrdbar=atof(dummy);
								else if((!strcasecmp(id1,"-2000001"))&&(!strcasecmp(id2,"1"))) param->BRgluinodrbard=atof(dummy);
								else if((!strcasecmp(id1,"1000002"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoulubar=atof(dummy);
								else if((!strcasecmp(id1,"-1000002"))&&(!strcasecmp(id2,"2"))) param->BRgluinoulbaru=atof(dummy);
								else if((!strcasecmp(id1,"2000002"))&&(!strcasecmp(id2,"-2"))) param->BRgluinourubar=atof(dummy);
								else if((!strcasecmp(id1,"-2000002"))&&(!strcasecmp(id2,"2"))) param->BRgluinourbaru=atof(dummy);
								else if((!strcasecmp(id1,"1000003"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoslsbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000003"))&&(!strcasecmp(id2,"3"))) param->BRgluinoslbars=atof(dummy);
								else if((!strcasecmp(id1,"2000003"))&&(!strcasecmp(id2,"-3"))) param->BRgluinosrsbar=atof(dummy);
								else if((!strcasecmp(id1,"-2000003"))&&(!strcasecmp(id2,"3"))) param->BRgluinosrbars=atof(dummy);
								else if((!strcasecmp(id1,"1000004"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoclcbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000004"))&&(!strcasecmp(id2,"4"))) param->BRgluinoclbarc=atof(dummy);
								else if((!strcasecmp(id1,"2000004"))&&(!strcasecmp(id2,"-4"))) param->BRgluinocrcbar=atof(dummy);
								else if((!strcasecmp(id1,"-2000004"))&&(!strcasecmp(id2,"4"))) param->BRgluinocrbarc=atof(dummy);
								else if((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"-5"))) param->BRgluinob1bbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000005"))&&(!strcasecmp(id2,"5"))) param->BRgluinob1barb=atof(dummy);
								else if((!strcasecmp(id1,"2000005"))&&(!strcasecmp(id2,"-5"))) param->BRgluinob2bbar=atof(dummy);
								else if((!strcasecmp(id1,"-2000005"))&&(!strcasecmp(id2,"5"))) param->BRgluinob2barb=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"-6"))) param->BRgluinot1tbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000006"))&&(!strcasecmp(id2,"6"))) param->BRgluinot1bart=atof(dummy);
								else if((!strcasecmp(id1,"2000006"))&&(!strcasecmp(id2,"-6"))) param->BRgluinot2tbar=atof(dummy);
								else if((!strcasecmp(id1,"-2000006"))&&(!strcasecmp(id2,"6"))) param->BRgluinot2bart=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"21"))) param->BRgluinoo1g=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"21"))) param->BRgluinoo2g=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"21"))) param->BRgluinoo3g=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"21"))) param->BRgluinoo4g=atof(dummy);
							}
							else 
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1"))&&(!strcasecmp(id2,"-1"))) param->BRgluinoo1ddbar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoo1uubar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"3"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoo1ssbar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoo1ccbar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"5"))&&(!strcasecmp(id2,"-5"))) param->BRgluinoo1bbbar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"6"))&&(!strcasecmp(id2,"-6"))) param->BRgluinoo1ttbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1"))&&(!strcasecmp(id2,"-1"))) param->BRgluinoo2ddbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoo2uubar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"3"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoo2ssbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoo2ccbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"5"))&&(!strcasecmp(id2,"-5"))) param->BRgluinoo2bbbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"6"))&&(!strcasecmp(id2,"-6"))) param->BRgluinoo2ttbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1"))&&(!strcasecmp(id2,"-1"))) param->BRgluinoo3ddbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoo3uubar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"3"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoo3ssbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoo3ccbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"5"))&&(!strcasecmp(id2,"-5"))) param->BRgluinoo3bbbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"6"))&&(!strcasecmp(id2,"-6"))) param->BRgluinoo3ttbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1"))&&(!strcasecmp(id2,"-1"))) param->BRgluinoo4ddbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoo4uubar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"3"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoo4ssbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoo4ccbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"5"))&&(!strcasecmp(id2,"-5"))) param->BRgluinoo4bbbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"6"))&&(!strcasecmp(id2,"-6"))) param->BRgluinoo4ttbar=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id2,"-1"))) param->BRgluinoc1dubar=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"1"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoc1udbar=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoc1scbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"3"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoc1csbar=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"6"))&&(!strcasecmp(id2,"-5"))) param->BRgluinoc1btbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"5"))&&(!strcasecmp(id2,"-6"))) param->BRgluinoc1tbbar=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id2,"-1"))) param->BRgluinoc2dubar=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"1"))&&(!strcasecmp(id2,"-2"))) param->BRgluinoc2udbar=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id2,"-3"))) param->BRgluinoc2scbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"3"))&&(!strcasecmp(id2,"-4"))) param->BRgluinoc2csbar=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"6"))&&(!strcasecmp(id2,"-5"))) param->BRgluinoc2btbar=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"5"))&&(!strcasecmp(id2,"-6"))) param->BRgluinoc2tbbar=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"-5"))&&(!strcasecmp(id2,"-24"))) param->BRgluinot1barW=atof(dummy);
								else if((!strcasecmp(id1,"-1000006"))&&(!strcasecmp(id2,"24"))&&(!strcasecmp(id2,"5"))) param->BRgluinot1W=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"-5"))&&(!strcasecmp(id2,"-37"))) param->BRgluinot1barH=atof(dummy);
								else if((!strcasecmp(id1,"-1000006"))&&(!strcasecmp(id2,"37"))&&(!strcasecmp(id2,"5"))) param->BRgluinot1H=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000006:
				{
					fscanf(lecture,"%lf",&param->width_t1);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								 if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"6"))) param->BRt1o1t=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"6"))) param->BRt1o2t=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"6"))) param->BRt1o3t=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"6"))) param->BRt1o4t=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"5"))) param->BRt1c1b=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"5"))) param->BRt1c2b=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"4"))) param->BRt1o1c=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"2"))) param->BRt1o1u=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"4"))) param->BRt1gluinoc=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000006:
				{
					fscanf(lecture,"%lf",&param->width_t2);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"6"))) param->BRt2o1t=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"6"))) param->BRt2o2t=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"6"))) param->BRt2o3t=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"6"))) param->BRt2o4t=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"5"))) param->BRt2c1b=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"5"))) param->BRt2c2b=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"25"))) param->BRt2t1h=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"23"))) param->BRt2t1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"24"))) param->BRt2b1W=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"4"))) param->BRt2o1c=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"2"))) param->BRt2o1u=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"4"))) param->BRt2gluinoc=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000005:
				{
					fscanf(lecture,"%lf",&param->width_b1);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"5"))) param->BRb1o1b=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"5"))) param->BRb1o2b=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"5"))) param->BRb1o3b=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"5"))) param->BRb1o4b=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"6"))) param->BRb1c1t=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"6"))) param->BRb1c2t=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"5"))) param->BRb1gluinob=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"-24"))) param->BRb1t1W=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000005:
				{
					fscanf(lecture,"%lf",&param->width_b2);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"5"))) param->BRb2o1b=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"5"))) param->BRb2o2b=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"5"))) param->BRb2o3b=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"5"))) param->BRb2o4b=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"6"))) param->BRb2c1t=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"6"))) param->BRb2c2t=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"5"))) param->BRb2gluinob=atof(dummy);
								else if((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"25"))) param->BRb2b1h=atof(dummy);
								else if((!strcasecmp(id1,"1000005"))&&(!strcasecmp(id2,"23"))) param->BRb2b1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000006"))&&(!strcasecmp(id2,"-24"))) param->BRb2t1W=atof(dummy);
								else if((!strcasecmp(id1,"2000006"))&&(!strcasecmp(id2,"-24"))) param->BRb2t2W=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000002:
				{
					fscanf(lecture,"%lf",&param->width_ul);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"2"))) param->BRulo1u=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"2"))) param->BRulo2u=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"2"))) param->BRulo3u=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"2"))) param->BRulo4u=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"1"))) param->BRulc1d=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"1"))) param->BRulc2d=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"2"))) param->BRulgluinou=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000002:
				{
					fscanf(lecture,"%lf",&param->width_ur);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"2"))) param->BRuro1u=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"2"))) param->BRuro2u=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"2"))) param->BRuro3u=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"2"))) param->BRuro4u=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"2"))) param->BRurgluinou=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"1"))) param->BRurc1d=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"1"))) param->BRurc2d=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"2"))) param->BRurgluinou=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000001:
				{
					fscanf(lecture,"%lf",&param->width_dl);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1"))) param->BRdlo1d=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1"))) param->BRdlo2d=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1"))) param->BRdlo3d=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1"))) param->BRdlo4d=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"2"))) param->BRdlc1u=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"2"))) param->BRdlc2u=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"1"))) param->BRdlgluinod=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000001:
				{
					fscanf(lecture,"%lf",&param->width_dr);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"1"))) param->BRdro1d=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"1"))) param->BRdro2d=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"1"))) param->BRdro3d=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"1"))) param->BRdro4d=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"2"))) param->BRdrc1u=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"2"))) param->BRdrc2u=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"1"))) param->BRdrgluinod=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000004:
				{
					fscanf(lecture,"%lf",&param->width_ul);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"4"))) param->BRclo1c=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"4"))) param->BRclo2c=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"4"))) param->BRclo3c=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"4"))) param->BRclo4c=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"3"))) param->BRclc1s=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"3"))) param->BRclc2s=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"4"))) param->BRclgluinoc=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000004:
				{
					fscanf(lecture,"%lf",&param->width_cr);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"4"))) param->BRcro1c=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"4"))) param->BRcro2c=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"4"))) param->BRcro3c=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"4"))) param->BRcro4c=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"3"))) param->BRcrc1s=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"3"))) param->BRcrc2s=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"4"))) param->BRcrgluinoc=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000003:
				{
					fscanf(lecture,"%lf",&param->width_dl);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"3"))) param->BRslo1s=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"3"))) param->BRslo2s=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"3"))) param->BRslo3s=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"3"))) param->BRslo4s=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"4"))) param->BRslc1c=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"4"))) param->BRslc2c=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"3"))) param->BRslgluinos=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000003:
				{
					fscanf(lecture,"%lf",&param->width_dr);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"3"))) param->BRsro1s=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"3"))) param->BRsro2s=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"3"))) param->BRsro3s=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"3"))) param->BRsro4s=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"4"))) param->BRsrc1c=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"4"))) param->BRsrc2c=atof(dummy);
								else if((!strcasecmp(id1,"1000021"))&&(!strcasecmp(id2,"3"))) param->BRsrgluinos=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000011:
				{
					fscanf(lecture,"%lf",&param->width_el);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"11"))) param->BRelo1e=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"11"))) param->BRelo2e=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"11"))) param->BRelo3e=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"11"))) param->BRelo4e=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"12"))) param->BRelc1nue=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"12"))) param->BRelc2nue=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000011:
				{
					fscanf(lecture,"%lf",&param->width_er);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"11"))) param->BRero1e=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"11"))) param->BRero2e=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"11"))) param->BRero3e=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"11"))) param->BRero4e=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"12"))) param->BRerc1nue=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"12"))) param->BRerc2nue=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000013:
				{
					fscanf(lecture,"%lf",&param->width_ml);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"13"))) param->BRmlo1m=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"13"))) param->BRmlo2m=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"13"))) param->BRmlo3m=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"13"))) param->BRmlo4m=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"14"))) param->BRmlc1num=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"14"))) param->BRmlc2num=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000013:
				{
					fscanf(lecture,"%lf",&param->width_mr);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"13"))) param->BRmro1m=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"13"))) param->BRmro2m=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"13"))) param->BRmro3m=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"13"))) param->BRmro4m=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"14"))) param->BRmrc1num=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"14"))) param->BRmrc2num=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000015:
				{
					fscanf(lecture,"%lf",&param->width_tau1);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"15"))) param->BRtau1o1tau=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"15"))) param->BRtau1o2tau=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"15"))) param->BRtau1o3tau=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"15"))) param->BRtau1o4tau=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"16"))) param->BRtau1c1nutau=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"16"))) param->BRtau1c2nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-37"))) param->BRtau1nutaulH=atof(dummy);
								else if((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-24"))) param->BRtau1nutaulW=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 2000015:
				{
					fscanf(lecture,"%lf",&param->width_tau2);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"15"))) param->BRtau2o1tau=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"15"))) param->BRtau2o2tau=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"15"))) param->BRtau2o3tau=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"15"))) param->BRtau2o4tau=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"16"))) param->BRtau2c1nutau=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"16"))) param->BRtau2c2nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-37"))) param->BRtau2nutaulH=atof(dummy);
								else if((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-24"))) param->BRtau2nutaulW=atof(dummy);
								else if((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"25"))) param->BRtau2tau1h=atof(dummy);
								else if((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"35"))) param->BRtau2tau1H=atof(dummy);
								else if((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"36"))) param->BRtau2tau1A=atof(dummy);
								else if((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"23"))) param->BRtau2tau1Z=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000012:
				{
					fscanf(lecture,"%lf",&param->width_nuel);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"12"))) param->BRnuelo1nue=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"12"))) param->BRnuelo2nue=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"12"))) param->BRnuelo3nue=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"12"))) param->BRnuelo4nue=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"11"))) param->BRnuelc1e=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"11"))) param->BRnuelc2e=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000014:
				{
					fscanf(lecture,"%lf",&param->width_numl);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"14"))) param->BRnumlo1num=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"14"))) param->BRnumlo2num=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"14"))) param->BRnumlo3num=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"14"))) param->BRnumlo4num=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"13"))) param->BRnumlc1m=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"13"))) param->BRnumlc2m=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000016:
				{
					fscanf(lecture,"%lf",&param->width_nutaul);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"16"))) param->BRnutaulo1nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"16"))) param->BRnutaulo2nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"16"))) param->BRnutaulo3nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"16"))) param->BRnutaulo4nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"15"))) param->BRnutaulc1tau=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"15"))) param->BRnutaulc2tau=atof(dummy);
								else if((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"-37"))) param->BRnutaultau1H=atof(dummy);
								else if((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"-37"))) param->BRnutaultau2H=atof(dummy);
								else if((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"-24"))) param->BRnutaultau1W=atof(dummy);
								else if((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"-24"))) param->BRnutaultau2W=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000024:
				{
					fscanf(lecture,"%lf",&param->width_c1);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if(((!strcasecmp(id1,"1000012"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"1000012"))&&(!strcasecmp(id1,"-11")))) param->BRc1nuele=atof(dummy);
								else if(((!strcasecmp(id1,"1000014"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"1000014"))&&(!strcasecmp(id2,"-13")))) param->BRc1numlm=atof(dummy);
								else if(((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id2,"1000016"))&&(!strcasecmp(id2,"-15")))) param->BRc1nutaultau=atof(dummy);
								else if(((!strcasecmp(id1,"-1000011"))&&(!strcasecmp(id2,"12")))||((!strcasecmp(id2,"-1000011"))&&(!strcasecmp(id2,"12")))) param->BRc1elnue=atof(dummy);
								else if(((!strcasecmp(id1,"-1000013"))&&(!strcasecmp(id2,"14")))||((!strcasecmp(id2,"-1000013"))&&(!strcasecmp(id2,"14")))) param->BRc1mlnum=atof(dummy);
								else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-1000015"))&&(!strcasecmp(id2,"16")))) param->BRc1tau1nutau=atof(dummy);
								else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-2000015"))&&(!strcasecmp(id2,"16")))) param->BRc1tau2nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"23"))) param->BRc1c1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"24"))) param->BRc1o1W=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"24"))) param->BRc1o2W=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"25"))) param->BRc1c1h=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"24"))) param->BRc1o1W=atof(dummy);
							}
							else
							{ 
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id1,"-1"))) param->BRc1o1udbar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id1,"-1"))) param->BRc1o1csbar=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-11"))&&(!strcasecmp(id1,"12"))) param->BRc1o1enue=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-13"))&&(!strcasecmp(id1,"14"))) param->BRc1o1mnum=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-15"))&&(!strcasecmp(id1,"16"))) param->BRc1o1taunutau=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id1,"-1"))) param->BRc1o2udbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id1,"-1"))) param->BRc1o2csbar=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"-11"))&&(!strcasecmp(id1,"12"))) param->BRc1o2enue=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"-13"))&&(!strcasecmp(id1,"14"))) param->BRc1o2mnum=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"-15"))&&(!strcasecmp(id1,"16"))) param->BRc1o2taunutau=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id1,"-1"))) param->BRc1o3udbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id1,"-1"))) param->BRc1o3csbar=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"-11"))&&(!strcasecmp(id1,"12"))) param->BRc1o3enue=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"-13"))&&(!strcasecmp(id1,"14"))) param->BRc1o3mnum=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"-15"))&&(!strcasecmp(id1,"16"))) param->BRc1o3taunutau=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"2"))&&(!strcasecmp(id1,"-1"))) param->BRc1o4udbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"4"))&&(!strcasecmp(id1,"-1"))) param->BRc1o4csbar=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"-11"))&&(!strcasecmp(id1,"12"))) param->BRc1o4enue=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"-13"))&&(!strcasecmp(id1,"14"))) param->BRc1o4mnum=atof(dummy);
								else if((!strcasecmp(id1,"1000035"))&&(!strcasecmp(id2,"-15"))&&(!strcasecmp(id1,"16"))) param->BRc1o4taunutau=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000037:
				{
					fscanf(lecture,"%lf",&param->width_c2);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if(((!strcasecmp(id1,"1000012"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"1000012"))&&(!strcasecmp(id1,"-11")))) param->BRc2nuele=atof(dummy);
								else if(((!strcasecmp(id1,"1000014"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"1000014"))&&(!strcasecmp(id2,"-13")))) param->BRc2numlm=atof(dummy);
								else if(((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id2,"1000016"))&&(!strcasecmp(id2,"-15")))) param->BRc2nutaultau=atof(dummy);
								else if(((!strcasecmp(id1,"-1000011"))&&(!strcasecmp(id2,"12")))||((!strcasecmp(id2,"-1000011"))&&(!strcasecmp(id2,"12")))) param->BRc2elnue=atof(dummy);
								else if(((!strcasecmp(id1,"-1000013"))&&(!strcasecmp(id2,"14")))||((!strcasecmp(id2,"-1000013"))&&(!strcasecmp(id2,"14")))) param->BRc2mlnum=atof(dummy);
								else if(((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-1000015"))&&(!strcasecmp(id2,"16")))) param->BRc2tau1nutau=atof(dummy);
								else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-2000015"))&&(!strcasecmp(id2,"16")))) param->BRc2tau2nutau=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"23"))) param->BRc2c1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"24"))) param->BRc2o1W=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"24"))) param->BRc2o2W=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"25"))) param->BRc2c1h=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000022:
				{
					fscanf(lecture,"%lf",&param->width_o1);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000023:
				{
					fscanf(lecture,"%lf",&param->width_o2);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{ 
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"23"))) param->BRo2o1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"25"))) param->BRo2o1h=atof(dummy);
								else if((!strcasecmp(id1,"1000015"))&&(!strcasecmp(id2,"-15"))) param->BRo2tau1taubar=atof(dummy);
								else if((!strcasecmp(id1,"-1000015"))&&(!strcasecmp(id2,"15"))) param->BRo2tau1bartau=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"22"))) param->BRo2o1gamma=atof(dummy);
								else if(((!strcasecmp(id1,"1000011"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"1000011"))&&(!strcasecmp(id1,"-11")))) param->BRo2elebar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000011"))&&(!strcasecmp(id2,"11")))||((!strcasecmp(id2,"-1000011"))&&(!strcasecmp(id1,"11")))) param->BRo2elbare=atof(dummy);
								else if(((!strcasecmp(id1,"2000011"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"2000011"))&&(!strcasecmp(id1,"-11")))) param->BRo2erebar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000011"))&&(!strcasecmp(id2,"11")))||((!strcasecmp(id2,"-2000011"))&&(!strcasecmp(id1,"11")))) param->BRo2erbare=atof(dummy);
								else if(((!strcasecmp(id1,"1000013"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"1000013"))&&(!strcasecmp(id1,"-13")))) param->BRo2mlmbar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000013"))&&(!strcasecmp(id2,"13")))||((!strcasecmp(id2,"-1000013"))&&(!strcasecmp(id1,"13")))) param->BRo2mlbarm=atof(dummy);
								else if(((!strcasecmp(id1,"2000013"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"2000013"))&&(!strcasecmp(id1,"-13")))) param->BRo2mrmbar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000013"))&&(!strcasecmp(id2,"13")))||((!strcasecmp(id2,"-2000013"))&&(!strcasecmp(id1,"13")))) param->BRo2mrbarm=atof(dummy);
								else if(((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id2,"2000015"))&&(!strcasecmp(id1,"-15")))) param->BRo2tau2taubar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"15")))||((!strcasecmp(id2,"-2000015"))&&(!strcasecmp(id1,"15")))) param->BRo2tau2bartau=atof(dummy);
								else if(((!strcasecmp(id1,"1000012"))&&(!strcasecmp(id2,"-12")))||((!strcasecmp(id2,"1000012"))&&(!strcasecmp(id1,"-12")))) param->BRo2nuelnuebar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000012"))&&(!strcasecmp(id2,"12")))||((!strcasecmp(id2,"-1000012"))&&(!strcasecmp(id1,"12")))) param->BRo2nuelbarnue=atof(dummy);
								else if(((!strcasecmp(id1,"1000014"))&&(!strcasecmp(id2,"-14")))||((!strcasecmp(id2,"1000014"))&&(!strcasecmp(id1,"-14")))) param->BRo2numlnumbar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000014"))&&(!strcasecmp(id2,"14")))||((!strcasecmp(id2,"-1000014"))&&(!strcasecmp(id1,"14")))) param->BRo2numlbarnum=atof(dummy);
								else if(((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-16")))||((!strcasecmp(id2,"1000016"))&&(!strcasecmp(id1,"-16")))) param->BRo2nutaulnutaubar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000016"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-1000016"))&&(!strcasecmp(id1,"16")))) param->BRo2nutaulbarnutau=atof(dummy);
							}
							else
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-2"))&&(!strcasecmp(id3,"2"))) param->BRo2o1ubaru=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-1"))&&(!strcasecmp(id3,"1"))) param->BRo2o1dbard=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-4"))&&(!strcasecmp(id3,"4"))) param->BRo2o1cbarc=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-3"))&&(!strcasecmp(id3,"3"))) param->BRo2o1sbars=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-5"))&&(!strcasecmp(id3,"5"))) param->BRo2o1bbarb=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-6"))&&(!strcasecmp(id3,"6"))) param->BRo2o1tbart=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-11"))&&(!strcasecmp(id3,"11"))) param->BRo2o1ebare=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-13"))&&(!strcasecmp(id3,"13"))) param->BRo2o1mbarm=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-15"))&&(!strcasecmp(id3,"15"))) param->BRo2o1taubartau=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-12"))&&(!strcasecmp(id3,"12"))) param->BRo2o1nuebarnue=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-14"))&&(!strcasecmp(id3,"14"))) param->BRo2o1numbarnum=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"-16"))&&(!strcasecmp(id3,"16"))) param->BRo2o1nutaubarnutau=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-2"))&&(!strcasecmp(id3,"1"))) param->BRo2c1ubard=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"-1"))&&(!strcasecmp(id3,"2"))) param->BRo2c1dbaru=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-4"))&&(!strcasecmp(id3,"3"))) param->BRo2c1cbars=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"-3"))&&(!strcasecmp(id3,"4"))) param->BRo2c1sbarc=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-6"))&&(!strcasecmp(id3,"5"))) param->BRo2c1tbarb=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"-5"))&&(!strcasecmp(id3,"6"))) param->BRo2c1bbart=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-12"))&&(!strcasecmp(id3,"11"))) param->BRo2c1nuebare=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"12"))&&(!strcasecmp(id3,"-11"))) param->BRo2c1nueebar=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-14"))&&(!strcasecmp(id3,"13"))) param->BRo2c1numbarm=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"14"))&&(!strcasecmp(id3,"-13"))) param->BRo2c1nummbar=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-16"))&&(!strcasecmp(id3,"15"))) param->BRo2c1nutaubartau=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"16"))&&(!strcasecmp(id3,"-15"))) param->BRo2c1nutautaubar=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-2"))&&(!strcasecmp(id3,"1"))) param->BRo2c2ubard=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"-1"))&&(!strcasecmp(id3,"2"))) param->BRo2c2dbaru=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-4"))&&(!strcasecmp(id3,"3"))) param->BRo2c2cbars=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"-3"))&&(!strcasecmp(id3,"4"))) param->BRo2c2sbarc=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-6"))&&(!strcasecmp(id3,"5"))) param->BRo2c2tbarb=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"-5"))&&(!strcasecmp(id3,"6"))) param->BRo2c2bbart=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-12"))&&(!strcasecmp(id3,"11"))) param->BRo2c2nuebare=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"12"))&&(!strcasecmp(id3,"-11"))) param->BRo2c2nueebar=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-14"))&&(!strcasecmp(id3,"13"))) param->BRo2c2numbarm=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"14"))&&(!strcasecmp(id3,"-13"))) param->BRo2c2nummbar=atof(dummy);
								else if((!strcasecmp(id1,"1000037"))&&(!strcasecmp(id2,"-16"))&&(!strcasecmp(id3,"15"))) param->BRo2c2nutaubartau=atof(dummy);
								else if((!strcasecmp(id1,"-1000037"))&&(!strcasecmp(id2,"16"))&&(!strcasecmp(id3,"-15"))) param->BRo2c2nutautaubar=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000025:
				{
					fscanf(lecture,"%lf",&param->width_o3);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"23"))) param->BRo3o1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"23"))) param->BRo3o2Z=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-24"))) param->BRo3c1W=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"24"))) param->BRo3c1barW=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"22"))) param->BRo3o1gamma=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"22"))) param->BRo3o2gamma=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"25"))) param->BRo3o1h=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"25"))) param->BRo3o2h=atof(dummy);
								else if(((!strcasecmp(id1,"1000011"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"1000011"))&&(!strcasecmp(id1,"-11")))) param->BRo3elebar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000011"))&&(!strcasecmp(id2,"11")))||((!strcasecmp(id2,"-1000011"))&&(!strcasecmp(id1,"11")))) param->BRo3elbare=atof(dummy);
								else if(((!strcasecmp(id1,"2000011"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"2000011"))&&(!strcasecmp(id1,"-11")))) param->BRo3erebar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000011"))&&(!strcasecmp(id2,"11")))||((!strcasecmp(id2,"-2000011"))&&(!strcasecmp(id1,"11")))) param->BRo3erbare=atof(dummy);
								else if(((!strcasecmp(id1,"1000013"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"1000013"))&&(!strcasecmp(id1,"-13")))) param->BRo3mlmbar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000013"))&&(!strcasecmp(id2,"13")))||((!strcasecmp(id2,"-1000013"))&&(!strcasecmp(id1,"13")))) param->BRo3mlbarm=atof(dummy);
								else if(((!strcasecmp(id1,"2000013"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"2000013"))&&(!strcasecmp(id1,"-13")))) param->BRo3mrmbar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000013"))&&(!strcasecmp(id2,"13")))||((!strcasecmp(id2,"-2000013"))&&(!strcasecmp(id1,"13")))) param->BRo3mrbarm=atof(dummy);
								else if(((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id2,"2000015"))&&(!strcasecmp(id1,"-15")))) param->BRo3tau2taubar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"15")))||((!strcasecmp(id2,"-2000015"))&&(!strcasecmp(id1,"15")))) param->BRo3tau2bartau=atof(dummy);
								else if(((!strcasecmp(id1,"1000012"))&&(!strcasecmp(id2,"-12")))||((!strcasecmp(id2,"1000012"))&&(!strcasecmp(id1,"-12")))) param->BRo3nuelnuebar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000012"))&&(!strcasecmp(id2,"12")))||((!strcasecmp(id2,"-1000012"))&&(!strcasecmp(id1,"12")))) param->BRo3nuelbarnue=atof(dummy);
								else if(((!strcasecmp(id1,"1000014"))&&(!strcasecmp(id2,"-14")))||((!strcasecmp(id2,"1000014"))&&(!strcasecmp(id1,"-14")))) param->BRo3numlnumbar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000014"))&&(!strcasecmp(id2,"14")))||((!strcasecmp(id2,"-1000014"))&&(!strcasecmp(id1,"14")))) param->BRo3numlbarnum=atof(dummy);
								else if(((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-16")))||((!strcasecmp(id2,"1000016"))&&(!strcasecmp(id1,"-16")))) param->BRo3nutaulnutaubar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000016"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-1000016"))&&(!strcasecmp(id1,"16")))) param->BRo3nutaulbarnutau=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
				case 1000035:
				{
					fscanf(lecture,"%lf",&param->width_o4);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(((atof(dummy)*(atoi(dummy)-atof(dummy)))!=0.)||((atoi(dummy)==1)&&(strcasecmp(dummy,"1"))))
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(!strcasecmp(nda,"3")) fscanf(lecture,"%s",id3); else sprintf(id3,"0");
							
							if(!strcasecmp(id3,"0"))
							{
								if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"23"))) param->BRo4o1Z=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"23"))) param->BRo4o2Z=atof(dummy);
								else if((!strcasecmp(id1,"1000024"))&&(!strcasecmp(id2,"-24"))) param->BRo4c1W=atof(dummy);
								else if((!strcasecmp(id1,"-1000024"))&&(!strcasecmp(id2,"24"))) param->BRo4c1barW=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"22"))) param->BRo4o1gamma=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"22"))) param->BRo4o2gamma=atof(dummy);
								else if((!strcasecmp(id1,"1000025"))&&(!strcasecmp(id2,"22"))) param->BRo4o3gamma=atof(dummy);
								else if((!strcasecmp(id1,"1000022"))&&(!strcasecmp(id2,"25"))) param->BRo4o1h=atof(dummy);
								else if((!strcasecmp(id1,"1000023"))&&(!strcasecmp(id2,"25"))) param->BRo4o2h=atof(dummy);
								else if(((!strcasecmp(id1,"1000011"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"1000011"))&&(!strcasecmp(id1,"-11")))) param->BRo4elebar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000011"))&&(!strcasecmp(id2,"11")))||((!strcasecmp(id2,"-1000011"))&&(!strcasecmp(id1,"11")))) param->BRo4elbare=atof(dummy);
								else if(((!strcasecmp(id1,"2000011"))&&(!strcasecmp(id2,"-11")))||((!strcasecmp(id2,"2000011"))&&(!strcasecmp(id1,"-11")))) param->BRo4erebar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000011"))&&(!strcasecmp(id2,"11")))||((!strcasecmp(id2,"-2000011"))&&(!strcasecmp(id1,"11")))) param->BRo4erbare=atof(dummy);
								else if(((!strcasecmp(id1,"1000013"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"1000013"))&&(!strcasecmp(id1,"-13")))) param->BRo4mlmbar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000013"))&&(!strcasecmp(id2,"13")))||((!strcasecmp(id2,"-1000013"))&&(!strcasecmp(id1,"13")))) param->BRo4mlbarm=atof(dummy);
								else if(((!strcasecmp(id1,"2000013"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id2,"2000013"))&&(!strcasecmp(id1,"-13")))) param->BRo4mrmbar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000013"))&&(!strcasecmp(id2,"13")))||((!strcasecmp(id2,"-2000013"))&&(!strcasecmp(id1,"13")))) param->BRo4mrbarm=atof(dummy);
								else if(((!strcasecmp(id1,"2000015"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id2,"2000015"))&&(!strcasecmp(id1,"-15")))) param->BRo4tau2taubar=atof(dummy);
								else if(((!strcasecmp(id1,"-2000015"))&&(!strcasecmp(id2,"15")))||((!strcasecmp(id2,"-2000015"))&&(!strcasecmp(id1,"15")))) param->BRo4tau2bartau=atof(dummy);
								else if(((!strcasecmp(id1,"1000012"))&&(!strcasecmp(id2,"-12")))||((!strcasecmp(id2,"1000012"))&&(!strcasecmp(id1,"-12")))) param->BRo4nuelnuebar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000012"))&&(!strcasecmp(id2,"12")))||((!strcasecmp(id2,"-1000012"))&&(!strcasecmp(id1,"12")))) param->BRo4nuelbarnue=atof(dummy);
								else if(((!strcasecmp(id1,"1000014"))&&(!strcasecmp(id2,"-14")))||((!strcasecmp(id2,"1000014"))&&(!strcasecmp(id1,"-14")))) param->BRo4numlnumbar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000014"))&&(!strcasecmp(id2,"14")))||((!strcasecmp(id2,"-1000014"))&&(!strcasecmp(id1,"14")))) param->BRo4numlbarnum=atof(dummy);
								else if(((!strcasecmp(id1,"1000016"))&&(!strcasecmp(id2,"-16")))||((!strcasecmp(id2,"1000016"))&&(!strcasecmp(id1,"-16")))) param->BRo4nutaulnutaubar=atof(dummy);
								else if(((!strcasecmp(id1,"-1000016"))&&(!strcasecmp(id2,"16")))||((!strcasecmp(id2,"-1000016"))&&(!strcasecmp(id1,"16")))) param->BRo4nutaulbarnutau=atof(dummy);
							}
						}
					}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					break;
				}
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);
	
	return 1;	
}
