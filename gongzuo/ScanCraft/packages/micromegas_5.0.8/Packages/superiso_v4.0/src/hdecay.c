#include "include.h"
#include "include_dm.h"
#include "hdecay.h"

int Hdecay(char name[], struct parameters* param)
{
	FILE *tmp,*input,*output;
	char dummy[500],dum[5],tmp_char[500],namedir[300];
	int fail;
	char *curdir;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.hdtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"%s/slha.in",namedir);
	input=fopen(name,"r");
	output=fopen(tmp_char,"w");

	while(EOF != fscanf(input,"%c",dummy))
	{
		if(!strncasecmp("\n",dummy,1)) 
		{
			fprintf(output,"%c",dummy[0]);
			if(EOF != fscanf(input,"%c",dummy))
	
			if(!strncasecmp("d",dummy,1)) 
			{
				dum[0]=dummy[0];
				fscanf(input,"%c",dummy);
				
			if(!strncasecmp("e",dummy,1)) 
			{
				dum[1]=dummy[0];
				fscanf(input,"%c",dummy);
			
			if(!strncasecmp("c",dummy,1)) 
			{
				dum[2]=dummy[0];
				fscanf(input,"%c",dummy);
			
			if(!strncasecmp("a",dummy,1)) 
			{
				dum[3]=dummy[0];
				fscanf(input,"%c",dummy);
							
			if(!strncasecmp("y",dummy,1)) 
			{
				dum[4]=dummy[0];
				
				if(!strncasecmp("decay",dum,5)) 
				{
					while(EOF != fscanf(input,"%c",dummy))
					if(!strncasecmp("\n",dummy,1))
					{
						fscanf(input,"%c",dummy);
		
					if(!strncasecmp("b",dummy,1))
					{
						dum[0]=dummy[0];
						fscanf(input,"%c",dummy);
						
					if(!strncasecmp("l",dummy,1))
					{
						dum[1]=dummy[0];
						fscanf(input,"%c",dummy);
	
					if(!strncasecmp("o",dummy,1))
					{
						dum[2]=dummy[0];
						fscanf(input,"%c",dummy);
						
					if(!strncasecmp("c",dummy,1))
					{
						dum[3]=dummy[0];
						fscanf(input,"%c",dummy);
	
					if(!strncasecmp("k",dummy,1))
					{
						dum[4]=dummy[0];
					
					if(!strncasecmp("block",dum,5))
					{
						fprintf(output,"%s",dum);
						break;
					}
					}
					}
					}
					}
					}
					}
				}
				else fprintf(output,"%s",dum); 
			}
		} else fprintf(output,"%c%c%c%c",dum[0],dum[1],dum[2],dummy[0]);
		} else fprintf(output,"%c%c%c",dum[0],dum[1],dummy[0]);
		} else fprintf(output,"%c%c",dum[0],dummy[0]);
		} else fprintf(output,"%c",dummy[0]);
		} else fprintf(output,"%c",dummy[0]);
	}
	fclose(input);
	fclose(output);
	
	chdir(namedir);
	
	tmp=fopen("hdecay.in","w");
	
	fprintf(tmp,"SLHAIN   = 1\n");
	fprintf(tmp,"SLHAOUT  = 1\n");
	fprintf(tmp,"COUPVAR  = 1\n");
	fprintf(tmp,"HIGGS    = 5\n");
	fprintf(tmp,"OMIT ELW = 1\n");
	fprintf(tmp,"SM4      = 0\n");
	fprintf(tmp,"FERMPHOB = 0\n");
	fprintf(tmp,"2HDM     = 0\n");
	fprintf(tmp,"MODEL    = 1\n");
	fprintf(tmp,"TGBET    = %.5e\n",param->tan_beta);
	fprintf(tmp,"MABEG    = %.5e\n",param->mass_A0);
	fprintf(tmp,"MAEND    = %.5e\n",param->mass_A0);
	fprintf(tmp,"NMA      = 1\n");
	fprintf(tmp,"********************* hMSSM (MODEL = 10) *********************************\n");
	fprintf(tmp,"MHL      = 125.D0\n");
	fprintf(tmp,"**************************************************************************\n");
	fprintf(tmp,"ALS(MZ)  = %.5e\n",param->alphas_MZ);
	fprintf(tmp,"MSBAR(2) = %.5e\n",param->mass_s);
	fprintf(tmp,"MCBAR(3) = %.5e\n",param->mass_c);
	fprintf(tmp,"MBBAR(MB)= %.5e\n",param->mass_b);
	fprintf(tmp,"MT       = %.5e\n",param->mass_top_pole);
	fprintf(tmp,"MTAU     = %.5e\n",param->mass_tau);
	fprintf(tmp,"MMUON    = %.5e\n",param->mass_mu);
	fprintf(tmp,"1/ALPHA  = %.5e\n",137.0359895e0);
	fprintf(tmp,"GF       = %.5e\n",param->Gfermi);
	fprintf(tmp,"GAMW     = %.5e\n",param->width_W);
	fprintf(tmp,"GAMZ     = %.5e\n",param->width_Z);
	fprintf(tmp,"MZ       = %.5e\n",param->mass_Z);
	fprintf(tmp,"MW       = %.5e\n",param->mass_W);
	fprintf(tmp,"VTB      = %.5e\n",cabs(param->Vtb));
	fprintf(tmp,"VTS      = %.5e\n",cabs(param->Vts));
	fprintf(tmp,"VTD      = %.5e\n",cabs(param->Vtd));
	fprintf(tmp,"VCB      = %.5e\n",cabs(param->Vcb));
	fprintf(tmp,"VCS      = %.5e\n",cabs(param->Vcs));
	fprintf(tmp,"VCD      = %.5e\n",cabs(param->Vcd));
	fprintf(tmp,"VUB      = %.5e\n",cabs(param->Vub));
	fprintf(tmp,"VUS      = %.5e\n",cabs(param->Vus));
	fprintf(tmp,"VUD      = %.5e\n",cabs(param->Vud));
	fprintf(tmp,"********************* 4TH GENERATION *************************************\n");
	fprintf(tmp,"  SCENARIO FOR ELW. CORRECTIONS TO H -> GG (EVERYTHING IN GEV):\n");
	fprintf(tmp,"  GG_ELW = 1: MTP = 500    MBP = 450    MNUP = 375    MEP = 450\n");
	fprintf(tmp,"  GG_ELW = 2: MBP = MNUP = MEP = 600    MTP = MBP+50*(1+LOG(M_H/115)/5)\n");
	fprintf(tmp,"\n");
	fprintf(tmp,"GG_ELW   = 1\n");
	fprintf(tmp,"MTP      = 500.D0\n");
	fprintf(tmp,"MBP      = 450.D0\n");
	fprintf(tmp,"MNUP     = 375.D0\n");
	fprintf(tmp,"MEP      = 450.D0\n");
	fprintf(tmp,"************************** 2 Higgs Doublet Model *************************\n");
	fprintf(tmp,"  TYPE: 1 (I), 2 (II), 3 (lepton-specific), 4 (flipped)\n");
	fprintf(tmp,"  PARAM: 1 (masses), 2 (lambda_i)\n");
	fprintf(tmp,"\n");
	fprintf(tmp,"PARAM    = 2\n");
	fprintf(tmp,"TYPE     = 1\n");
	fprintf(tmp,"********************\n");
	fprintf(tmp,"TGBET2HDM= 1.0D0\n");
	fprintf(tmp,"M_12^2   = 25600.D0\n");
	fprintf(tmp,"******************** PARAM=1:\n");
	fprintf(tmp,"ALPHA_H  = -0.14D0\n");
	fprintf(tmp,"MHL      = 125.D0\n");
	fprintf(tmp,"MHH      = 210.D0\n");
	fprintf(tmp,"MHA      = 130.D0\n");
	fprintf(tmp,"MH+-     = 130.D0\n");
	fprintf(tmp,"******************** PARAM=2:\n");
	fprintf(tmp,"LAMBDA1  = 2.6885665050462264D0\n");
	fprintf(tmp,"LAMBDA2  = 0.000156876030254505681D0\n");
	fprintf(tmp,"LAMBDA3  = 0.46295674052962260D0\n");
	fprintf(tmp,"LAMBDA4  = 0.96605498373771792D0\n");
	fprintf(tmp,"LAMBDA5  = -0.88138084173680198D0\n");
	fprintf(tmp,"**************************************************************************\n");
	fprintf(tmp,"SUSYSCALE= %.5e\n",param->MSOFT_Q);
	fprintf(tmp,"MU       = %.5e\n",param->mu_Q);
	fprintf(tmp,"M2       = %.5e\n",param->M2_Q);
	fprintf(tmp,"MGLUINO  = %.5e\n",param->mass_gluino);
	fprintf(tmp,"MSL1     = %.5e\n",param->MeL_Q);
	fprintf(tmp,"MER1     = %.5e\n",param->MeR_Q);
	fprintf(tmp,"MQL1     = %.5e\n",param->MqL1_Q);
	fprintf(tmp,"MUR1     = %.5e\n",param->MuR_Q);
	fprintf(tmp,"MDR1     = %.5e\n",param->MdR_Q);
	fprintf(tmp,"MSL      = %.5e\n",param->MtauL_Q);
	fprintf(tmp,"MER      = %.5e\n",param->MtauR_Q);
	fprintf(tmp,"MSQ      = %.5e\n",param->MqL3_Q);
	fprintf(tmp,"MUR      = %.5e\n",param->MtR_Q);
	fprintf(tmp,"MDR      = %.5e\n",param->MbR_Q);
	fprintf(tmp,"AL       = %.5e\n",param->A_tau);
	fprintf(tmp,"AU       = %.5e\n",param->A_t);
	fprintf(tmp,"AD       = %.5e\n",param->A_b);
	fprintf(tmp,"ON-SHELL = 0\n");
	fprintf(tmp,"ON-SH-WZ = 0\n");
	fprintf(tmp,"IPOLE    = 0\n");
	fprintf(tmp,"OFF-SUSY = 0\n");
	fprintf(tmp,"INDIDEC  = 0\n");
	fprintf(tmp,"NF-GG    = 5\n");
	fprintf(tmp,"IGOLD    = 0\n");
	fprintf(tmp,"MPLANCK  = %.5e\n",2.4e18);
	fprintf(tmp,"MGOLD    = %.5e\n",1.e-13);
	fprintf(tmp,"******************* VARIATION OF HIGGS COUPLINGS *************************\n");
	fprintf(tmp,"ELWK     = 1\n");
	fprintf(tmp,"CW       = 1.D0\n");
	fprintf(tmp,"CZ       = 1.D0\n");
	fprintf(tmp,"Ctau     = 1.D0\n");
	fprintf(tmp,"Cmu      = 1.D0\n");
	fprintf(tmp,"Ct       = 1.D0\n");
	fprintf(tmp,"Cb       = 1.D0\n");
	fprintf(tmp,"Cc       = 1.D0\n");
	fprintf(tmp,"Cs       = 1.D0\n");
	fprintf(tmp,"Cgaga    = 0.D0\n");
	fprintf(tmp,"Cgg      = 0.D0\n");
	fprintf(tmp,"CZga     = 0.D0\n");
	fprintf(tmp,"********************* 4TH GENERATION *************************************\n");
	fprintf(tmp,"Ctp      = 0.D0\n");
	fprintf(tmp,"Cbp      = 0.D0\n");
	fprintf(tmp,"Cnup     = 0.D0\n");
	fprintf(tmp,"Cep      = 0.D0\n");
	fclose(tmp);
	
	sprintf(tmp_char,"%s > tmp",HDECAY);

	fail=system(tmp_char);

	if((fail)||(!test_file("slha.out")))	
	{	
		chdir(curdir);
		sprintf(tmp_char,"rm -rf %s",namedir);
 		system(tmp_char);
		printf("Problem with Hdecay...\n");
		return 2;
	}
	
	SM_Decays_Reader("slha.out",param);
	Higgs_Decays_Reader("slha.out",param);

	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.) 
	return 0; else return 2;
}

/* -------------------------------- */

int HdecaySM(char name[], double mh, struct parameters* param)
{
	FILE *tmp;
	char tmp_char[500],namedir[300];
	int fail;
	char *curdir;
	
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.hdtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);
	
	tmp=fopen("hdecay.in","w");
	
	fprintf(tmp,"SLHAIN   = 0\n");
	fprintf(tmp,"SLHAOUT  = 1\n");
	fprintf(tmp,"COUPVAR  = 1\n");
	fprintf(tmp,"HIGGS    = 0\n");
	fprintf(tmp,"OMIT ELW = 1\n");
	fprintf(tmp,"SM4      = 0\n");
	fprintf(tmp,"FERMPHOB = 0\n");
	fprintf(tmp,"2HDM     = 0\n");
	fprintf(tmp,"MODEL    = 1\n");
	fprintf(tmp,"TGBET    = %.5e\n",param->tan_beta);
	fprintf(tmp,"MABEG    = %.5e\n",param->mass_h0);
	fprintf(tmp,"MAEND    = %.5e\n",param->mass_h0);
	fprintf(tmp,"NMA      = 1\n");
	fprintf(tmp,"********************* hMSSM (MODEL = 10) *********************************\n");
	fprintf(tmp,"MHL      = 125.D0\n");
	fprintf(tmp,"**************************************************************************\n");
	fprintf(tmp,"ALS(MZ)  = %.5e\n",param->alphas_MZ);
	fprintf(tmp,"MSBAR(2) = %.5e\n",param->mass_s);
	fprintf(tmp,"MCBAR(3) = %.5e\n",param->mass_c);
	fprintf(tmp,"MBBAR(MB)= %.5e\n",param->mass_b);
	fprintf(tmp,"MT       = %.5e\n",param->mass_top_pole);
	fprintf(tmp,"MTAU     = %.5e\n",param->mass_tau);
	fprintf(tmp,"MMUON    = %.5e\n",param->mass_mu);
	fprintf(tmp,"1/ALPHA  = %.5e\n",137.0359895e0);
	fprintf(tmp,"GF       = %.5e\n",param->Gfermi);
	fprintf(tmp,"GAMW     = %.5e\n",param->width_W);
	fprintf(tmp,"GAMZ     = %.5e\n",param->width_Z);
	fprintf(tmp,"MZ       = %.5e\n",param->mass_Z);
	fprintf(tmp,"MW       = %.5e\n",param->mass_W);
	fprintf(tmp,"VTB      = %.5e\n",cabs(param->Vtb));
	fprintf(tmp,"VTS      = %.5e\n",cabs(param->Vts));
	fprintf(tmp,"VTD      = %.5e\n",cabs(param->Vtd));
	fprintf(tmp,"VCB      = %.5e\n",cabs(param->Vcb));
	fprintf(tmp,"VCS      = %.5e\n",cabs(param->Vcs));
	fprintf(tmp,"VCD      = %.5e\n",cabs(param->Vcd));
	fprintf(tmp,"VUB      = %.5e\n",cabs(param->Vub));
	fprintf(tmp,"VUS      = %.5e\n",cabs(param->Vus));
	fprintf(tmp,"VUD      = %.5e\n",cabs(param->Vud));
	fprintf(tmp,"********************* 4TH GENERATION *************************************\n");
	fprintf(tmp,"  SCENARIO FOR ELW. CORRECTIONS TO H -> GG (EVERYTHING IN GEV):\n");
	fprintf(tmp,"  GG_ELW = 1: MTP = 500    MBP = 450    MNUP = 375    MEP = 450\n");
	fprintf(tmp,"  GG_ELW = 2: MBP = MNUP = MEP = 600    MTP = MBP+50*(1+LOG(M_H/115)/5)\n");
	fprintf(tmp,"\n");
	fprintf(tmp,"GG_ELW   = 1\n");
	fprintf(tmp,"MTP      = 500.D0\n");
	fprintf(tmp,"MBP      = 450.D0\n");
	fprintf(tmp,"MNUP     = 375.D0\n");
	fprintf(tmp,"MEP      = 450.D0\n");
	fprintf(tmp,"************************** 2 Higgs Doublet Model *************************\n");
	fprintf(tmp,"  TYPE: 1 (I), 2 (II), 3 (lepton-specific), 4 (flipped)\n");
	fprintf(tmp,"  PARAM: 1 (masses), 2 (lambda_i)\n");
	fprintf(tmp,"\n");
	fprintf(tmp,"PARAM    = 2\n");
	fprintf(tmp,"TYPE     = 1\n");
	fprintf(tmp,"********************\n");
	fprintf(tmp,"TGBET2HDM= 1.0D0\n");
	fprintf(tmp,"M_12^2   = 25600.D0\n");
	fprintf(tmp,"******************** PARAM=1:\n");
	fprintf(tmp,"ALPHA_H  = -0.14D0\n");
	fprintf(tmp,"MHL      = 125.D0\n");
	fprintf(tmp,"MHH      = 210.D0\n");
	fprintf(tmp,"MHA      = 130.D0\n");
	fprintf(tmp,"MH+-     = 130.D0\n");
	fprintf(tmp,"******************** PARAM=2:\n");
	fprintf(tmp,"LAMBDA1  = 2.6885665050462264D0\n");
	fprintf(tmp,"LAMBDA2  = 0.000156876030254505681D0\n");
	fprintf(tmp,"LAMBDA3  = 0.46295674052962260D0\n");
	fprintf(tmp,"LAMBDA4  = 0.96605498373771792D0\n");
	fprintf(tmp,"LAMBDA5  = -0.88138084173680198D0\n");
	fprintf(tmp,"**************************************************************************\n");
	fprintf(tmp,"SUSYSCALE= %.5e\n",param->MSOFT_Q);
	fprintf(tmp,"MU       = %.5e\n",param->mu_Q);
	fprintf(tmp,"M2       = %.5e\n",param->M2_Q);
	fprintf(tmp,"MGLUINO  = %.5e\n",param->mass_gluino);
	fprintf(tmp,"MSL1     = %.5e\n",param->MeL_Q);
	fprintf(tmp,"MER1     = %.5e\n",param->MeR_Q);
	fprintf(tmp,"MQL1     = %.5e\n",param->MqL1_Q);
	fprintf(tmp,"MUR1     = %.5e\n",param->MuR_Q);
	fprintf(tmp,"MDR1     = %.5e\n",param->MdR_Q);
	fprintf(tmp,"MSL      = %.5e\n",param->MtauL_Q);
	fprintf(tmp,"MER      = %.5e\n",param->MtauR_Q);
	fprintf(tmp,"MSQ      = %.5e\n",param->MqL3_Q);
	fprintf(tmp,"MUR      = %.5e\n",param->MtR_Q);
	fprintf(tmp,"MDR      = %.5e\n",param->MbR_Q);
	fprintf(tmp,"AL       = %.5e\n",param->A_tau);
	fprintf(tmp,"AU       = %.5e\n",param->A_t);
	fprintf(tmp,"AD       = %.5e\n",param->A_b);
	fprintf(tmp,"ON-SHELL = 0\n");
	fprintf(tmp,"ON-SH-WZ = 0\n");
	fprintf(tmp,"IPOLE    = 0\n");
	fprintf(tmp,"OFF-SUSY = 0\n");
	fprintf(tmp,"INDIDEC  = 0\n");
	fprintf(tmp,"NF-GG    = 5\n");
	fprintf(tmp,"IGOLD    = 0\n");
	fprintf(tmp,"MPLANCK  = %.5e\n",2.4e18);
	fprintf(tmp,"MGOLD    = %.5e\n",1.e-13);
	fprintf(tmp,"******************* VARIATION OF HIGGS COUPLINGS *************************\n");
	fprintf(tmp,"ELWK     = 1\n");
	fprintf(tmp,"CW       = 1.D0\n");
	fprintf(tmp,"CZ       = 1.D0\n");
	fprintf(tmp,"Ctau     = 1.D0\n");
	fprintf(tmp,"Cmu      = 1.D0\n");
	fprintf(tmp,"Ct       = 1.D0\n");
	fprintf(tmp,"Cb       = 1.D0\n");
	fprintf(tmp,"Cc       = 1.D0\n");
	fprintf(tmp,"Cs       = 1.D0\n");
	fprintf(tmp,"Cgaga    = 0.D0\n");
	fprintf(tmp,"Cgg      = 0.D0\n");
	fprintf(tmp,"CZga     = 0.D0\n");
	fprintf(tmp,"********************* 4TH GENERATION *************************************\n");
	fprintf(tmp,"Ctp      = 0.D0\n");
	fprintf(tmp,"Cbp      = 0.D0\n");
	fprintf(tmp,"Cnup     = 0.D0\n");
	fprintf(tmp,"Cep      = 0.D0\n");
	
	fclose(tmp);
	
	sprintf(tmp_char,"%s",HDECAY);

	fail=system(tmp_char);

	if((fail)||(!test_file("slha.out")))	
	{	
		chdir(curdir);
		sprintf(tmp_char,"rm -rf %s",namedir);
 		system(tmp_char);
		return 0;
	}
	
	FILE *lecture;
	char dummy[500],nda[500],id1[500],id2[500];
	
	lecture = fopen("slha.out","r");
	
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 36: fscanf(lecture,"%lf",&param->mass_h0SM); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);	
		}
		else if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(lecture,"%s",dummy); 
			switch(abs(atoi(dummy)))
			{
				case 25: 
				{
					fscanf(lecture,"%lf",&param->width_h0SM);
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(atof(dummy)>0.) 
						{
							fscanf(lecture,"%s",nda);
							fscanf(lecture,"%s",id1);
							fscanf(lecture,"%s",id2);
							if(((!strcasecmp(id1,"5"))&&(!strcasecmp(id2,"-5")))||((!strcasecmp(id1,"-5"))&&(!strcasecmp(id2,"5")))) param->BRh0bb_SM=atof(dummy);
							else if(((!strcasecmp(id1,"15"))&&(!strcasecmp(id2,"-15")))||((!strcasecmp(id1,"-15"))&&(!strcasecmp(id2,"15")))) param->BRh0tautau_SM=atof(dummy);
							else if(((!strcasecmp(id1,"13"))&&(!strcasecmp(id2,"-13")))||((!strcasecmp(id1,"-13"))&&(!strcasecmp(id2,"13")))) param->BRh0mumu_SM=atof(dummy);
							else if(((!strcasecmp(id1,"3"))&&(!strcasecmp(id2,"-3")))||((!strcasecmp(id1,"-3"))&&(!strcasecmp(id2,"3")))) param->BRh0ss_SM=atof(dummy);
							else if(((!strcasecmp(id1,"4"))&&(!strcasecmp(id2,"-4")))||((!strcasecmp(id1,"-4"))&&(!strcasecmp(id2,"4")))) param->BRh0cc_SM=atof(dummy);
							else if(((!strcasecmp(id1,"6"))&&(!strcasecmp(id2,"-6")))||((!strcasecmp(id1,"-6"))&&(!strcasecmp(id2,"6")))) param->BRh0tt_SM=atof(dummy);
							else if(((!strcasecmp(id1,"24"))&&(!strcasecmp(id2,"-24")))||((!strcasecmp(id1,"-24"))&&(!strcasecmp(id2,"24")))) param->BRh0WW_SM=atof(dummy);
							else if(((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"23")))||((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"22")))) param->BRh0gaZ_SM=atof(dummy);
							else if((!strcasecmp(id1,"22"))&&(!strcasecmp(id2,"22"))) param->BRh0gaga_SM=atof(dummy);
							else if((!strcasecmp(id1,"21"))&&(!strcasecmp(id2,"21"))) param->BRh0gg_SM=atof(dummy);
							else if((!strcasecmp(id1,"23"))&&(!strcasecmp(id2,"23"))) param->BRh0ZZ_SM=atof(dummy);
						}
					if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
					}	
					break;
				}
				
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);

	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}
