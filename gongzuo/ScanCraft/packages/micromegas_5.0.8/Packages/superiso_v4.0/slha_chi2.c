#define NBOBSMAX 500

#include "src/include.h"

/*--------------------------------------------------------*/
/* Calculation of the observables using a given SLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	int test;

  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the SLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}

	printf("\n");
	
	printf("SuperIso v4.0 - F. Mahmoudi\n\n");
	printf("SLHA input file\n\n");

	test=test_slha(name);
	
	if(test>0)
	{
		if(test==2) printf("WARNING: only tested in the MFV scenario!\n\n");
	}
	else 
	{
		if(test==-1) printf("Invalid point\n\n");
		else if(test==-2) printf("Model not yet implemented\n\n");
		else if(test==-3) printf("Invalid SLHA file\n\n");
		else if(test==-4) printf("SLHA file absent\n\n");
		exit(1);
	}

	struct parameters param;
	int ie=0;

	double values[NBOBSMAX];
	char* names[NBOBSMAX];
	
	int nmax;

	double complex C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3];
	double complex C0eb[11],C1eb[11],C0ew[11],C1ew[11],C2ew[11],C2eb[11],Cpeb[11],CQpeb[3],CQ0eb[3],CQ1eb[3];
	double obsKstar[Nobs_BKsll+1];
	double obsK[Nobs_BKll+1];
	double obsphi[Nobs_Bsphill+1];

	double mu_W,mu_b,mu_spec;	
	double lambda_h=0.5;

	Init_param(&param);
	
	Les_Houches_Reader(name,&param);
	
	//param.fullFF=0; /* uncomment to use soft FF approach for the inclusive b->s ll observables - full FF approach by default */
	//param.likelihoodBKstarmumu=0; /* 0 = LHCb method of moments / 1 = LHCb likelihood fit for B -> K* mu+ mu- observables */

	mu_W=2.*param.mass_W;	
	mu_b=param.mass_b_pole/2.;
	mu_spec=sqrt(lambda_h*param.mass_b);

	CW_calculator(2,C0w,C1w,C2w,mu_W,&param); /* 2 = muon */

	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	values[ie]=delta0(C0b,C0spec,C1b,C1spec,Cpb,&param,mu_b,mu_spec,lambda_h);
	names[ie++]="delta0";
	
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	values[ie]=bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,&param);
	names[ie++]="bsgamma";

	mu_W=param.mass_W;	
	mu_b=param.mass_b_pole;

	CW_calculator(2,C0w,C1w,C2w,mu_W,&param); /* 2 = muon */
	CW_calculator(1,C0ew,C1ew,C2ew,mu_W,&param); /* 1 = electron */

	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);

	C_calculator_base1(C0ew,C1ew,C2ew,mu_W,C0eb,C1eb,C2eb,mu_b,&param);
	Cprime_calculator(1,Cpeb,CQpeb,mu_W,mu_b,&param);
	CQ_calculator(1,CQ0eb,CQ1eb,mu_W,mu_b,&param);
	
	values[ie]=Bsmumu_untag(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
	names[ie++]="Bsmumu_untag";
		
//	values[ie]=Bsll_untag(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
//	names[ie++]="Bsee_untag";
	
	values[ie]=Bdmumu(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
	names[ie++]="Bdmumu";

	values[ie]=BRBXsll_lowq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
  	names[ie++]="BRBXsmumu_lowq2";
 	values[ie]=BRBXsll_highq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
  	names[ie++]="BRBXsmumu_highq2";
//	values[ie]=BRBXsll(2,0.045/param.mass_b_1S/param.mass_b_1S,0.999,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
//	names[ie++]="BRBXsmumu_full";
	
	values[ie]=BRBXsll_lowq2(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
  	names[ie++]="BRBXsee_lowq2";
	values[ie]=BRBXsll_highq2(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
  	names[ie++]="BRBXsee_highq2";
//	values[ie]=BRBXsll(1,0.045/param.mass_b_1S/param.mass_b_1S,0.999,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
// 	names[ie++]="BRBXsee_full";

	values[ie]=BR_BKstargamma(0,C0b,C1b,C2b,Cpb,&param,mu_b);
 	names[ie++]="BR_B0Kstar0gamma";

	values[ie]=BRBKstarll(2,0,0.1,0.98,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(0.98-0.1);
 	names[ie++]="dG_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_0.1-0.98";
/*	values[ie]=obsKstar[18];
 	names[ie++]="P5p_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[5];
	names[ie++]="P1_B0Kstar0mumu_0.1_0.98";
	values[ie]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	names[ie++]="ATRe_B0Kstar0mumu_0.1_0.98";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_0.1-0.98";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_0.1-0.98";*/

if(param.likelihoodBKstarmumu==1)
{
	values[ie]=BRBKstarll(2,0,1.1,2.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.5-1.1);
 	names[ie++]="dG_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_1.1-2.5";
/*	values[ie]=obsKstar[18];
	names[ie++]="P5p_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_1.1-2.5";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_1.1-2.5";*/
 
	values[ie]=BRBKstarll(2,0,2.5,4.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-2.5);
	names[ie++]="dG_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_2.5-4";
/*	values[ie]=obsKstar[18];
 	names[ie++]="P5p_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_2.5-4";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_2.5-4";*/
 
	values[ie]=BRBKstarll(2,0,4.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
	names[ie++]="dG_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_4-6";
/*	values[ie]=obsKstar[18];
 	names[ie++]="P5p_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_4-6";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_4-6";*/

	values[ie]=BRBKstarll(2,0,6.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
 	names[ie++]="dG_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_6-8";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_6-8";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_6-8";*/

/*	values[ie]=BRBKstarll(2,0,11.,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
 	names[ie++]="dG_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_11-12.5";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_11-12.5";*/

	values[ie]=BRBKstarll(2,0,15.,17.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
 	names[ie++]="dG_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[28];
	names[ie++]="S7_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_15-17";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_15-17";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_15-17";*/

	values[ie]=BRBKstarll(2,0,17.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-17.);
	names[ie++]="dG_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_17-19";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_17-19";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_17-19";*/

/*	values[ie]=BRBKstarll(2,0,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
 	names[ie++]="dG_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_1.1-6";
	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_1.1-6";
 	values[ie]=obsKstar[5];
	names[ie++]="P1_B0Kstar0mumu_1.1_6";
	values[ie]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	names[ie++]="ATRe_B0Kstar0mumu_1.1_6";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_1.1_6";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_1.1_6";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_1.1_6";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_1.1_6";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_1.1_6";*/
  	
/*	values[ie]=BRBKstarll(2,0,15.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-15);
  	names[ie++]="dG_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_15-19";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_15-19";*/
}
else
{
	values[ie]=BRBKstarll(2,0,1.1,2.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.5-1.1);
 	names[ie++]="dG_B0Kstar0mumu_1.1-2.5";
	values[ie]=BRBKstarll(2,0,2.5,4.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-2.5);
	names[ie++]="dG_B0Kstar0mumu_2.5-4";
	values[ie]=BRBKstarll(2,0,4.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
	names[ie++]="dG_B0Kstar0mumu_4-6";
	values[ie]=BRBKstarll(2,0,6.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
 	names[ie++]="dG_B0Kstar0mumu_6-8";
/*	values[ie]=BRBKstarll(2,0,11.,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
 	names[ie++]="dG_B0Kstar0mumu_11-12.5";*/
	values[ie]=BRBKstarll(2,0,15.,17.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
 	names[ie++]="dG_B0Kstar0mumu_15-17";
	values[ie]=BRBKstarll(2,0,17.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-17.);
	names[ie++]="dG_B0Kstar0mumu_17-19";
/*	values[ie]=BRBKstarll(2,0,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
 	names[ie++]="dG_B0Kstar0mumu_1.1-6";
	values[ie]=BRBKstarll(2,0,15.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-15);
  	names[ie++]="dG_B0Kstar0mumu_15-19";*/

	values[ie]=BRBKstarll(2,0,1.1,2.0,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-1.1);
// 	names[ie++]="dG_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_1.1-2";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_1.1-2";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_1.1-2";*/

	values[ie]=BRBKstarll(2,0,2.,3.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(3.-2.);
// 	names[ie++]="dG_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_2-3";
/*	values[ie]=obsKstar[18];
 	names[ie++]="P5p_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_2-3";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_2-3";*/

	values[ie]=BRBKstarll(2,0,3.,4.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-3.);
 //	names[ie++]="dG_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_3-4";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_3-4";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_3-4";*/

	values[ie]=BRBKstarll(2,0,4.,5.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(5.-4.);
// 	names[ie++]="dG_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_4-5";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_4-5";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_4-5";*/

	values[ie]=BRBKstarll(2,0,5.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-5.);
// 	names[ie++]="dG_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_5-6";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_5-6";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_5-6";*/

	values[ie]=BRBKstarll(2,0,6.,7.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(7.-6.);
// 	names[ie++]="dG_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_6-7";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_6-7";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_6-7";*/

	values[ie]=BRBKstarll(2,0,7.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-7.);
// 	names[ie++]="dG_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_7-8";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_7-8";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_7-8";*/

/*	values[ie]=BRBKstarll(2,0,11.,11.75,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(11.75-11.);
 	names[ie++]="dG_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_11-11.75";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_11-11.75";*/

/*	values[ie]=BRBKstarll(2,0,11.75,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.75);
 	names[ie++]="dG_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_11.75-12.5";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_11.75-12.5";*/

	values[ie]=BRBKstarll(2,0,15.,16.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(16.-15.);
// 	names[ie++]="dG_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[28];
	names[ie++]="S7_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_15-16";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_15-16";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_15-16";*/

	values[ie]=BRBKstarll(2,0,16.,17.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-16.);
// 	names[ie++]="dG_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[2];
  	names[ie++]="FL_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[28];
	names[ie++]="S7_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_16-17";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_16-17";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_16-17";*/

	values[ie]=BRBKstarll(2,0,17.,18.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(18.-17.);
// 	names[ie++]="dG_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_17-18";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_17-18";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_17-18";*/

	values[ie]=BRBKstarll(2,0,18.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-18.);
// 	names[ie++]="dG_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[1];
  	names[ie++]="AFB_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[25];
  	names[ie++]="S3_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[26];
  	names[ie++]="S4_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[27];
  	names[ie++]="S5_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[28];
  	names[ie++]="S7_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[29];
  	names[ie++]="S8_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[30];
  	names[ie++]="S9_B0Kstar0mumu_18-19";
/*	values[ie]=obsKstar[18];
   	names[ie++]="P5p_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[5];
   	names[ie++]="P1_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[14];
   	names[ie++]="P2_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[15];
   	names[ie++]="P3_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[17];
   	names[ie++]="P4p_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[19];
   	names[ie++]="P6p_B0Kstar0mumu_18-19";
	values[ie]=obsKstar[21];
   	names[ie++]="P8p_B0Kstar0mumu_18-19";*/
}

/*	values[ie]=BRBKstarll(2,1,0.1,2.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-0.1);
  	names[ie++]="dG_B+Kstar+mumu_0.1-2";
	values[ie]=BRBKstarll(2,1,2.,4.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-2.);
  	names[ie++]="dG_B+Kstar+mumu_2-4";
	values[ie]=BRBKstarll(2,1,4.,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
  	names[ie++]="dG_B+Kstar+mumu_4-6";
	values[ie]=BRBKstarll(2,1,6.,8.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
  	names[ie++]="dG_B+Kstar+mumu_6-8";
	values[ie]=BRBKstarll(2,1,11.,12.5,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
  	names[ie++]="dG_B+Kstar+mumu_11-12.5";
	values[ie]=BRBKstarll(2,1,15.,17.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
  	names[ie++]="dG_B+Kstar+mumu_15-17";
	values[ie]=BRBKstarll(2,1,17.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-17.);
  	names[ie++]="dG_B+Kstar+mumu_17-19";*/
	
	values[ie]=BRBKstarll(2,1,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
  	names[ie++]="dG_B+Kstar+mumu_1.1-6";
	values[ie]=BRBKstarll(2,1,15.,19.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-15.);
  	names[ie++]="dG_B+Kstar+mumu_15-19";
	
	values[ie]=BRBKstarll(1,0,0.1,pow(param.m_Bd-param.m_Kstar0,2.)*0.999,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b);
  	names[ie++]="BR_B0Kstar0ee_full";

//	values[ie]=BRBKstarll(1,0,0.0009,1.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(1.-0.0009);
// 	names[ie++]="dG_B0Kstar0ee_0.0009-1";
                
/*	values[ie]=BRBKstarll(1,0,0.002,1.12,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(1.12-0.002);
  	names[ie++]="dG_B0Kstar0ee_0.002-1.12";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0ee_0.002_1.12";
	values[ie]=obsKstar[5];
	names[ie++]="AT2_B0Kstar0ee_0.002_1.12";
	values[ie]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	names[ie++]="ATRe_B0Kstar0ee_0.002_1.12";*/
                
/*	values[ie]=BRBKstarll(1,0,1.1,6.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(6.-1.1);
  	names[ie++]="dG_B0Kstar0ee_1.1-6";
	values[ie]=obsKstar[2];
	names[ie++]="FL_B0Kstar0ee_1.1_6";
	values[ie]=obsKstar[5];
	names[ie++]="AT2_B0Kstar0ee_1.1_6";
	values[ie]=(4.*obsKstar[1])/(3.*(1.-obsKstar[2]));
	names[ie++]="ATRe_B0Kstar0ee_1.1_6";*/

//	values[ie]=BRBKstarll(1,0,1.,6.,obsKstar,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(6.-1.);
// 	names[ie++]="dG_B0Kstar0ee_1-6";

	values[ie]=BRBKstarll(2,0,0.045,1.1,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/BRBKstarll(1,0,0.045,1.1,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)-1.;
	names[ie++]="RKstar_0.045-1.1-1";
	values[ie]=BRBKstarll(2,0,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/BRBKstarll(1,0,1.1,6.,obsKstar,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)-1.;
	names[ie++]="RKstar_1.1-6-1";

/*	values[ie]=BRBKll(2,0,0.1,2.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-0.1);
  	names[ie++]="dG_B0K0mumu_0.1-2";
	values[ie]=BRBKll(2,0,2.,4.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-2.);
  	names[ie++]="dG_B0K0mumu_2-4";
	values[ie]=BRBKll(2,0,4.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-4.);
  	names[ie++]="dG_B0K0mumu_4-6";
	values[ie]=BRBKll(2,0,6.,8.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-6.);
  	names[ie++]="dG_B0K0mumu_6-8";
	values[ie]=BRBKll(2,0,11.,12.5,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
  	names[ie++]="dG_B0K0mumu_11-12.5";
	values[ie]=BRBKll(2,0,15.,17.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
  	names[ie++]="dG_B0K0mumu_15-17";
	values[ie]=BRBKll(2,0,17.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-17.);
  	names[ie++]="dG_B0K0mumu_17-22";*/
	
	values[ie]=BRBKll(2,0,1.1,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
  	names[ie++]="dG_B0K0mumu_1.1-6";
	values[ie]=BRBKll(2,0,15.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-15.);
  	names[ie++]="dG_B0K0mumu_15-22";


/*	values[ie]=BRBKll(2,1,0.1,0.98,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(0.98-0.1);	
 	names[ie++]="dG_B+K+mumu_0.1-0.98";
	values[ie]=BRBKll(2,1,1.1,2.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-1.1);	
 	names[ie++]="dG_B+K+mumu_1.1-2";
	values[ie]=BRBKll(2,1,2.,3.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(3.-2.);
  	names[ie++]="dG_B+K+mumu_2-3";
	values[ie]=BRBKll(2,1,3,4.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(4.-3.);
  	names[ie++]="dG_B+K+mumu_3-4";
	values[ie]=BRBKll(2,1,4.,5.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(5.-4.);
  	names[ie++]="dG_B+K+mumu_4-5";
	values[ie]=BRBKll(2,1,5.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-5.);
  	names[ie++]="dG_B+K+mumu_5-6";
	values[ie]=BRBKll(2,1,6.,7.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(7.-6.);
  	names[ie++]="dG_B+K+mumu_6-7";
	values[ie]=BRBKll(2,1,7.,8.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-7.);
  	names[ie++]="dG_B+K+mumu_7-8";
	values[ie]=BRBKll(2,1,11.,11.8,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(11.8-11.);
  	names[ie++]="dG_B+K+mumu_11-11.8";
	values[ie]=BRBKll(2,1,11.8,12.5,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.8);
  	names[ie++]="dG_B+K+mumu_11.8-12.5";
	values[ie]=BRBKll(2,1,15.,16.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(16.-15.);
  	names[ie++]="dG_B+K+mumu_15-16";
	values[ie]=BRBKll(2,1,16.,17.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-16.);
  	names[ie++]="dG_B+K+mumu_16-17";
	values[ie]=BRBKll(2,1,17.,18.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(18.-17.);
  	names[ie++]="dG_B+K+mumu_17-18";
	values[ie]=BRBKll(2,1,18.,19.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(19.-18.);
  	names[ie++]="dG_B+K+mumu_18-19";
	values[ie]=BRBKll(2,1,19.,20.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(20.-19.);
  	names[ie++]="dG_B+K+mumu_19-20";
	values[ie]=BRBKll(2,1,20.,21.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(21.-20.);
  	names[ie++]="dG_B+K+mumu_20-21";
	values[ie]=BRBKll(2,1,21.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-21.);
  	names[ie++]="dG_B+K+mumu_21-22";*/

	values[ie]=BRBKll(2,1,1.1,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.1);
  	names[ie++]="dG_B+K+mumu_1.1-6";
	values[ie]=BRBKll(2,1,15.,22.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(22.-15.);
  	names[ie++]="dG_B+K+mumu_15-22";

//	values[ie]=BRBKll(2,1,1.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.);
// 	names[ie++]="dG_B+K+mumu_1-6";

//	values[ie]=BRBKll(1,1,1.,6.,obsK,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)/(6.-1.);
// 	names[ie++]="dG_B+K+ee_1-6";

	values[ie]=BRBKll(2,1,1.,6.,obsK,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/BRBKll(1,1,1.,6.,obsK,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,&param,mu_b)-1.;
  	names[ie++]="RK_1-6-1";

	values[ie]=BRBsphill(2,0.1,2.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(2.-0.1);
 	names[ie++]="dG_Bsphimumu_0.1-2";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_0.1-2";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_0.1-2";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_0.1-2";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_0.1-2";

	values[ie]=BRBsphill(2,2.,5.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(5.-2.);
 	names[ie++]="dG_Bsphimumu_2-5";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_2-5";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_2-5";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_2-5";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_2-5";

	values[ie]=BRBsphill(2,5.,8.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(8.-5.);
 	names[ie++]="dG_Bsphimumu_5-8";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_5-8";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_5-8";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_5-8";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_5-8";

/*	values[ie]=BRBsphill(2,11.,12.5,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(12.5-11.);
 	names[ie++]="dG_Bsphimumu_11-12.5";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_11-12.5";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_11-12.5";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_11-12.5";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_11-12.5";*/

	values[ie]=BRBsphill(2,15.,17.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(17.-15.);
 	names[ie++]="dG_Bsphimumu_15-17";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_15-17";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_15-17";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_15-17";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_15-17";

	values[ie]=BRBsphill(2,17.,pow(param.m_Bs-param.m_phi,2.)*0.999,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(pow(param.m_Bs-param.m_phi,2.)*0.999-17.);
 	names[ie++]="dG_Bsphimumu_17-19";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_17-19";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_17-19";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_17-19";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_17-19";

/*	values[ie]=BRBsphill(2,1.,6.,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(6.-1.);
 	names[ie++]="dG_Bsphimumu_1-6";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_1-6";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_1-6";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_1-6";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_1-6";*/

/*	values[ie]=BRBsphill(2,15.,pow(param.m_Bs-param.m_phi,2.)*0.999,obsphi,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b)/(pow(param.m_Bs-param.m_phi,2.)*0.999-15.);
 	names[ie++]="dG_Bsphimumu_15-19";
	values[ie]=obsphi[1];
   	names[ie++]="FL_Bsphimumu_15-19";
	values[ie]=obsphi[3];
   	names[ie++]="S3_Bsphimumu_15-19";
	values[ie]=obsphi[4];
   	names[ie++]="S4_Bsphimumu_15-19";
	values[ie]=obsphi[6];
   	names[ie++]="S7_Bsphimumu_15-19";*/
   	  	   	
   	nmax=ie;
   	
	for(ie=0;ie<nmax;ie++) printf("%d\t%-25s:\t %.3e\n",ie,names[ie],values[ie]);
		
	int nbobs;
	printf("chi2 = %.2f\n",chi2(names,nmax,&param,&nbobs));
	printf("n_obs=%d\n\n",nbobs);
	
	return 1;
}
