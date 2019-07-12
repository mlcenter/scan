#include "src/include.h"
/*--------------------------------------------------------*/
/* Calculation of the observables using a given FLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{ 
	struct parameters param;
	struct flhaparam paramflha;
	char name[50];
		
  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the FLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}
	
	printf("\n");
	
	printf("SuperIso v4.0 - F. Mahmoudi\n\n");
	printf("FLHA input file\n\n");

	Init_param(&param);
	
	Les_Houches_Reader(name,&param);
	slha_adjust(&param);
	param.SM=1;

	FLHA_Reader(name,&paramflha);
	
	printf("\n");
	
	printf("Observable\t\t\tValue\n\n");

	double complex C0b[11],C1b[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3];
	double complex C0eb[11],C1eb[11],C0ew[11],C1ew[11],C2ew[11],C2eb[11],Cpeb[11],CQpeb[3],CQ0eb[3],CQ1eb[3];
	CQ0b[1]=CQ0b[2]=CQ1b[1]=CQ1b[2]=CQpb[1]=CQpb[2]=0.;
	CQ0eb[1]=CQ0eb[2]=CQ1eb[1]=CQ1eb[2]=CQpeb[1]=CQpeb[2]=0.;

	double complex C0bSM[11],C1bSM[11],C0wSM[11],C1wSM[11],C2wSM[11],C2bSM[11],CpbSM[11],CQpbSM[3],CQ0bSM[3],CQ1bSM[3];
	double complex C0ebSM[11],C1ebSM[11],C0ewSM[11],C1ewSM[11],C2ewSM[11],C2ebSM[11],CpebSM[11],CQpebSM[3],CQ0ebSM[3],CQ1ebSM[3];
	CQ0bSM[1]=CQ0bSM[2]=CQ1bSM[1]=CQ1bSM[2]=CQpbSM[1]=CQpbSM[2]=0.;
	CQ0ebSM[1]=CQ0ebSM[2]=CQ1ebSM[1]=CQ1ebSM[2]=CQpebSM[1]=CQpebSM[2]=0.;

	double mu_W=2.*param.mass_W;
	if(paramflha.Q>50.) mu_W=paramflha.Q;

	double mu_b=param.mass_b_pole/2.;
		
	if(paramflha.dC7[0]!=0.||paramflha.dC8[0]!=0.||paramflha.dC9[0]!=0.||paramflha.dC10[0]!=0.||paramflha.dC9e[0]!=0.||paramflha.dC10e[0]!=0.||paramflha.dC7p[0]!=0.||paramflha.dC8p[0]!=0.||paramflha.dC9p[0]!=0.||paramflha.dC10p[0]!=0.||paramflha.dC9pe[0]!=0.||paramflha.dC10pe[0]!=0.)
	{	
		CW_calculator(2,C0w,C1w,C2w,mu_W,&param);	
		C0w[7]+=paramflha.dC7[0];
		C1w[7]+=paramflha.dC7[1];
		C2w[7]+=paramflha.dC7[2];
		C0w[8]+=paramflha.dC8[0];
		C1w[8]+=paramflha.dC8[1];
		C2w[8]+=paramflha.dC8[2];
		C0w[9]+=paramflha.dC9[0];
		C1w[9]+=paramflha.dC9[1];
		C2w[9]+=paramflha.dC9[2];
		C0w[10]+=paramflha.dC10[0];
		C1w[10]+=paramflha.dC10[1];
		C2w[10]+=paramflha.dC10[2];

		CW_calculator(1,C0ew,C1ew,C2ew,mu_W,&param);
		C0ew[9]+=paramflha.dC9e[0];
		C1ew[9]+=paramflha.dC9e[1];
		C2ew[9]+=paramflha.dC9e[2];
		C0ew[10]+=paramflha.dC10e[0];
		C1ew[10]+=paramflha.dC10e[1];
		C2ew[10]+=paramflha.dC10e[2];
		
		C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
		C_calculator_base1(C0ew,C1ew,C2ew,mu_W,C0eb,C1eb,C2eb,mu_b,&param);

		Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
		Cpb[7]+=paramflha.dC7p[0];
		Cpb[8]+=paramflha.dC8p[0];
		Cpb[9]+=paramflha.dC9p[0];
		Cpb[10]+=paramflha.dC10p[0];
		
		Cprime_calculator(1,Cpeb,CQpeb,mu_W,mu_b,&param);
		Cpeb[9]+=paramflha.dC9pe[0];
		Cpeb[10]+=paramflha.dC10pe[0];
	}
	else
	{
		C0w[7]=paramflha.C7[0];
		C1w[7]=paramflha.C7[1];
		C2w[7]=paramflha.C7[2];
		C0w[8]=paramflha.C8[0];
		C1w[8]=paramflha.C8[1];
		C2w[8]=paramflha.C8[2];
		C0w[9]=paramflha.C9[0];
		C1w[9]=paramflha.C9[1];
		C2w[9]=paramflha.C9[2];
		C0w[10]=paramflha.C10[0];
		C1w[10]=paramflha.C10[1];
		C2w[10]=paramflha.C10[2];

		C0ew[9]=paramflha.C9e[0];
		C1ew[9]=paramflha.C9e[1];
		C2ew[9]=paramflha.C9e[2];
		C0ew[10]=paramflha.C10e[0];
		C1ew[10]=paramflha.C10e[1];
		C2ew[10]=paramflha.C10e[2];
		
		C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
		
		Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
		Cpb[7]=paramflha.C7p[0];
		Cpb[8]=paramflha.C8p[0];
		Cpb[9]=paramflha.C9p[0];
		Cpb[10]=paramflha.C10p[0];
		
		C_calculator_base1(C0ew,C1ew,C2ew,mu_W,C0eb,C1eb,C2eb,mu_b,&param);
		
		Cprime_calculator(1,Cpeb,CQpeb,mu_W,mu_b,&param);
		Cpeb[9]=paramflha.C9pe[0];
		Cpeb[10]=paramflha.C10pe[0];
	}
	
	CW_calculator(2,C0wSM,C1wSM,C2wSM,mu_W,&param);
	CW_calculator(1,C0ewSM,C1ewSM,C2ewSM,mu_W,&param);
	C_calculator_base1(C0wSM,C1wSM,C2wSM,mu_W,C0bSM,C1bSM,C2bSM,mu_b,&param);	
	C_calculator_base1(C0ewSM,C1ewSM,C2ewSM,mu_W,C0ebSM,C1ebSM,C2ebSM,mu_b,&param);

	Cprime_calculator(2,CpbSM,CQpbSM,mu_W,mu_b,&param);
	Cprime_calculator(1,CpebSM,CQpebSM,mu_W,mu_b,&param);

	double alphas_mub=alphas_running(mu_b,param.mass_top_pole,param.mass_b_pole,&param);
	
	double complex C9=C0b[9]+alphas_mub/4./pi*(C1b[9]+alphas_mub/4./pi*C2b[9]);
	double complex C9SM=C0bSM[9]+alphas_mub/4./pi*(C1bSM[9]+alphas_mub/4./pi*C2bSM[9]);
	double complex C10=C0b[10]+alphas_mub/4./pi*(C1b[10]+alphas_mub/4./pi*C2b[10]);
	double complex C10SM=C0bSM[10]+alphas_mub/4./pi*(C1bSM[10]+alphas_mub/4./pi*C2bSM[10]);

	double complex C9e=C0eb[9]+alphas_mub/4./pi*(C1eb[9]+alphas_mub/4./pi*C2eb[9]);
	double complex C9eSM=C0ebSM[9]+alphas_mub/4./pi*(C1ebSM[9]+alphas_mub/4./pi*C2ebSM[9]);
	double complex C10e=C0eb[10]+alphas_mub/4./pi*(C1eb[10]+alphas_mub/4./pi*C2eb[10]);
	double complex C10eSM=C0ebSM[10]+alphas_mub/4./pi*(C1ebSM[10]+alphas_mub/4./pi*C2ebSM[10]);
	
	printf("deltaC9mu\t\t\t%.3e\n",creal(C9-C9SM));
	printf("deltaC10mu\t\t\t%.3e\n\n",creal(C10-C10SM));
	
	printf("deltaC9e\t\t\t%.3e\n",creal(C9e-C9eSM));
	printf("deltaC10e\t\t\t%.3e\n\n",creal(C10-C10SM));
	
	printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,&param));
	
	return 1;
}
