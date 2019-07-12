#include "src/include.h"
/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/*--------------------------------------------------------*/

int main()
{ 
	struct parameters param;
		
	Init_param(&param);
	
	slha_adjust(&param);
	param.SM=1;
	
	printf("\n");
	
	printf("SuperIso v4.0 - F. Mahmoudi\n\n");
	printf("Standard Model predictions\n\n");
	printf("Observable\t\t\tValue\n\n");

	double complex C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3];
	CQ0b[1]=CQ0b[2]=CQ1b[1]=CQ1b[2]=CQpb[1]=CQpb[2]=0.;
	double obs[Nobs_BKsll+1];

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,&param));
	
	double lambda_h=0.5;
	double mu_spec=sqrt(lambda_h*param.mass_b);		
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0(C0b,C0spec,C1b,C1spec,Cpb,&param,mu_b,mu_spec,lambda_h));
	
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	
	printf("BR(B0->K* gamma)\t\t%.3e\n",BR_BKstargamma(0,C0b,C1b,C2b,Cpb,&param,mu_b));
	printf("BR(B+->K* gamma)\t\t%.3e\n\n",BR_BKstargamma(1,C0b,C1b,C2b,Cpb,&param,mu_b));
	
	printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));

 	printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu(&param));
    printf("R(B->tau nu)\t\t\t%.3e\n\n",RBtaunu(&param));
    printf("BR(B->D tau nu)\t\t\t%.3e\n",BRBDlnu_full(3,1,obs,&param));
    printf("AFB(B->D tau nu)\t\t%.3e\n",obs[1]);
    printf("Ptau(B->D tau nu)\t\t%.3e\n",obs[2]);
	printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n\n",BDtaunu_BDenu(&param));
    printf("BR(B->D* tau nu)\t\t%.3e\n",BRBDstarlnu_full(3,1,obs,&param));
    printf("AFB(B->D* tau nu)\t\t%.3e\n",obs[1]);
    printf("Ptau(B->D* tau nu)\t\t%.3e\n",obs[2]);
    printf("PD*(B->D* tau nu)\t\t%.3e\n",obs[3]);
	printf("BR(B->D* tau nu)/BR(B->D* e nu)\t%.3e\n\n",BDstartaunu_BDstarenu(&param));
	printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu(&param));
	printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu(&param));
	printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu(&param));
	printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu(&param));
	printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23(&param));
									
	return 1;
}
