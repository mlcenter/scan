#include "include.h"
#include "include_dm.h"
#include "sdecay.h"

int Sdecay(char name[], struct parameters* param)
{
	FILE *tmp,*tmp2;
	char tmp_char[500],namedir[300];
	int fail;
	char *curdir;
	int dum;
	
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.sdtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"%s/%s",namedir,"susyhit.in");
	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"* SUspect-SdecaY-Hdecay-InTerface options:\n");
	fprintf(tmp,"------------------------------------------\n");
	fprintf(tmp,"  (1) link SuSpect-Sdecay-Hdecay and take the input parameters from SuSpect\n");
	fprintf(tmp,"  (2) link Sdecay-Hdecay, not SuSpect, take the input parameters from any SLHA\n");
	fprintf(tmp,"      input, called slhaspectrum.in or slhaspectrumFV.in in the FV case\n");
	fprintf(tmp,"2\n");
	fprintf(tmp,"* Choice of the output, SLHA format (1) or simple (0):\n");
	fprintf(tmp,"------------------------------------------------------\n");
	fprintf(tmp,"1\n");
	fprintf(tmp,"* HDECAY input parameters:\n");
	fprintf(tmp,"--------------------------\n");
	fprintf(tmp,"MSBAR(1) = %.5e\n",param->mass_s);
	fprintf(tmp,"MC       = %.5e\n",param->mass_c);
	fprintf(tmp,"MMUON    = %.5e\n",param->mass_mu);
	fprintf(tmp,"1/ALPHA  = %.5e\n",137.0359895e0);
	fprintf(tmp,"GAMW     = %.5e\n",param->width_W);
	fprintf(tmp,"GAMZ     = %.5e\n",param->width_Z);
	fprintf(tmp,"VUS      = %.5e\n",cabs(param->Vus));
	fprintf(tmp,"VCB      = %.5e\n",cabs(param->Vcb));
	fprintf(tmp,"VUB/VCB  = %.5e\n",cabs(param->Vub/param->Vcb));
	fprintf(tmp,"* FCNC stop decays for SDECAY:\n");
	fprintf(tmp,"----------------------------------------\n");
	fprintf(tmp,"  (0) No flavour violation at tree-level.\n");
	fprintf(tmp,"  (1) FCNC stop decays in SDECAY. Works only with output in SLHA format, and input parameters have to be provided in an SLHA2 file.\n");
	fprintf(tmp,"0\n");
	fprintf(tmp,"  (0) Light stop four-body decay widths summed up.\n");
	fprintf(tmp,"  (1) Stop four-body decay separated by final states.\n");
	fprintf(tmp,"0\n");

	fclose(tmp);
	
	tmp=fopen(name,"r");
	sprintf(tmp_char,"%s/%s",namedir,"slhaspectrum.in");
	tmp2=fopen(tmp_char,"w");
	while((dum=getc(tmp))!=EOF) putc(dum,tmp2);
	fclose(tmp);
	fclose(tmp2);

	chdir(namedir);
	
	sprintf(tmp_char,"%s > tmp",SDECAY);

	fail=system(tmp_char);

	if((fail)||(!test_file("susyhit_slha.out")))	
	{	
		chdir(curdir);
		sprintf(tmp_char,"rm -rf %s",namedir);
 		system(tmp_char);
		printf("Sdecay failed\n");
		return 0;
	}
	
	SUSY_Decays_Reader("susyhit_slha.out",param);

	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}
