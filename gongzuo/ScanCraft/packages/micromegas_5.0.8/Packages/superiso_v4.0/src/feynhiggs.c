#include "include.h"
#include "include_dm.h"
#include "feynhiggs.h"


int FeynHiggs(char name[], struct parameters* param)
{
	FILE *input,*output;
	char dummy[500],dum[5],namemod[200];
	int fail;
	
	if((param->width_h0!=0.)||(param->width_H0!=0.)||(param->width_A0!=0.)||(param->width_H!=0.))
	{
		sprintf(namemod,"%s.fhmod",name);
	
		input=fopen(name,"r");
		output=fopen(namemod,"w");
		
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
	}
	
	char tmp_char[500];
	
	sprintf(tmp_char,"%s %s 43423012 > %s.fhtmp", FEYNHIGGS, name, name);
	fail=system(tmp_char);
		
	sprintf(tmp_char,"%s.fh-001", name);
	
	
	if((fail)||(!test_file(tmp_char)))
	{	
		sprintf(tmp_char,"rm -f %s.fhtmp %s.fh-001 %s.fhmod", name, name, name);
 		system(tmp_char);
		printf("Problem with FeynHiggs, using NLO level...\n");
		return FeynHiggsNLO(name,param);
	}
		
	SM_Decays_Reader(tmp_char,param);
	Higgs_Decays_Reader(tmp_char,param);

	/*input = fopen(tmp_char,"r");
	
	while(EOF != fscanf(input,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(input,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 25: fscanf(input,"%lf",&param->mass_h0); break;
					case 26: fscanf(input,"%lf",&param->mass_h0); break;
					case 35: fscanf(input,"%lf",&param->mass_H0); break;
					case 36: fscanf(input,"%lf",&param->mass_A0); break;
					case 37: fscanf(input,"%lf",&param->mass_H); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(input,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
		}
		else if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(input,"%s",dummy); 
			switch(abs(atoi(dummy)))
			{
				case 25: fscanf(input,"%lf",&param->width_h0); break;
				case 35: fscanf(input,"%lf",&param->width_H0); break;
				case 36: fscanf(input,"%lf",&param->width_A0); break;
				case 37: fscanf(input,"%lf",&param->width_H); break;
			}		
		}
	}
	fclose(input);
	*/
		
	sprintf(tmp_char,"rm -f %s.fhtmp %s.fh-001 %s.fhmod", name, name, name);
  	system(tmp_char);

	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.)	
	{	
		printf("Problem with FeynHiggs, using NLO level...\n");
		return FeynHiggsNLO(name,param);
	}
	else return 1;
}

/* -------------------------------- */

int FeynHiggsNLO(char name[], struct parameters* param)
{
	FILE *input,*output;
	char dummy[500],dum[5],namemod[200];
	int fail;
	
	if((param->width_h0!=0.)||(param->width_H0!=0.)||(param->width_A0!=0.)||(param->width_H!=0.))
	{
		sprintf(namemod,"%s.fhmod",name);
	
		input=fopen(name,"r");
		output=fopen(namemod,"w");
		
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
	}
	
	char tmp_char[500];
	
	sprintf(tmp_char,"%s %s 43413012 > %s.fhtmp", FEYNHIGGS, name, name);
	fail=system(tmp_char);
	
	sprintf(tmp_char,"%s.fh-001", name);
	
	if((fail)||(!test_file(tmp_char)))
	{	
		sprintf(tmp_char,"rm -f %s.fhtmp %s.fh-001 %s.fhmod", name, name, name);
  		system(tmp_char);
		printf("Problem with FeynHiggs, using Tree level...\n");
		return FeynHiggsTree(name,param);
	}
	
	SM_Decays_Reader(tmp_char,param);
	Higgs_Decays_Reader(tmp_char,param);
	
	/*input = fopen(tmp_char,"r");
	
	while(EOF != fscanf(input,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(input,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				switch(atoi(dummy)*((atoi(dummy)-atof(dummy))==0.))
				{
					case 25: fscanf(input,"%lf",&param->mass_h0); break;
					case 26: fscanf(input,"%lf",&param->mass_h0); break;
					case 35: fscanf(input,"%lf",&param->mass_H0); break;
					case 36: fscanf(input,"%lf",&param->mass_A0); break;
					case 37: fscanf(input,"%lf",&param->mass_H); break;
				}
			}	
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(input,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
		}
		else if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(input,"%s",dummy); 
			switch(abs(atoi(dummy)))
			{
				case 25: fscanf(input,"%lf",&param->width_h0); break;
				case 35: fscanf(input,"%lf",&param->width_H0); break;
				case 36: fscanf(input,"%lf",&param->width_A0); break;
				case 37: fscanf(input,"%lf",&param->width_H); break;
			}		
		}
	}
	fclose(input);
	*/
		
	sprintf(tmp_char,"rm -f %s.fhtmp %s.fh-001 %s.fhmod", name, name, name);
 	system(tmp_char);

	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.)	
	{	
		printf("Problem with FeynHiggs, using Tree level...\n");
		return FeynHiggsTree(name,param);
	}
	else return 2;
}

/* -------------------------------- */

int FeynHiggsTree(char name[], struct parameters* param)
{
	FILE *input,*output;
	char dummy[500],dum[5],namemod[200];
	int fail;
	
	if((param->width_h0!=0.)||(param->width_H0!=0.)||(param->width_A0!=0.)||(param->width_H!=0.))
	{
		sprintf(namemod,"%s.fhmod",name);
	
		input=fopen(name,"r");
		output=fopen(namemod,"w");
		
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
	}
	
	char tmp_char[500];
	
	sprintf(tmp_char,"%s %s 43403012 > %s.fhtmp", FEYNHIGGS, name, name);
	
	fail=system(tmp_char);
	
	sprintf(tmp_char,"%s.fh-001", name);
	if(!test_file(tmp_char)) 
	{
		sprintf(tmp_char,"rm -f %s.fhtmp %s.fh-001 %s.fhmod", name, name, name);
  		system(tmp_char);
		return 0;
	}

	SM_Decays_Reader(tmp_char,param);
	Higgs_Decays_Reader(tmp_char,param);

	/*input = fopen(tmp_char,"r");

	while(EOF != fscanf(input,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"DECAY"))
		{
			fscanf(input,"%s",dummy); 
			switch(abs(atoi(dummy)))
			{
				case 25: fscanf(input,"%lf",&param->width_h0); break;
				case 26: fscanf(input,"%lf",&param->width_h0); break;
				case 35: fscanf(input,"%lf",&param->width_H0); break;
				case 36: fscanf(input,"%lf",&param->width_A0); break;
				case 37: fscanf(input,"%lf",&param->width_H); break;
			}		
		}
	}
	fclose(input);
	*/
		
	sprintf(tmp_char,"rm -f %s.fhtmp %s.fh-001 %s.fhmod", name, name, name);
  	system(tmp_char);

	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.) return 0; else return 3;
}
