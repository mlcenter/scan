central[0]=1.2e-2;
errors[0]=5.1e-2;
names[0]="delta0";

central[1]=3.32e-4;
errors[1]=0.15e-4;
names[1]="bsgamma";

central[2]=3.0e-9;
errors[2]=sqrt(stat2BRBsll*0.6e-9*0.6e-9+syst2BRBsll*0.3e-9*0.3e-9);
names[2]="Bsmumu_untag";

central[3]=0.;
errors[3]=1.7e-7;
names[3]="Bsee_untag";

central[4]=1.5e-10;
errors[4]=sqrt(stat2BRBsll*1.2e-10*1.2e-10+syst2BRBsll*0.2e-10*0.2e-10);
names[4]="Bdmumu";

central[5]=0.66e-6;
errors[5]=0.88e-6;
names[5]="BRBXsmumu_lowq2";

central[6]=0.6e-6;
errors[6]=0.31e-6;
names[6]="BRBXsmumu_highq2";

central[7]=-0.;
errors[7]=-0.;
names[7]="BRBXsmumu_full";

central[8]=1.93e-6;
errors[8]=0.55e-6;
names[8]="BRBXsee_lowq2";

central[9]=0.56e-6;
errors[9]=0.19e-6;
names[9]="BRBXsee_highq2";

central[10]=-0.;
errors[10]=-0.;
names[10]="BRBXsee_full";

central[11]=1.016e-7;
errors[11]=sqrt(stat2BRB0Kstar0*0.073e-7*0.073e-7+syst2BRB0Kstar0*0.029e-7*0.029e-7+other2BRB0Kstar0*0.069e-7*0.069e-7);
names[11]="dG_B0Kstar0mumu_0.1-0.98";

if(likelihood==1)
{
central[12]=0.263;
errors[12]=sqrt(stat2AngB0Kstar0*0.045*0.045+syst2AngB0Kstar0*0.017*0.017);
names[12]="FL_B0Kstar0mumu_0.1-0.98";

central[13]=-0.003;
errors[13]=sqrt(stat2AngB0Kstar0*0.058*0.058+syst2AngB0Kstar0*0.009*0.009);
names[13]="AFB_B0Kstar0mumu_0.1-0.98";

central[14]=-0.036;
errors[14]=sqrt(stat2AngB0Kstar0*0.063*0.063+syst2AngB0Kstar0*0.005*0.005);
names[14]="S3_B0Kstar0mumu_0.1-0.98";

central[15]=-0.082;
errors[15]=sqrt(stat2AngB0Kstar0*0.069*0.069+syst2AngB0Kstar0*0.009*0.009);
names[15]="S4_B0Kstar0mumu_0.1-0.98";

central[16]=0.17;
errors[16]=sqrt(stat2AngB0Kstar0*0.059*0.059+syst2AngB0Kstar0*0.018*0.018);
names[16]="S5_B0Kstar0mumu_0.1-0.98";

central[17]=0.015;
errors[17]=sqrt(stat2AngB0Kstar0*0.059*0.059+syst2AngB0Kstar0*0.006*0.006);
names[17]="S7_B0Kstar0mumu_0.1-0.98";

central[18]=-0.079;
errors[18]=sqrt(stat2AngB0Kstar0*0.076*0.076+syst2AngB0Kstar0*0.007*0.007);
names[18]="S8_B0Kstar0mumu_0.1-0.98";

central[19]=-0.083;
errors[19]=sqrt(stat2AngB0Kstar0*0.058*0.058+syst2AngB0Kstar0*0.004*0.004);
names[19]="S9_B0Kstar0mumu_0.1-0.98";

central[20]=0.387;
errors[20]=sqrt(stat2AngB0Kstar0*0.133*0.133+syst2AngB0Kstar0*0.052*0.052);
names[20]="P5p_B0Kstar0mumu_0.1-0.98";
}
else
{
central[12]=0.242;
errors[12]=sqrt(stat2AngB0Kstar0*0.058*0.058+syst2AngB0Kstar0*0.026*0.026);
names[12]="FL_B0Kstar0mumu_0.1-0.98";

central[13]=-0.138;
errors[13]=sqrt(stat2AngB0Kstar0*0.095*0.095+syst2AngB0Kstar0*0.072*0.072);
names[13]="AFB_B0Kstar0mumu_0.1-0.98";

central[14]=-0.014;
errors[14]=sqrt(stat2AngB0Kstar0*0.060*0.060+syst2AngB0Kstar0*0.008*0.008);
names[14]="S3_B0Kstar0mumu_0.1-0.98";

central[15]=-0.039;
errors[15]=sqrt(stat2AngB0Kstar0*0.091*0.091+syst2AngB0Kstar0*0.015*0.015);
names[15]="S4_B0Kstar0mumu_0.1-0.98";

central[16]=0.129;
errors[16]=sqrt(stat2AngB0Kstar0*0.068*0.068+syst2AngB0Kstar0*0.011*0.011);
names[16]="S5_B0Kstar0mumu_0.1-0.98";

central[17]=0.038;
errors[17]=sqrt(stat2AngB0Kstar0*0.063*0.063+syst2AngB0Kstar0*0.009*0.009);
names[17]="S7_B0Kstar0mumu_0.1-0.98";

central[18]=-0.063;
errors[18]=sqrt(stat2AngB0Kstar0*0.080*0.080+syst2AngB0Kstar0*0.009*0.009);
names[18]="S8_B0Kstar0mumu_0.1-0.98";

central[19]=-0.113;
errors[19]=sqrt(stat2AngB0Kstar0*0.063*0.063+syst2AngB0Kstar0*0.004*0.004);
names[19]="S9_B0Kstar0mumu_0.1-0.98";

central[20]=0.300;
errors[20]=sqrt(stat2AngB0Kstar0*0.171*0.171+syst2AngB0Kstar0*0.023*0.023);
names[20]="P5p_B0Kstar0mumu_0.1-0.98";
}

central[21]=0.326e-7;
errors[21]=sqrt(stat2BRB0Kstar0*0.032e-7*0.032e-7+syst2BRB0Kstar0*0.010e-7*0.010e-7+other2BRB0Kstar0*0.022e-7*0.022e-7);
names[21]="dG_B0Kstar0mumu_1.1-2.5";

central[22]=0.66;
errors[22]=sqrt(stat2AngB0Kstar0*0.083*0.083+syst2AngB0Kstar0*0.022*0.022);
names[22]="FL_B0Kstar0mumu_1.1-2.5";

central[23]=-0.191;
errors[23]=sqrt(stat2AngB0Kstar0*0.080*0.080+syst2AngB0Kstar0*0.012*0.012);
names[23]="AFB_B0Kstar0mumu_1.1-2.5";

central[24]=-0.077;
errors[24]=sqrt(stat2AngB0Kstar0*0.105*0.105+syst2AngB0Kstar0*0.005*0.005);
names[24]="S3_B0Kstar0mumu_1.1-2.5";

central[25]=0.077;
errors[25]=sqrt(stat2AngB0Kstar0*0.113*0.113+syst2AngB0Kstar0*0.005*0.005);
names[25]="S4_B0Kstar0mumu_1.1-2.5";

central[26]=0.137;
errors[26]=sqrt(stat2AngB0Kstar0*0.099*0.099+syst2AngB0Kstar0*0.009*0.009);
names[26]="S5_B0Kstar0mumu_1.1-2.5";

central[27]=-0.219;
errors[27]=sqrt(stat2AngB0Kstar0*0.104*0.104+syst2AngB0Kstar0*0.004*0.004);
names[27]="S7_B0Kstar0mumu_1.1-2.5";

central[28]=0.098;
errors[28]=sqrt(stat2AngB0Kstar0*0.123*0.123+syst2AngB0Kstar0*0.005*0.005);
names[28]="S8_B0Kstar0mumu_1.1-2.5";

central[29]=-0.119;
errors[29]=sqrt(stat2AngB0Kstar0*0.104*0.104+syst2AngB0Kstar0*0.005*0.005);
names[29]="S9_B0Kstar0mumu_1.1-2.5";

central[30]=0.289;
errors[30]=sqrt(stat2AngB0Kstar0*0.220*0.220+syst2AngB0Kstar0*0.023*0.023);
names[30]="P5p_B0Kstar0mumu_1.1-2.5";

central[31]=0.334e-7;
errors[31]=sqrt(stat2BRB0Kstar0*0.033e-7*0.033e-7+syst2BRB0Kstar0*0.009e-7*0.009e-7+other2BRB0Kstar0*0.023e-7*0.023e-7);
names[31]="dG_B0Kstar0mumu_2.5-4";

central[32]=0.876;
errors[32]=sqrt(stat2AngB0Kstar0*0.109*0.109+syst2AngB0Kstar0*0.017*0.017);
names[32]="FL_B0Kstar0mumu_2.5-4";

central[33]=-0.118;
errors[33]=sqrt(stat2AngB0Kstar0*0.090*0.090+syst2AngB0Kstar0*0.007*0.007);
names[33]="AFB_B0Kstar0mumu_2.5-4";

central[34]=0.035;
errors[34]=sqrt(stat2AngB0Kstar0*0.098*0.098+syst2AngB0Kstar0*0.007*0.007);
names[34]="S3_B0Kstar0mumu_2.5-4";

central[35]=0.234;
errors[35]=sqrt(stat2AngB0Kstar0*0.144*0.144+syst2AngB0Kstar0*0.006*0.006);
names[35]="S4_B0Kstar0mumu_2.5-4";

central[36]=-0.022;
errors[36]=sqrt(stat2AngB0Kstar0*0.11*0.11+syst2AngB0Kstar0*0.008*0.008);
names[36]="S5_B0Kstar0mumu_2.5-4";

central[37]=0.068;
errors[37]=sqrt(stat2AngB0Kstar0*0.120*0.120+syst2AngB0Kstar0*0.005*0.005);
names[37]="S7_B0Kstar0mumu_2.5-4";

central[38]=-0.03;
errors[38]=sqrt(stat2AngB0Kstar0*0.131*0.131+syst2AngB0Kstar0*0.006*0.006);
names[38]="S8_B0Kstar0mumu_2.5-4";

central[39]=-0.092;
errors[39]=sqrt(stat2AngB0Kstar0*0.125*0.125+syst2AngB0Kstar0*0.007*0.007);
names[39]="S9_B0Kstar0mumu_2.5-4";

central[40]=-0.066;
errors[40]=sqrt(stat2AngB0Kstar0*0.364*0.364+syst2AngB0Kstar0*0.023*0.023);
names[40]="P5p_B0Kstar0mumu_2.5-4";

central[41]=0.354e-7;
errors[41]=sqrt(stat2BRB0Kstar0*0.027e-7*0.027e-7+syst2BRB0Kstar0*0.009e-7*0.009e-7+other2BRB0Kstar0*0.024e-7*0.024e-7);
names[41]="dG_B0Kstar0mumu_4-6";

central[42]=0.611;
errors[42]=sqrt(stat2AngB0Kstar0*0.053*0.053+syst2AngB0Kstar0*0.017*0.017);
names[42]="FL_B0Kstar0mumu_4-6";

central[43]=0.025;
errors[43]=sqrt(stat2AngB0Kstar0*0.052*0.052+syst2AngB0Kstar0*0.004*0.004);
names[43]="AFB_B0Kstar0mumu_4-6";

central[44]=0.035;
errors[44]=sqrt(stat2AngB0Kstar0*0.069*0.069+syst2AngB0Kstar0*0.007*0.007);
names[44]="S3_B0Kstar0mumu_4-6";

central[45]=0.219;
errors[45]=sqrt(stat2AngB0Kstar0*0.086*0.086+syst2AngB0Kstar0*0.008*0.008);
names[45]="S4_B0Kstar0mumu_4-6";

central[46]=-0.146;
errors[46]=sqrt(stat2AngB0Kstar0*0.078*0.078+syst2AngB0Kstar0*0.011*0.011);
names[46]="S5_B0Kstar0mumu_4-6";

central[47]=-0.016;
errors[47]=sqrt(stat2AngB0Kstar0*0.081*0.081+syst2AngB0Kstar0*0.004*0.004);
names[47]="S7_B0Kstar0mumu_4-6";

central[48]=-0.167;
errors[48]=sqrt(stat2AngB0Kstar0*0.094*0.094+syst2AngB0Kstar0*0.004*0.004);
names[48]="S8_B0Kstar0mumu_4-6";

central[49]=-0.032;
errors[49]=sqrt(stat2AngB0Kstar0*0.071*0.071+syst2AngB0Kstar0*0.004*0.004);
names[49]="S9_B0Kstar0mumu_4-6";

central[50]=-0.3;
errors[50]=sqrt(stat2AngB0Kstar0*0.159*0.159+syst2AngB0Kstar0*0.023*0.023);
names[50]="P5p_B0Kstar0mumu_4-6";

central[51]=0.429e-7;
errors[51]=sqrt(stat2BRB0Kstar0*0.028e-7*0.028e-7+syst2BRB0Kstar0*0.010e-7*0.010e-7+other2BRB0Kstar0*0.029e-7*0.029e-7);
names[51]="dG_B0Kstar0mumu_6-8";

central[52]=0.579;
errors[52]=sqrt(stat2AngB0Kstar0*0.046*0.046+syst2AngB0Kstar0*0.015*0.015);
names[52]="FL_B0Kstar0mumu_6-8";

central[53]=0.152;
errors[53]=sqrt(stat2AngB0Kstar0*0.041*0.041+syst2AngB0Kstar0*0.008*0.008);
names[53]="AFB_B0Kstar0mumu_6-8";

central[54]=-0.042;
errors[54]=sqrt(stat2AngB0Kstar0*0.059*0.059+syst2AngB0Kstar0*0.011*0.011);
names[54]="S3_B0Kstar0mumu_6-8";

central[55]=0.296;
errors[55]=sqrt(stat2AngB0Kstar0*0.067*0.067+syst2AngB0Kstar0*0.011*0.011);
names[55]="S4_B0Kstar0mumu_6-8";

central[56]=-0.249;
errors[56]=sqrt(stat2AngB0Kstar0*0.060*0.060+syst2AngB0Kstar0*0.012*0.012);
names[56]="S5_B0Kstar0mumu_6-8";

central[57]=-0.047;
errors[57]=sqrt(stat2AngB0Kstar0*0.068*0.068+syst2AngB0Kstar0*0.003*0.003);
names[57]="S7_B0Kstar0mumu_6-8";

central[58]=0.085;
errors[58]=sqrt(stat2AngB0Kstar0*0.072*0.072+syst2AngB0Kstar0*0.006*0.006);
names[58]="S8_B0Kstar0mumu_6-8";

central[59]=-0.024;
errors[59]=sqrt(stat2AngB0Kstar0*0.060*0.060+syst2AngB0Kstar0*0.005*0.005);
names[59]="S9_B0Kstar0mumu_6-8";

central[60]=-0.505;
errors[60]=sqrt(stat2AngB0Kstar0*0.122*0.122+syst2AngB0Kstar0*0.024*0.024);
names[60]="P5p_B0Kstar0mumu_6-8";

central[61]=0.487e-7;
errors[61]=sqrt(stat2BRB0Kstar0*0.032e-7*0.032e-7+syst2BRB0Kstar0*0.012e-7*0.012e-7+other2BRB0Kstar0*0.033e-7*0.033e-7);
names[61]="dG_B0Kstar0mumu_11-12.5";

central[62]=0.493;
errors[62]=sqrt(stat2AngB0Kstar0*0.049*0.049+syst2AngB0Kstar0*0.013*0.013);
names[62]="FL_B0Kstar0mumu_11-12.5";

central[63]=0.318;
errors[63]=sqrt(stat2AngB0Kstar0*0.044*0.044+syst2AngB0Kstar0*0.009*0.009);
names[63]="AFB_B0Kstar0mumu_11-12.5";

central[64]=-0.189;
errors[64]=sqrt(stat2AngB0Kstar0*0.058*0.058+syst2AngB0Kstar0*0.005*0.005);
names[64]="S3_B0Kstar0mumu_11-12.5";

central[65]=0.283;
errors[65]=sqrt(stat2AngB0Kstar0*0.095*0.095+syst2AngB0Kstar0*0.009*0.009);
names[65]="S4_B0Kstar0mumu_11-12.5";

central[66]=-0.327;
errors[66]=sqrt(stat2AngB0Kstar0*0.079*0.079+syst2AngB0Kstar0*0.009*0.009);
names[66]="S5_B0Kstar0mumu_11-12.5";

central[67]=-0.141;
errors[67]=sqrt(stat2AngB0Kstar0*0.074*0.074+syst2AngB0Kstar0*0.005*0.005);
names[67]="S7_B0Kstar0mumu_11-12.5";

central[68]=0.007;
errors[68]=sqrt(stat2AngB0Kstar0*0.072*0.072+syst2AngB0Kstar0*0.005*0.005);
names[68]="S8_B0Kstar0mumu_11-12.5";

central[69]=-0.004;
errors[69]=sqrt(stat2AngB0Kstar0*0.073*0.073+syst2AngB0Kstar0*0.006*0.006);
names[69]="S9_B0Kstar0mumu_11-12.5";

central[70]=-0.655;
errors[70]=sqrt(stat2AngB0Kstar0*0.160*0.160+syst2AngB0Kstar0*0.015*0.015);
names[70]="P5p_B0Kstar0mumu_11-12.5";

central[71]=0.534e-7;
errors[71]=sqrt(stat2BRB0Kstar0*0.037e-7*0.037e-7+syst2BRB0Kstar0*0.020e-7*0.020e-7+other2BRB0Kstar0*0.036e-7*0.036e-7);
names[71]="dG_B0Kstar0mumu_15-17";

central[72]=0.349;
errors[72]=sqrt(stat2AngB0Kstar0*0.039*0.039+syst2AngB0Kstar0*0.009*0.009);
names[72]="FL_B0Kstar0mumu_15-17";

central[73]=0.411;
errors[73]=sqrt(stat2AngB0Kstar0*0.041*0.041+syst2AngB0Kstar0*0.008*0.008);
names[73]="AFB_B0Kstar0mumu_15-17";

central[74]=-0.142;
errors[74]=sqrt(stat2AngB0Kstar0*0.049*0.049+syst2AngB0Kstar0*0.007*0.007);
names[74]="S3_B0Kstar0mumu_15-17";

central[75]=0.321;
errors[75]=sqrt(stat2AngB0Kstar0*0.074*0.074+syst2AngB0Kstar0*0.007*0.007);
names[75]="S4_B0Kstar0mumu_15-17";

central[76]=-0.316;
errors[76]=sqrt(stat2AngB0Kstar0*0.057*0.057+syst2AngB0Kstar0*0.009*0.009);
names[76]="S5_B0Kstar0mumu_15-17";

central[77]=0.061;
errors[77]=sqrt(stat2AngB0Kstar0*0.058*0.058+syst2AngB0Kstar0*0.005*0.005);
names[77]="S7_B0Kstar0mumu_15-17";

central[78]=-0.003;
errors[78]=sqrt(stat2AngB0Kstar0*0.061*0.061+syst2AngB0Kstar0*0.003*0.003);
names[78]="S8_B0Kstar0mumu_15-17";

central[79]=-0.019;
errors[79]=sqrt(stat2AngB0Kstar0*0.056*0.056+syst2AngB0Kstar0*0.004*0.004);
names[79]="S9_B0Kstar0mumu_15-17";

central[80]=-0.662;
errors[80]=sqrt(stat2AngB0Kstar0*0.127*0.127+syst2AngB0Kstar0*0.017*0.017);
names[80]="P5p_B0Kstar0mumu_15-17";

central[81]=0.355e-7;
errors[81]=sqrt(stat2BRB0Kstar0*0.027e-7*0.027e-7+syst2BRB0Kstar0*0.017e-7*0.017e-7+other2BRB0Kstar0*0.024e-7*0.024e-7);
names[81]="dG_B0Kstar0mumu_17-19";

central[82]=0.354;
errors[82]=sqrt(stat2AngB0Kstar0*0.049*0.049+syst2AngB0Kstar0*0.025*0.025);
names[82]="FL_B0Kstar0mumu_17-19";

central[83]=0.305;
errors[83]=sqrt(stat2AngB0Kstar0*0.049*0.049+syst2AngB0Kstar0*0.013*0.013);
names[83]="AFB_B0Kstar0mumu_17-19";

central[84]=-0.188;
errors[84]=sqrt(stat2AngB0Kstar0*0.084*0.084+syst2AngB0Kstar0*0.017*0.017);
names[84]="S3_B0Kstar0mumu_17-19";

central[85]=0.266;
errors[85]=sqrt(stat2AngB0Kstar0*0.072*0.072+syst2AngB0Kstar0*0.01*0.01);
names[85]="S4_B0Kstar0mumu_17-19";

central[86]=-0.323;
errors[86]=sqrt(stat2AngB0Kstar0*0.072*0.072+syst2AngB0Kstar0*0.009*0.009);
names[86]="S5_B0Kstar0mumu_17-19";

central[87]=0.044;
errors[87]=sqrt(stat2AngB0Kstar0*0.073*0.073+syst2AngB0Kstar0*0.013*0.013);
names[87]="S7_B0Kstar0mumu_17-19";

central[88]=-0.013;
errors[88]=sqrt(stat2AngB0Kstar0*0.071*0.071+syst2AngB0Kstar0*0.005*0.005);
names[88]="S8_B0Kstar0mumu_17-19";

central[89]=-0.094;
errors[89]=sqrt(stat2AngB0Kstar0*0.067*0.067+syst2AngB0Kstar0*0.004*0.004);
names[89]="S9_B0Kstar0mumu_17-19";

central[90]=-0.676;
errors[90]=sqrt(stat2AngB0Kstar0*0.152*0.152+syst2AngB0Kstar0*0.017*0.017);
names[90]="P5p_B0Kstar0mumu_17-19";

central[91]=0.342e-7;
errors[91]=sqrt(stat2BRB0Kstar0*0.017e-7*0.017e-7+syst2BRB0Kstar0*0.009e-7*0.009e-7+other2BRB0Kstar0*0.023e-7*0.023e-7);
names[91]="dG_B0Kstar0mumu_1.1-6";

central[92]=0.690;
errors[92]=sqrt(stat2AngB0Kstar0*0.036*0.036+syst2AngB0Kstar0*0.017*0.017);
names[92]="FL_B0Kstar0mumu_1.1-6";

central[93]=-0.075;
errors[93]=sqrt(stat2AngB0Kstar0*0.034*0.034+syst2AngB0Kstar0*0.007*0.007);
names[93]="AFB_B0Kstar0mumu_1.1-6";

central[94]=0.012;
errors[94]=sqrt(stat2AngB0Kstar0*0.038*0.038+syst2AngB0Kstar0*0.004*0.004);
names[94]="S3_B0Kstar0mumu_1.1-6";

central[95]=0.155;
errors[95]=sqrt(stat2AngB0Kstar0*0.057*0.057+syst2AngB0Kstar0*0.004*0.004);
names[95]="S4_B0Kstar0mumu_1.1-6";

central[96]=-0.023;
errors[96]=sqrt(stat2AngB0Kstar0*0.050*0.050+syst2AngB0Kstar0*0.005*0.005);
names[96]="S5_B0Kstar0mumu_1.1-6";

central[97]=-0.077;
errors[97]=sqrt(stat2AngB0Kstar0*0.050*0.050+syst2AngB0Kstar0*0.006*0.006);
names[97]="S7_B0Kstar0mumu_1.1-6";

central[98]=-0.028;
errors[98]=sqrt(stat2AngB0Kstar0*0.058*0.058+syst2AngB0Kstar0*0.008*0.008);
names[98]="S8_B0Kstar0mumu_1.1-6";

central[99]=-0.064;
errors[99]=sqrt(stat2AngB0Kstar0*0.042*0.042+syst2AngB0Kstar0*0.004*0.004);
names[99]="S9_B0Kstar0mumu_1.1-6";

central[100]=-0.049;
errors[100]=sqrt(stat2AngB0Kstar0*0.108*0.108+syst2AngB0Kstar0*0.014*0.014);
names[100]="P5p_B0Kstar0mumu_1.1-6";

central[101]=0.436e-7;
errors[101]=sqrt(stat2BRB0Kstar0*0.019e-7*0.019e-7+syst2BRB0Kstar0*0.007e-7*0.007e-7+other2BRB0Kstar0*0.030e-7*0.030e-7);
names[101]="dG_B0Kstar0mumu_15-19";

if(likelihood==1)
{
central[102]=0.344;
errors[102]=sqrt(stat2AngB0Kstar0*0.030*0.030+syst2AngB0Kstar0*0.008*0.008);
names[102]="FL_B0Kstar0mumu_15-19";

central[103]=0.355;
errors[103]=sqrt(stat2AngB0Kstar0*0.027*0.027+syst2AngB0Kstar0*0.009*0.009);
names[103]="AFB_B0Kstar0mumu_15-19";

central[104]=-0.163;
errors[104]=sqrt(stat2AngB0Kstar0*0.033*0.033+syst2AngB0Kstar0*0.009*0.009);
names[104]="S3_B0Kstar0mumu_15-19";

central[105]=0.284;
errors[105]=sqrt(stat2AngB0Kstar0*0.041*0.041+syst2AngB0Kstar0*0.007*0.007);
names[105]="S4_B0Kstar0mumu_15-19";

central[106]=-0.325;
errors[106]=sqrt(stat2AngB0Kstar0*0.037*0.037+syst2AngB0Kstar0*0.009*0.009);
names[106]="S5_B0Kstar0mumu_15-19";

central[107]=0.048;
errors[107]=sqrt(stat2AngB0Kstar0*0.043*0.043+syst2AngB0Kstar0*0.006*0.006);
names[107]="S7_B0Kstar0mumu_15-19";

central[108]=-0.028;
errors[108]=sqrt(stat2AngB0Kstar0*0.045*0.045+syst2AngB0Kstar0*0.003*0.003);
names[108]="S8_B0Kstar0mumu_15-19";

central[109]=-0.053;
errors[109]=sqrt(stat2AngB0Kstar0*0.039*0.039+syst2AngB0Kstar0*0.002*0.002);
names[109]="S9_B0Kstar0mumu_15-19";

central[110]=-0.684;
errors[110]=sqrt(stat2AngB0Kstar0*0.081*0.081+syst2AngB0Kstar0*0.020*0.020);
names[110]="P5p_B0Kstar0mumu_15-19";
}
else
{
central[102]=0.357;
errors[102]=sqrt(stat2AngB0Kstar0*0.035*0.035+syst2AngB0Kstar0*0.011*0.011);
names[102]="FL_B0Kstar0mumu_15-19";

central[103]=0.367;
errors[103]=sqrt(stat2AngB0Kstar0*0.037*0.037+syst2AngB0Kstar0*0.007*0.007);
names[103]="AFB_B0Kstar0mumu_15-19";

central[104]=-0.135;
errors[104]=sqrt(stat2AngB0Kstar0*0.050*0.050+syst2AngB0Kstar0*0.012*0.012);
names[104]="S3_B0Kstar0mumu_15-19";

central[105]=0.314;
errors[105]=sqrt(stat2AngB0Kstar0*0.054*0.054+syst2AngB0Kstar0*0.027*0.027);
names[105]="S4_B0Kstar0mumu_15-19";

central[106]=-0.335;
errors[106]=sqrt(stat2AngB0Kstar0*0.047*0.047+syst2AngB0Kstar0*0.007*0.007);
names[106]="S5_B0Kstar0mumu_15-19";

central[107]=0.066;
errors[107]=sqrt(stat2AngB0Kstar0*0.049*0.049+syst2AngB0Kstar0*0.014*0.014);
names[107]="S7_B0Kstar0mumu_15-19";

central[108]=-0.024;
errors[108]=sqrt(stat2AngB0Kstar0*0.048*0.048+syst2AngB0Kstar0*0.009*0.009);
names[108]="S8_B0Kstar0mumu_15-19";

central[109]=-0.056;
errors[109]=sqrt(stat2AngB0Kstar0*0.047*0.047+syst2AngB0Kstar0*0.014*0.014);
names[109]="S9_B0Kstar0mumu_15-19";

central[110]=-0.709;
errors[110]=sqrt(stat2AngB0Kstar0*0.093*0.093+syst2AngB0Kstar0*0.016*0.016);
names[110]="P5p_B0Kstar0mumu_15-19";
}

central[111]=0.;
errors[111]=0.;
names[111]="dG_B0Kstar0mumu_1.1-2";

central[112]=0.768;
errors[112]=sqrt(stat2AngB0Kstar0*0.141*0.141+syst2AngB0Kstar0*0.025*0.025);
names[112]="FL_B0Kstar0mumu_1.1-2";

central[113]=-0.333;
errors[113]=sqrt(stat2AngB0Kstar0*0.130*0.130+syst2AngB0Kstar0*0.012*0.012);
names[113]="AFB_B0Kstar0mumu_1.1-2";

central[114]=0.065;
errors[114]=sqrt(stat2AngB0Kstar0*0.137*0.137+syst2AngB0Kstar0*0.007*0.007);
names[114]="S3_B0Kstar0mumu_1.1-2";

central[115]=-0.127;
errors[115]=sqrt(stat2AngB0Kstar0*0.190*0.190+syst2AngB0Kstar0*0.027*0.027);
names[115]="S4_B0Kstar0mumu_1.1-2";

central[116]=0.286;
errors[116]=sqrt(stat2AngB0Kstar0*0.172*0.172+syst2AngB0Kstar0*0.009*0.009);
names[116]="S5_B0Kstar0mumu_1.1-2";

central[117]=-0.293;
errors[117]=sqrt(stat2AngB0Kstar0*0.180*0.180+syst2AngB0Kstar0*0.005*0.005);
names[117]="S7_B0Kstar0mumu_1.1-2";

central[118]=0.114;
errors[118]=sqrt(stat2AngB0Kstar0*0.196*0.196+syst2AngB0Kstar0*0.006*0.006);
names[118]="S8_B0Kstar0mumu_1.1-2";

central[119]=-0.110;
errors[119]=sqrt(stat2AngB0Kstar0*0.140*0.140+syst2AngB0Kstar0*0.001*0.001);
names[119]="S9_B0Kstar0mumu_1.1-2";

central[120]=0.606;
errors[120]=sqrt(stat2AngB0Kstar0*0.769*0.769+syst2AngB0Kstar0*0.017*0.017);
names[120]="P5p_B0Kstar0mumu_1.1-2";

central[121]=0.;
errors[121]=0.;
names[121]="dG_B0Kstar0mumu_2-3";

central[122]=0.690;
errors[122]=sqrt(stat2AngB0Kstar0*0.113*0.113+syst2AngB0Kstar0*0.023*0.023);
names[122]="FL_B0Kstar0mumu_2-3";

central[123]=-0.158;
errors[123]=sqrt(stat2AngB0Kstar0*0.090*0.090+syst2AngB0Kstar0*0.008*0.008);
names[123]="AFB_B0Kstar0mumu_2-3";

central[124]=0.006;
errors[124]=sqrt(stat2AngB0Kstar0*0.100*0.100+syst2AngB0Kstar0*0.007*0.007);
names[124]="S3_B0Kstar0mumu_2-3";

central[125]=0.339;
errors[125]=sqrt(stat2AngB0Kstar0*0.140*0.140+syst2AngB0Kstar0*0.041*0.041);
names[125]="S4_B0Kstar0mumu_2-3";

central[126]=0.206;
errors[126]=sqrt(stat2AngB0Kstar0*0.131*0.131+syst2AngB0Kstar0*0.009*0.009);
names[126]="S5_B0Kstar0mumu_2-3";

central[127]=-0.252;
errors[127]=sqrt(stat2AngB0Kstar0*0.151*0.151+syst2AngB0Kstar0*0.002*0.002);
names[127]="S7_B0Kstar0mumu_2-3";

central[128]=0.176;
errors[128]=sqrt(stat2AngB0Kstar0*0.165*0.165+syst2AngB0Kstar0*0.006*0.006);
names[128]="S8_B0Kstar0mumu_2-3";

central[129]=0.000;
errors[129]=sqrt(stat2AngB0Kstar0*0.102*0.102+syst2AngB0Kstar0*0.003*0.003);
names[129]="S9_B0Kstar0mumu_2-3";

central[130]=0.461;
errors[130]=sqrt(stat2AngB0Kstar0*0.313*0.313+syst2AngB0Kstar0*0.019*0.019);
names[130]="P5p_B0Kstar0mumu_2-3";

central[131]=0.;
errors[131]=0.;
names[131]="dG_B0Kstar0mumu_3-4";

central[132]=0.873;
errors[132]=sqrt(stat2AngB0Kstar0*0.154*0.154+syst2AngB0Kstar0*0.023*0.023);
names[132]="FL_B0Kstar0mumu_3-4";

central[133]=-0.041;
errors[133]=sqrt(stat2AngB0Kstar0*0.091*0.091+syst2AngB0Kstar0*0.002*0.002);
names[133]="AFB_B0Kstar0mumu_3-4";

central[134]=0.078;
errors[134]=sqrt(stat2AngB0Kstar0*0.131*0.131+syst2AngB0Kstar0*0.008*0.008);
names[134]="S3_B0Kstar0mumu_3-4";

central[135]=0.046;
errors[135]=sqrt(stat2AngB0Kstar0*0.196*0.196+syst2AngB0Kstar0*0.046*0.046);
names[135]="S4_B0Kstar0mumu_3-4";

central[136]=-0.110;
errors[136]=sqrt(stat2AngB0Kstar0*0.169*0.169+syst2AngB0Kstar0*0.004*0.004);
names[136]="S5_B0Kstar0mumu_3-4";

central[137]=0.171;
errors[137]=sqrt(stat2AngB0Kstar0*0.175*0.175+syst2AngB0Kstar0*0.002*0.002);
names[137]="S7_B0Kstar0mumu_3-4";

central[138]=-0.097;
errors[138]=sqrt(stat2AngB0Kstar0*0.189*0.189+syst2AngB0Kstar0*0.002*0.002);
names[138]="S8_B0Kstar0mumu_3-4";

central[139]=-0.203;
errors[139]=sqrt(stat2AngB0Kstar0*0.132*0.132+syst2AngB0Kstar0*0.002*0.002);
names[139]="S9_B0Kstar0mumu_3-4";

central[140]=-0.295;
errors[140]=sqrt(stat2AngB0Kstar0*7.112*7.112+syst2AngB0Kstar0*0.023*0.023);
names[140]="P5p_B0Kstar0mumu_3-4";

central[141]=0.;
errors[141]=0.;
names[141]="dG_B0Kstar0mumu_4-5";

central[142]=0.899;
errors[142]=sqrt(stat2AngB0Kstar0*0.106*0.106+syst2AngB0Kstar0*0.023*0.023);
names[142]="FL_B0Kstar0mumu_4-5";

central[143]=0.052;
errors[143]=sqrt(stat2AngB0Kstar0*0.080*0.080+syst2AngB0Kstar0*0.004*0.004);
names[143]="AFB_B0Kstar0mumu_4-5";

central[144]=0.200;
errors[144]=sqrt(stat2AngB0Kstar0*0.101*0.101+syst2AngB0Kstar0*0.007*0.007);
names[144]="S3_B0Kstar0mumu_4-5";

central[145]=0.148;
errors[145]=sqrt(stat2AngB0Kstar0*0.154*0.154+syst2AngB0Kstar0*0.047*0.047);
names[145]="S4_B0Kstar0mumu_4-5";

central[146]=-0.306;
errors[146]=sqrt(stat2AngB0Kstar0*0.141*0.141+syst2AngB0Kstar0*0.004*0.004);
names[146]="S5_B0Kstar0mumu_4-5";

central[147]=-0.082;
errors[147]=sqrt(stat2AngB0Kstar0*0.129*0.129+syst2AngB0Kstar0*0.001*0.001);
names[147]="S7_B0Kstar0mumu_4-5";

central[148]=-0.107;
errors[148]=sqrt(stat2AngB0Kstar0*0.146*0.146+syst2AngB0Kstar0*0.003*0.003);
names[148]="S8_B0Kstar0mumu_4-5";

central[149]=0.181;
errors[149]=sqrt(stat2AngB0Kstar0*0.105*0.105+syst2AngB0Kstar0*0.001*0.001);
names[149]="S9_B0Kstar0mumu_4-5";

central[150]=-0.799;
errors[150]=sqrt(stat2AngB0Kstar0*18.19*18.19+syst2AngB0Kstar0*0.022*0.022);
names[150]="P5p_B0Kstar0mumu_4-5";

central[151]=0.;
errors[151]=0.;
names[151]="dG_B0Kstar0mumu_5-6";

central[152]=0.644;
errors[152]=sqrt(stat2AngB0Kstar0*0.130*0.130+syst2AngB0Kstar0*0.025*0.025);
names[152]="FL_B0Kstar0mumu_5-6";

central[153]=0.057;
errors[153]=sqrt(stat2AngB0Kstar0*0.094*0.094+syst2AngB0Kstar0*0.006*0.006);
names[153]="AFB_B0Kstar0mumu_5-6";

central[154]=-0.122;
errors[154]=sqrt(stat2AngB0Kstar0*0.126*0.126+syst2AngB0Kstar0*0.009*0.009);
names[154]="S3_B0Kstar0mumu_5-6";

central[155]=0.273;
errors[155]=sqrt(stat2AngB0Kstar0*0.184*0.184+syst2AngB0Kstar0*0.048*0.048);
names[155]="S4_B0Kstar0mumu_5-6";

central[156]=-0.095;
errors[156]=sqrt(stat2AngB0Kstar0*0.142*0.142+syst2AngB0Kstar0*0.004*0.004);
names[156]="S5_B0Kstar0mumu_5-6";

central[157]=0.038;
errors[157]=sqrt(stat2AngB0Kstar0*0.135*0.135+syst2AngB0Kstar0*0.002*0.002);
names[157]="S7_B0Kstar0mumu_5-6";

central[158]=0.037;
errors[158]=sqrt(stat2AngB0Kstar0*0.160*0.160+syst2AngB0Kstar0*0.003*0.003);
names[158]="S8_B0Kstar0mumu_5-6";

central[159]=-0.080;
errors[159]=sqrt(stat2AngB0Kstar0*0.120*0.120+syst2AngB0Kstar0*0.001*0.001);
names[159]="S9_B0Kstar0mumu_5-6";

central[160]=-0.197;
errors[160]=sqrt(stat2AngB0Kstar0*0.334*0.334+syst2AngB0Kstar0*0.018*0.018);
names[160]="P5p_B0Kstar0mumu_5-6";

central[161]=0.;
errors[161]=0.;
names[161]="dG_B0Kstar0mumu_6-7";

central[162]=0.644;
errors[162]=sqrt(stat2AngB0Kstar0*0.089*0.089+syst2AngB0Kstar0*0.025*0.025);
names[162]="FL_B0Kstar0mumu_6-7";

central[163]=0.058;
errors[163]=sqrt(stat2AngB0Kstar0*0.064*0.064+syst2AngB0Kstar0*0.009*0.009);
names[163]="AFB_B0Kstar0mumu_6-7";

central[164]=-0.069;
errors[164]=sqrt(stat2AngB0Kstar0*0.091*0.091+syst2AngB0Kstar0*0.004*0.004);
names[164]="S3_B0Kstar0mumu_6-7";

central[165]=0.311;
errors[165]=sqrt(stat2AngB0Kstar0*0.118*0.118+syst2AngB0Kstar0*0.052*0.052);
names[165]="S4_B0Kstar0mumu_6-7";

central[166]=-0.339;
errors[166]=sqrt(stat2AngB0Kstar0*0.114*0.114+syst2AngB0Kstar0*0.008*0.008);
names[166]="S5_B0Kstar0mumu_6-7";

central[167]=0.009;
errors[167]=sqrt(stat2AngB0Kstar0*0.124*0.124+syst2AngB0Kstar0*0.004*0.004);
names[167]="S7_B0Kstar0mumu_6-7";

central[168]=-0.080;
errors[168]=sqrt(stat2AngB0Kstar0*0.131*0.131+syst2AngB0Kstar0*0.002*0.002);
names[168]="S8_B0Kstar0mumu_6-7";

central[169]=0.061;
errors[169]=sqrt(stat2AngB0Kstar0*0.091*0.091+syst2AngB0Kstar0*0.001*0.001);
names[169]="S9_B0Kstar0mumu_6-7";

central[170]=-0.713;
errors[170]=sqrt(stat2AngB0Kstar0*0.268*0.268+syst2AngB0Kstar0*0.015*0.015);
names[170]="P5p_B0Kstar0mumu_6-7";

central[171]=0.;
errors[171]=0.;
names[171]="dG_B0Kstar0mumu_7-8";

central[172]=0.609;
errors[172]=sqrt(stat2AngB0Kstar0*0.103*0.103+syst2AngB0Kstar0*0.025*0.025);
names[172]="FL_B0Kstar0mumu_7-8";

central[173]=0.241;
errors[173]=sqrt(stat2AngB0Kstar0*0.080*0.080+syst2AngB0Kstar0*0.0012*0.012);
names[173]="AFB_B0Kstar0mumu_7-8";

central[174]=-0.054;
errors[174]=sqrt(stat2AngB0Kstar0*0.099*0.099+syst2AngB0Kstar0*0.005*0.005);
names[174]="S3_B0Kstar0mumu_7-8";

central[175]=0.236;
errors[175]=sqrt(stat2AngB0Kstar0*0.136*0.136+syst2AngB0Kstar0*0.058*0.058);
names[175]="S4_B0Kstar0mumu_7-8";

central[176]=-0.386;
errors[176]=sqrt(stat2AngB0Kstar0*0.135*0.135+syst2AngB0Kstar0*0.007*0.007);
names[176]="S5_B0Kstar0mumu_7-8";

central[177]=-0.094;
errors[177]=sqrt(stat2AngB0Kstar0*0.130*0.130+syst2AngB0Kstar0*0.003*0.003);
names[177]="S7_B0Kstar0mumu_7-8";

central[178]=0.295;
errors[178]=sqrt(stat2AngB0Kstar0*0.139*0.139+syst2AngB0Kstar0*0.002*0.002);
names[178]="S8_B0Kstar0mumu_7-8";

central[179]=0.030;
errors[179]=sqrt(stat2AngB0Kstar0*0.100*0.100+syst2AngB0Kstar0*0.001*0.001);
names[179]="S9_B0Kstar0mumu_7-8";

central[180]=-0.808;
errors[180]=sqrt(stat2AngB0Kstar0*0.303*0.303+syst2AngB0Kstar0*0.010*0.010);
names[180]="P5p_B0Kstar0mumu_7-8";

central[181]=0.;
errors[181]=0.;
names[181]="dG_B0Kstar0mumu_11-11.75";

central[182]=0.502;
errors[182]=sqrt(stat2AngB0Kstar0*0.090*0.090+syst2AngB0Kstar0*0.022*0.022);
names[182]="FL_B0Kstar0mumu_11-11.75";

central[183]=0.370;
errors[183]=sqrt(stat2AngB0Kstar0*0.076*0.076+syst2AngB0Kstar0*0.015*0.015);
names[183]="AFB_B0Kstar0mumu_11-11.75";

central[184]=-0.217;
errors[184]=sqrt(stat2AngB0Kstar0*0.090*0.090+syst2AngB0Kstar0*0.008*0.008);
names[184]="S3_B0Kstar0mumu_11-11.75";

central[185]=0.252;
errors[185]=sqrt(stat2AngB0Kstar0*0.113*0.113+syst2AngB0Kstar0*0.063*0.063);
names[185]="S4_B0Kstar0mumu_11-11.75";

central[186]=-0.235;
errors[186]=sqrt(stat2AngB0Kstar0*0.115*0.115+syst2AngB0Kstar0*0.013*0.013);
names[186]="S5_B0Kstar0mumu_11-11.75";

central[187]=-0.110;
errors[187]=sqrt(stat2AngB0Kstar0*0.114*0.114+syst2AngB0Kstar0*0.002*0.002);
names[187]="S7_B0Kstar0mumu_11-11.75";

central[188]=0.079;
errors[188]=sqrt(stat2AngB0Kstar0*0.122*0.122+syst2AngB0Kstar0*0.003*0.003);
names[188]="S8_B0Kstar0mumu_11-11.75";

central[189]=-0.084;
errors[189]=sqrt(stat2AngB0Kstar0*0.102*0.102+syst2AngB0Kstar0*0.003*0.003);
names[189]="S9_B0Kstar0mumu_11-11.75";

central[190]=-0.485;
errors[190]=sqrt(stat2AngB0Kstar0*0.224*0.224+syst2AngB0Kstar0*0.028*0.028);
names[190]="P5p_B0Kstar0mumu_11-11.75";

central[191]=0.;
errors[191]=0.;
names[191]="dG_B0Kstar0mumu_11.75-12.5";

central[192]=0.734;
errors[192]=sqrt(stat2AngB0Kstar0*0.107*0.107+syst2AngB0Kstar0*0.018*0.018);
names[192]="FL_B0Kstar0mumu_11.75-12.5";

central[193]=0.293;
errors[193]=sqrt(stat2AngB0Kstar0*0.064*0.064+syst2AngB0Kstar0*0.014*0.014);
names[193]="AFB_B0Kstar0mumu_11.75-12.5";

central[194]=-0.157;
errors[194]=sqrt(stat2AngB0Kstar0*0.098*0.098+syst2AngB0Kstar0*0.008*0.008);
names[194]="S3_B0Kstar0mumu_11.75-12.5";

central[195]=0.309;
errors[195]=sqrt(stat2AngB0Kstar0*0.111*0.111+syst2AngB0Kstar0*0.056*0.056);
names[195]="S4_B0Kstar0mumu_11.75-12.5";

central[196]=-0.366;
errors[196]=sqrt(stat2AngB0Kstar0*0.112*0.112+syst2AngB0Kstar0*0.012*0.012);
names[196]="S5_B0Kstar0mumu_11.75-12.5";

central[197]=-0.212;
errors[197]=sqrt(stat2AngB0Kstar0*0.118*0.118+syst2AngB0Kstar0*0.002*0.002);
names[197]="S7_B0Kstar0mumu_11.75-12.5";

central[198]=0.090;
errors[198]=sqrt(stat2AngB0Kstar0*0.111*0.111+syst2AngB0Kstar0*0.003*0.003);
names[198]="S8_B0Kstar0mumu_11.75-12.5";

central[199]=0.030;
errors[199]=sqrt(stat2AngB0Kstar0*0.093*0.093+syst2AngB0Kstar0*0.002*0.002);
names[199]="S9_B0Kstar0mumu_11.75-12.5";

central[200]=-0.827;
errors[200]=sqrt(stat2AngB0Kstar0*0.357*0.357+syst2AngB0Kstar0*0.026*0.026);
names[200]="P5p_B0Kstar0mumu_11.75-12.5";

central[201]=0.;
errors[201]=0.;
names[201]="dG_B0Kstar0mumu_15-16";

central[202]=0.385;
errors[202]=sqrt(stat2AngB0Kstar0*0.067*0.067+syst2AngB0Kstar0*0.013*0.013);
names[202]="FL_B0Kstar0mumu_15-16";

central[203]=0.396;
errors[203]=sqrt(stat2AngB0Kstar0*0.068*0.068+syst2AngB0Kstar0*0.009*0.009);
names[203]="AFB_B0Kstar0mumu_15-16";

central[204]=-0.060;
errors[204]=sqrt(stat2AngB0Kstar0*0.088*0.088+syst2AngB0Kstar0*0.006*0.006);
names[204]="S3_B0Kstar0mumu_15-16";

central[205]=0.321;
errors[205]=sqrt(stat2AngB0Kstar0*0.099*0.099+syst2AngB0Kstar0*0.007*0.007);
names[205]="S4_B0Kstar0mumu_15-16";

central[206]=-0.360;
errors[206]=sqrt(stat2AngB0Kstar0*0.092*0.092+syst2AngB0Kstar0*0.006*0.006);
names[206]="S5_B0Kstar0mumu_15-16";

central[207]=0.040;
errors[207]=sqrt(stat2AngB0Kstar0*0.092*0.092+syst2AngB0Kstar0*0.002*0.002);
names[207]="S7_B0Kstar0mumu_15-16";

central[208]=0.057;
errors[208]=sqrt(stat2AngB0Kstar0*0.095*0.095+syst2AngB0Kstar0*0.005*0.005);
names[208]="S8_B0Kstar0mumu_15-16";

central[209]=-0.054;
errors[209]=sqrt(stat2AngB0Kstar0*0.087*0.087+syst2AngB0Kstar0*0.005*0.005);
names[209]="S9_B0Kstar0mumu_15-16";

central[210]=-0.758;
errors[210]=sqrt(stat2AngB0Kstar0*0.179*0.179+syst2AngB0Kstar0*0.013*0.013);
names[210]="P5p_B0Kstar0mumu_15-16";

central[211]=0.;
errors[211]=0.;
names[211]="dG_B0Kstar0mumu_16-17";

central[212]=0.295;
errors[212]=sqrt(stat2AngB0Kstar0*0.062*0.062+syst2AngB0Kstar0*0.013*0.013);
names[212]="FL_B0Kstar0mumu_16-17";

central[213]=0.451;
errors[213]=sqrt(stat2AngB0Kstar0*0.071*0.071+syst2AngB0Kstar0*0.007*0.007);
names[213]="AFB_B0Kstar0mumu_16-17";

central[214]=-0.250;
errors[214]=sqrt(stat2AngB0Kstar0*0.092*0.092+syst2AngB0Kstar0*0.007*0.007);
names[214]="S3_B0Kstar0mumu_16-17";

central[215]=0.246;
errors[215]=sqrt(stat2AngB0Kstar0*0.096*0.096+syst2AngB0Kstar0*0.029*0.029);
names[215]="S4_B0Kstar0mumu_16-17";

central[216]=-0.254;
errors[216]=sqrt(stat2AngB0Kstar0*0.081*0.081+syst2AngB0Kstar0*0.010*0.010);
names[216]="S5_B0Kstar0mumu_16-17";

central[217]=0.144;
errors[217]=sqrt(stat2AngB0Kstar0*0.091*0.091+syst2AngB0Kstar0*0.005*0.005);
names[217]="S7_B0Kstar0mumu_16-17";

central[218]=-0.055;
errors[218]=sqrt(stat2AngB0Kstar0*0.090*0.090+syst2AngB0Kstar0*0.005*0.005);
names[218]="S8_B0Kstar0mumu_16-17";

central[219]=-0.014;
errors[219]=sqrt(stat2AngB0Kstar0*0.086*0.086+syst2AngB0Kstar0*0.004*0.004);
names[219]="S9_B0Kstar0mumu_16-17";

central[220]=-0.567;
errors[220]=sqrt(stat2AngB0Kstar0*0.186*0.186+syst2AngB0Kstar0*0.014*0.014);
names[220]="P5p_B0Kstar0mumu_16-17";

central[221]=0.;
errors[221]=0.;
names[221]="dG_B0Kstar0mumu_17-18";

central[222]=0.363;
errors[222]=sqrt(stat2AngB0Kstar0*0.073*0.073+syst2AngB0Kstar0*0.017*0.017);
names[222]="FL_B0Kstar0mumu_17-18";

central[223]=0.274;
errors[223]=sqrt(stat2AngB0Kstar0*0.069*0.069+syst2AngB0Kstar0*0.008*0.008);
names[223]="AFB_B0Kstar0mumu_17-18";

central[224]=-0.099;
errors[224]=sqrt(stat2AngB0Kstar0*0.092*0.092+syst2AngB0Kstar0*0.011*0.011);
names[224]="S3_B0Kstar0mumu_17-18";

central[225]=0.229;
errors[225]=sqrt(stat2AngB0Kstar0*0.096*0.096+syst2AngB0Kstar0*0.045*0.045);
names[225]="S4_B0Kstar0mumu_17-18";

central[226]=-0.305;
errors[226]=sqrt(stat2AngB0Kstar0*0.088*0.088+syst2AngB0Kstar0*0.015*0.015);
names[226]="S5_B0Kstar0mumu_17-18";

central[227]=0.022;
errors[227]=sqrt(stat2AngB0Kstar0*0.094*0.094+syst2AngB0Kstar0*0.011*0.011);
names[227]="S7_B0Kstar0mumu_17-18";

central[228]=0.007;
errors[228]=sqrt(stat2AngB0Kstar0*0.098*0.098+syst2AngB0Kstar0*0.001*0.001);
names[228]="S8_B0Kstar0mumu_17-18";

central[229]=-0.090;
errors[229]=sqrt(stat2AngB0Kstar0*0.095*0.095+syst2AngB0Kstar0*0.002*0.002);
names[229]="S9_B0Kstar0mumu_17-18";

central[230]=-0.646;
errors[230]=sqrt(stat2AngB0Kstar0*0.190*0.190+syst2AngB0Kstar0*0.027*0.027);
names[230]="P5p_B0Kstar0mumu_17-18";

central[231]=0.;
errors[231]=0.;
names[231]="dG_B0Kstar0mumu_18-19";

central[232]=0.421;
errors[232]=sqrt(stat2AngB0Kstar0*0.100*0.100+syst2AngB0Kstar0*0.013*0.013);
names[232]="FL_B0Kstar0mumu_18-19";

central[233]=0.354;
errors[233]=sqrt(stat2AngB0Kstar0*0.111*0.111+syst2AngB0Kstar0*0.012*0.012);
names[233]="AFB_B0Kstar0mumu_18-19";

central[234]=-0.131;
errors[234]=sqrt(stat2AngB0Kstar0*0.130*0.130+syst2AngB0Kstar0*0.012*0.012);
names[234]="S3_B0Kstar0mumu_18-19";

central[235]=0.607;
errors[235]=sqrt(stat2AngB0Kstar0*0.170*0.170+syst2AngB0Kstar0*0.059*0.059);
names[235]="S4_B0Kstar0mumu_18-19";

central[236]=-0.534;
errors[236]=sqrt(stat2AngB0Kstar0*0.150*0.150+syst2AngB0Kstar0*0.015*0.015);
names[236]="S5_B0Kstar0mumu_18-19";

central[237]=0.058;
errors[237]=sqrt(stat2AngB0Kstar0*0.124*0.124+syst2AngB0Kstar0*0.006*0.006);
names[237]="S7_B0Kstar0mumu_18-19";

central[238]=-0.149;
errors[238]=sqrt(stat2AngB0Kstar0*0.139*0.139+syst2AngB0Kstar0*0.010*0.010);
names[238]="S8_B0Kstar0mumu_18-19";

central[239]=-0.079;
errors[239]=sqrt(stat2AngB0Kstar0*0.122*0.122+syst2AngB0Kstar0*0.007*0.007);
names[239]="S9_B0Kstar0mumu_18-19";

central[240]=-1.070;
errors[240]=sqrt(stat2AngB0Kstar0*0.349*0.349+syst2AngB0Kstar0*0.029*0.029);
names[240]="P5p_B0Kstar0mumu_18-19";

central[241]=59.2e-9;
errors[241]=sqrt(stat2BRBpKstarp*14.4*14.4+syst2BRBpKstarp*4.0*4.0)*1.e-9;
names[241]="dG_B+Kstar+mumu_0.1-2";

central[242]=55.9e-9;
errors[242]=sqrt(stat2BRBpKstarp*15.9*15.9+syst2BRBpKstarp*3.8*3.8)*1.e-9;
names[242]="dG_B+Kstar+mumu_2-4";

central[243]=24.9e-9;
errors[243]=sqrt(stat2BRBpKstarp*11.0*11.0+syst2BRBpKstarp*1.7*1.7)*1.e-9;
names[243]="dG_B+Kstar+mumu_4-6";

central[244]=33.0e-9;
errors[244]=sqrt(stat2BRBpKstarp*11.3*11.3+syst2BRBpKstarp*2.3*2.3)*1.e-9;
names[244]="dG_B+Kstar+mumu_6-8";

central[245]=82.8e-9;
errors[245]=sqrt(stat2BRBpKstarp*15.8*15.8+syst2BRBpKstarp*5.6*5.6)*1.e-9;
names[245]="dG_B+Kstar+mumu_11-12.5";

central[246]=64.4e-9;
errors[246]=sqrt(stat2BRBpKstarp*12.9*12.9+syst2BRBpKstarp*4.4*4.4)*1.e-9;
names[246]="dG_B+Kstar+mumu_15-17";

central[247]=11.6e-9;
errors[247]=sqrt(stat2BRBpKstarp*9.1*9.1+syst2BRBpKstarp*0.8*0.8)*1.e-9;
names[247]="dG_B+Kstar+mumu_17-19";

central[248]=36.6e-9;
errors[248]=sqrt(stat2BRBpKstarp*8.3*8.3+syst2BRBpKstarp*2.6*2.6)*1.e-9;
names[248]="dG_B+Kstar+mumu_1.1-6";

central[249]=39.5e-9;
errors[249]=sqrt(stat2BRBpKstarp*8.0*8.0+syst2BRBpKstarp*2.8*2.8)*1.e-9;
names[249]="dG_B+Kstar+mumu_15-19";

central[250]=1.03e-6;
errors[250]=0.19e-6;
names[250]="BR_B0Kstar0ee_full";

central[251]=3.1e-7;
errors[251]=sqrt(stat2BRBKstaree*0.9*0.9+syst2BRBKstaree*0.3*0.3+other2BRBKstaree*0.2*0.2)*1.e-7;
names[251]="dG_B0Kstar0ee_0.0009-1";

central[252]=0.;
errors[252]=0.;
names[252]="dG_B0Kstar0ee_0.002-1.12";

central[253]=0.16;
errors[253]=sqrt(stat2AngBKstaree*0.06*0.06+syst2AngBKstaree*0.03*0.03);
names[253]="FL_B0Kstar0ee_0.002_1.12";

central[254]=-0.23;
errors[254]=sqrt(stat2AngBKstaree*0.23*0.23+syst2AngBKstaree*0.05*0.05);
names[254]="AT2_B0Kstar0ee_0.002_1.12";

central[255]=0.10;
errors[255]=sqrt(stat2AngBKstaree*0.18*0.18+syst2AngBKstaree*0.05*0.05);
names[255]="ATRe_B0Kstar0ee_0.002_1.12";

central[256]=0.;
errors[256]=0.;
names[256]="dG_B0Kstar0ee_1.1-6";

central[257]=0.;
errors[257]=0.;
names[257]="FL_B0Kstar0ee_1.1_6";

central[258]=0.;
errors[258]=0.;
names[258]="P1_B0Kstar0ee_1.1_6";

central[259]=0.;
errors[259]=0.;
names[259]="ATRe_B0Kstar0ee_1.1_6";

central[260]=12.2e-9;
errors[260]=sqrt(stat2BRBK*5.9*5.9+syst2BRBK*0.6*0.6)*1.e-9;
names[260]="dG_B0K0mumu_0.1-2";

central[261]=18.7e-9;
errors[261]=sqrt(stat2BRBK*5.5*5.5+syst2BRBK*0.9*0.9)*1.e-9;
names[261]="dG_B0K0mumu_2-4";

central[262]=17.3e-9;
errors[262]=sqrt(stat2BRBK*5.3*5.3+syst2BRBK*0.9*0.9)*1.e-9;
names[262]="dG_B0K0mumu_4-6";

central[263]=27.0e-9;
errors[263]=sqrt(stat2BRBK*5.8*5.8+syst2BRBK*1.4*1.4)*1.e-9;
names[263]="dG_B0K0mumu_6-8";

central[264]=12.7e-9;
errors[264]=sqrt(stat2BRBK*4.5*4.5+syst2BRBK*0.6*0.6)*1.e-9;
names[264]="dG_B0K0mumu_11-12.5";

central[265]=14.3e-9;
errors[265]=sqrt(stat2BRBK*3.5*3.5+syst2BRBK*0.7*0.7)*1.e-9;
names[265]="dG_B0K0mumu_15-17";

central[266]=7.8e-9;
errors[266]=sqrt(stat2BRBK*1.7*1.7+syst2BRBK*0.4*0.4)*1.e-9;
names[266]="dG_B0K0mumu_17-22";

central[267]=18.7e-9;
errors[267]=sqrt(stat2BRBK*3.5*3.5+syst2BRBK*0.9*0.9)*1.e-9;
names[267]="dG_B0K0mumu_1.1-6";

central[268]=9.5e-9;
errors[268]=sqrt(stat2BRBK*1.6*1.6+syst2BRBK*0.5*0.5)*1.e-9;
names[268]="dG_B0K0mumu_15-22";

central[269]=33.2e-9;
errors[269]=sqrt(stat2BRBpKp*1.8*1.8+syst2BRBpKp*1.7*1.7)*1.e-9;
names[269]="dG_B+K+mumu_0.1-0.98";

central[270]=23.3e-9;
errors[270]=sqrt(stat2BRBpKp*1.5*1.5+syst2BRBpKp*1.2*1.2)*1.e-9;
names[270]="dG_B+K+mumu_1.1-2";

central[271]=28.2e-9;
errors[271]=sqrt(stat2BRBpKp*1.6*1.6+syst2BRBpKp*1.4*1.4)*1.e-9;
names[271]="dG_B+K+mumu_2-3";

central[272]=25.4e-9;
errors[272]=sqrt(stat2BRBpKp*1.5*1.5+syst2BRBpKp*1.3*1.3)*1.e-9;
names[272]="dG_B+K+mumu_3-4";

central[273]=22.1e-9;
errors[273]=sqrt(stat2BRBpKp*1.4*1.4+syst2BRBpKp*1.1*1.1)*1.e-9;
names[273]="dG_B+K+mumu_4-5";

central[274]=23.1e-9;
errors[274]=sqrt(stat2BRBpKp*1.4*1.4+syst2BRBpKp*1.2*1.2)*1.e-9;
names[274]="dG_B+K+mumu_5-6";

central[275]=24.5e-9;
errors[275]=sqrt(stat2BRBpKp*1.4*1.4+syst2BRBpKp*1.2*1.2)*1.e-9;
names[275]="dG_B+K+mumu_6-7";

central[276]=23.1e-9;
errors[276]=sqrt(stat2BRBpKp*1.4*1.4+syst2BRBpKp*1.2*1.2)*1.e-9;
names[276]="dG_B+K+mumu_7-8";

central[277]=17.7e-9;
errors[277]=sqrt(stat2BRBpKp*1.3*1.3+syst2BRBpKp*0.9*0.9)*1.e-9;
names[277]="dG_B+K+mumu_11-11.8";

central[278]=19.3e-9;
errors[278]=sqrt(stat2BRBpKp*1.2*1.2+syst2BRBpKp*1.0*1.0)*1.e-9;
names[278]="dG_B+K+mumu_11.8-12.5";

central[279]=16.1e-9;
errors[279]=sqrt(stat2BRBpKp*1.0*1.0+syst2BRBpKp*0.8*0.8)*1.e-9;
names[279]="dG_B+K+mumu_15-16";

central[280]=16.4e-9;
errors[280]=sqrt(stat2BRBpKp*1.0*1.0+syst2BRBpKp*0.8*0.8)*1.e-9;
names[280]="dG_B+K+mumu_16-17";

central[281]=20.6e-9;
errors[281]=sqrt(stat2BRBpKp*1.1*1.1+syst2BRBpKp*1.0*1.0)*1.e-9;
names[281]="dG_B+K+mumu_17-18";

central[282]=13.7e-9;
errors[282]=sqrt(stat2BRBpKp*1.0*1.0+syst2BRBpKp*0.7*0.7)*1.e-9;
names[282]="dG_B+K+mumu_18-19";

central[283]=7.4e-9;
errors[283]=sqrt(stat2BRBpKp*0.8*0.8+syst2BRBpKp*0.4*0.4)*1.e-9;
names[283]="dG_B+K+mumu_19-20";

central[284]=5.9e-9;
errors[284]=sqrt(stat2BRBpKp*0.7*0.7+syst2BRBpKp*0.3*0.3)*1.e-9;
names[284]="dG_B+K+mumu_20-21";

central[285]=4.3e-9;
errors[285]=sqrt(stat2BRBpKp*0.7*0.7+syst2BRBpKp*0.2*0.2)*1.e-9;
names[285]="dG_B+K+mumu_21-22";

central[286]=24.2e-9;
errors[286]=sqrt(stat2BRBpKp*0.7*0.7+syst2BRBpKp*1.2*1.2)*1.e-9;
names[286]="dG_B+K+mumu_1.1-6";

central[475]=0.03;
errors[475]=sqrt(stat2BRBpKp*0.03*0.03+syst2BRBpKp*0.02*0.02);
names[475]="FH_B+K+mumu_1.1-6";

central[287]=12.1e-9;
errors[287]=sqrt(stat2BRBpKp*0.4*0.4+syst2BRBpKp*0.6*0.6)*1.e-9;
names[287]="dG_B+K+mumu_15-22";

central[476]=0.035;
errors[476]=sqrt(stat2BRBpKp*0.035*0.035+syst2BRBpKp*0.02*0.02);
names[476]="FH_B+K+mumu_15-22";

central[288]=-0.;
errors[288]=-0.;
names[288]="dG_B+K+mumu_1-6";

central[289]=1.56e-7;
errors[289]=sqrt(statBRBpKpee*0.19*0.19+systBRBpKpee*0.06*0.06);
names[289]="dG_B+K+ee_1-6";

central[290]=0.745-1.;
errors[290]=sqrt(statRK*0.09*0.09+systRK*0.036*0.036);
names[290]="RK_1-6-1";

central[291]=5.85e-8;
errors[291]=sqrt(statBRBsPhi*0.73*0.73+systBRBsPhi*0.14*0.14+otherBRBsPhi*0.44*0.44)*1.e-8;
names[291]="dG_Bsphimumu_0.1-2";

central[292]=0.20;
errors[292]=sqrt(statAngBsPhi*0.09*0.09+systAngBsPhi*0.02*0.02);
names[292]="FL_Bsphimumu_0.1-2";

central[293]=-0.05;
errors[293]=sqrt(statAngBsPhi*0.13*0.13+systAngBsPhi*0.01*0.01);
names[293]="S3_Bsphimumu_0.1-2";

central[294]=-0.27;
errors[294]=sqrt(statAngBsPhi*0.28*0.28+systAngBsPhi*0.01*0.01);
names[294]="S4_Bsphimumu_0.1-2";

central[295]=0.04;
errors[295]=sqrt(statAngBsPhi*0.12*0.12);
names[295]="S7_Bsphimumu_0.1-2";

central[296]=2.56e-8;
errors[296]=sqrt(statBRBsPhi*0.42*0.42+systBRBsPhi*0.06*0.06+otherBRBsPhi*0.19*0.19)*1.e-8;
names[296]="dG_Bsphimumu_2-5";

central[297]=0.68;
errors[297]=sqrt(statAngBsPhi*0.16*0.16+systAngBsPhi*0.03*0.03);
names[297]="FL_Bsphimumu_2-5";

central[298]=-0.06;
errors[298]=sqrt(statAngBsPhi*0.23*0.23+systAngBsPhi*0.01*0.01);
names[298]="S3_Bsphimumu_2-5";

central[299]=0.47;
errors[299]=sqrt(statAngBsPhi*0.44*0.44+systAngBsPhi*0.01*0.01);
names[299]="S4_Bsphimumu_2-5";

central[300]=-0.03;
errors[300]=sqrt(statAngBsPhi*0.23*0.23+systAngBsPhi*0.01*0.01);
names[300]="S7_Bsphimumu_2-5";

central[301]=3.21e-8;
errors[301]=sqrt(statBRBsPhi*0.44*0.44+systBRBsPhi*0.08*0.08+otherBRBsPhi*0.24*0.24)*1.e-8;
names[301]="dG_Bsphimumu_5-8";

central[302]=0.54;
errors[302]=sqrt(statAngBsPhi*0.10*0.10+systAngBsPhi*0.02*0.02);
names[302]="FL_Bsphimumu_5-8";

central[303]=-0.10;
errors[303]=sqrt(statAngBsPhi*0.29*0.29+systAngBsPhi*0.01*0.01);
names[303]="S3_Bsphimumu_5-8";

central[304]=0.10;
errors[304]=sqrt(statAngBsPhi*0.18*0.18+systAngBsPhi*0.01*0.01);
names[304]="S4_Bsphimumu_5-8";

central[305]=0.04;
errors[305]=sqrt(statAngBsPhi*0.20*0.20+systAngBsPhi*0.01*0.01);
names[305]="S7_Bsphimumu_5-8";

central[306]=4.71e-8;
errors[306]=sqrt(statBRBsPhi*0.69*0.69+systBRBsPhi*0.15*0.15+otherBRBsPhi*0.36*0.36)*1.e-8;
names[306]="dG_Bsphimumu_11-12.5";

central[307]=0.29;
errors[307]=sqrt(statAngBsPhi*0.11*0.11+systAngBsPhi*0.04*0.04);
names[307]="FL_Bsphimumu_11-12.5";

central[308]=-0.19;
errors[308]=sqrt(statAngBsPhi*0.23*0.23+systAngBsPhi*0.01*0.01);
names[308]="S3_Bsphimumu_11-12.5";

central[309]=0.47;
errors[309]=sqrt(statAngBsPhi*0.29*0.29+systAngBsPhi*0.01*0.01);
names[309]="S4_Bsphimumu_11-12.5";

central[310]=0.;
errors[310]=sqrt(statAngBsPhi*0.17*0.17+systAngBsPhi*0.01*0.01);
names[310]="S7_Bsphimumu_11-12.5";

central[311]=4.52e-8;
errors[311]=sqrt(statBRBsPhi*0.57*0.57+systBRBsPhi*0.12*0.12+otherBRBsPhi*0.34*0.34)*1.e-8;
names[311]="dG_Bsphimumu_15-17";

central[312]=0.23;
errors[312]=sqrt(statAngBsPhi*0.09*0.09+systAngBsPhi*0.02*0.02);
names[312]="FL_Bsphimumu_15-17";

central[313]=-0.06;
errors[313]=sqrt(statAngBsPhi*0.19*0.19+systAngBsPhi*0.01*0.01);
names[313]="S3_Bsphimumu_15-17";

central[314]=0.03;
errors[314]=sqrt(statAngBsPhi*0.15*0.15+systAngBsPhi*0.01*0.01);
names[314]="S4_Bsphimumu_15-17";

central[315]=0.12;
errors[315]=sqrt(statAngBsPhi*0.16*0.16+systAngBsPhi*0.01*0.01);
names[315]="S7_Bsphimumu_15-17";

central[316]=3.96e-8;
errors[316]=sqrt(statBRBsPhi*0.57*0.57+systBRBsPhi*0.14*0.14+otherBRBsPhi*0.30*0.30)*1.e-8;
names[316]="dG_Bsphimumu_17-19";

central[317]=0.40;
errors[317]=sqrt(statAngBsPhi*0.15*0.15+systAngBsPhi*0.02*0.02);
names[317]="FL_Bsphimumu_17-19";

central[318]=-0.07;
errors[318]=sqrt(statAngBsPhi*0.27*0.27+systAngBsPhi*0.02*0.02);
names[318]="S3_Bsphimumu_17-19";

central[319]=0.39;
errors[319]=sqrt(statAngBsPhi*0.34*0.34+systAngBsPhi*0.02*0.02);
names[319]="S4_Bsphimumu_17-19";

central[320]=0.20;
errors[320]=sqrt(statAngBsPhi*0.29*0.29+systAngBsPhi*0.01*0.01);
names[320]="S7_Bsphimumu_17-19";

central[321]=2.58e-8;
errors[321]=sqrt(statBRBsPhi*0.33*0.33+systBRBsPhi*0.08*0.08+otherBRBsPhi*0.19*0.19)*1.e-8;
names[321]="dG_Bsphimumu_1-6";

central[322]=0.63;
errors[322]=sqrt(statAngBsPhi*0.09*0.09+systAngBsPhi*0.03*0.03);
names[322]="FL_Bsphimumu_1-6";

central[323]=-0.02;
errors[323]=sqrt(statAngBsPhi*0.13*0.13+systAngBsPhi*0.01*0.01);
names[323]="S3_Bsphimumu_1-6";

central[324]=0.19;
errors[324]=sqrt(statAngBsPhi*0.14*0.14+systAngBsPhi*0.01*0.01);
names[324]="S4_Bsphimumu_1-6";

central[325]=-0.03;
errors[325]=sqrt(statAngBsPhi*0.14*0.14);
names[325]="S7_Bsphimumu_1-6";

central[326]=4.04e-8;
errors[326]=sqrt(statBRBsPhi*0.39*0.39+systBRBsPhi*0.13*0.13+otherBRBsPhi*0.30*0.30)*1.e-8;
names[326]="dG_Bsphimumu_15-19";

central[327]=0.29;
errors[327]=sqrt(statAngBsPhi*0.07*0.07+systAngBsPhi*0.02*0.02);
names[327]="FL_Bsphimumu_15-19";

central[328]=-0.09;
errors[328]=sqrt(statAngBsPhi*0.12*0.12+systAngBsPhi*0.01*0.01);
names[328]="S3_Bsphimumu_15-19";

central[329]=0.14;
errors[329]=sqrt(statAngBsPhi*0.11*0.11+systAngBsPhi*0.01*0.01);
names[329]="S4_Bsphimumu_15-19";

central[330]=0.13;
errors[330]=sqrt(statAngBsPhi*0.11*0.11+systAngBsPhi*0.01*0.01);
names[330]="S7_Bsphimumu_15-19";

central[331]=0.080;
errors[331]=sqrt(stat2AngB0Kstar0*0.248*0.248+syst2AngB0Kstar0*0.044*0.044);
names[331]="P1_B0Kstar0mumu_1.1_6";

central[332]=-0.323;
errors[332]=0.159;
names[332]="ATRe_B0Kstar0mumu_1.1_6";

if(likelihood==1)
{
central[333]=-0.099;
errors[333]=sqrt(stat2AngB0Kstar0*0.168*0.168+syst2AngB0Kstar0*0.014*0.014);
names[333]="P1_B0Kstar0mumu_0.1_0.98";
}
else
{
central[333]=-0.038;
errors[333]=sqrt(stat2AngB0Kstar0*0.158*0.158+syst2AngB0Kstar0*0.020*0.020);
names[333]="P1_B0Kstar0mumu_0.1_0.98";
}

central[334]=0.;
errors[334]=0.;
names[334]="ATRe_B0Kstar0mumu_0.1_0.98";

central[335]=0.;
errors[335]=0.;
names[335]="dG_B0Kstar0ee_1-6";

if(likelihood==1)
{
central[336]=-0.003;
errors[336]=sqrt(stat2AngB0Kstar0*0.052*0.052+syst2AngB0Kstar0*0.007*0.007);
names[336]="P2_B0Kstar0mumu_0.1-0.98";

central[337]=0.113;
errors[337]=sqrt(stat2AngB0Kstar0*0.079*0.079+syst2AngB0Kstar0*0.006*0.006);
names[337]="P3_B0Kstar0mumu_0.1-0.98";

central[338]=0.185;
errors[338]=sqrt(stat2AngB0Kstar0*0.158*0.158+syst2AngB0Kstar0*0.023*0.023);
names[338]="P4p_B0Kstar0mumu_0.1-0.98";

central[339]=0.034;
errors[339]=sqrt(stat2AngB0Kstar0*0.135*0.135+syst2AngB0Kstar0*0.015*0.015);
names[339]="P6p_B0Kstar0mumu_0.1-0.98";

central[340]=0.180;
errors[340]=sqrt(stat2AngB0Kstar0*0.174*0.174+syst2AngB0Kstar0*0.007*0.007);
names[340]="P8p_B0Kstar0mumu_0.1-0.98";
}
else
{
central[336]=-0.119;
errors[336]=sqrt(stat2AngB0Kstar0*0.081*0.081+syst2AngB0Kstar0*0.063*0.063);
names[336]="P2_B0Kstar0mumu_0.1-0.98";

central[337]=0.147;
errors[337]=sqrt(stat2AngB0Kstar0*0.086*0.086+syst2AngB0Kstar0*0.005*0.005);
names[337]="P3_B0Kstar0mumu_0.1-0.98";

central[338]=0.086;
errors[338]=sqrt(stat2AngB0Kstar0*0.221*0.221+syst2AngB0Kstar0*0.026*0.026);
names[338]="P4p_B0Kstar0mumu_0.1-0.98";

central[339]=0.086;
errors[339]=sqrt(stat2AngB0Kstar0*0.152*0.152+syst2AngB0Kstar0*0.024*0.024);
names[339]="P6p_B0Kstar0mumu_0.1-0.98";

central[340]=0.143;
errors[340]=sqrt(stat2AngB0Kstar0*0.195*0.195+syst2AngB0Kstar0*0.022*0.022);
names[340]="P8p_B0Kstar0mumu_0.1-0.98";
}

central[341]=-0.451;
errors[341]=sqrt(stat2AngB0Kstar0*0.636*0.636+syst2AngB0Kstar0*0.038*0.038);
names[341]="P1_B0Kstar0mumu_1.1-2.5";

central[342]=-0.373;
errors[342]=sqrt(stat2AngB0Kstar0*0.199*0.199+syst2AngB0Kstar0*0.027*0.027);
names[342]="P2_B0Kstar0mumu_1.1-2.5";

central[343]=0.350;
errors[343]=sqrt(stat2AngB0Kstar0*0.330*0.330+syst2AngB0Kstar0*0.015*0.015);
names[343]="P3_B0Kstar0mumu_1.1-2.5";

central[344]=-0.163;
errors[344]=sqrt(stat2AngB0Kstar0*0.240*0.240+syst2AngB0Kstar0*0.021*0.021);
names[344]="P4p_B0Kstar0mumu_1.1-2.5";

central[345]=-0.463;
errors[345]=sqrt(stat2AngB0Kstar0*0.221*0.221+syst2AngB0Kstar0*0.012*0.012);
names[345]="P6p_B0Kstar0mumu_1.1-2.5";

central[346]=-0.208;
errors[346]=sqrt(stat2AngB0Kstar0*0.270*0.270+syst2AngB0Kstar0*0.024*0.024);
names[346]="P8p_B0Kstar0mumu_1.1-2.5";

central[347]=0.571;
errors[347]=sqrt(stat2AngB0Kstar0*2.404*2.404+syst2AngB0Kstar0*0.045*0.045);
names[347]="P1_B0Kstar0mumu_2.5-4";

central[348]=-0.636;
errors[348]=sqrt(stat2AngB0Kstar0*1.735*1.735+syst2AngB0Kstar0*0.015*0.015);
names[348]="P2_B0Kstar0mumu_2.5-4";

central[349]=0.745;
errors[349]=sqrt(stat2AngB0Kstar0*2.587*2.587+syst2AngB0Kstar0*0.030*0.030);
names[349]="P3_B0Kstar0mumu_2.5-4";

central[350]=-0.713;
errors[350]=sqrt(stat2AngB0Kstar0*1.305*1.305+syst2AngB0Kstar0*0.024*0.024);
names[350]="P4p_B0Kstar0mumu_2.5-4";

central[351]=0.205;
errors[351]=sqrt(stat2AngB0Kstar0*0.962*0.962+syst2AngB0Kstar0*0.013*0.013);
names[351]="P6p_B0Kstar0mumu_2.5-4";

central[352]=0.091;
errors[352]=sqrt(stat2AngB0Kstar0*0.650*0.650+syst2AngB0Kstar0*0.025*0.025);
names[352]="P8p_B0Kstar0mumu_2.5-4";

central[353]=0.180;
errors[353]=sqrt(stat2AngB0Kstar0*0.364*0.364+syst2AngB0Kstar0*0.027*0.027);
names[353]="P1_B0Kstar0mumu_4-6";

central[354]=0.042;
errors[354]=sqrt(stat2AngB0Kstar0*0.088*0.088+syst2AngB0Kstar0*0.011*0.011);
names[354]="P2_B0Kstar0mumu_4-6";

central[355]=0.083;
errors[355]=sqrt(stat2AngB0Kstar0*0.187*0.187+syst2AngB0Kstar0*0.023*0.023);
names[355]="P3_B0Kstar0mumu_4-6";

central[356]=-0.448;
errors[356]=sqrt(stat2AngB0Kstar0*0.172*0.172+syst2AngB0Kstar0*0.020*0.020);
names[356]="P4p_B0Kstar0mumu_4-6";

central[357]=-0.032;
errors[357]=sqrt(stat2AngB0Kstar0*0.167*0.167+syst2AngB0Kstar0*0.007*0.007);
names[357]="P6p_B0Kstar0mumu_4-6";

central[358]=0.342;
errors[358]=sqrt(stat2AngB0Kstar0*0.188*0.188+syst2AngB0Kstar0*0.009*0.009);
names[358]="P8p_B0Kstar0mumu_4-6";

central[359]=-0.199;
errors[359]=sqrt(stat2AngB0Kstar0*0.281*0.281+syst2AngB0Kstar0*0.025*0.025);
names[359]="P1_B0Kstar0mumu_6-8";

central[360]=0.241;
errors[360]=sqrt(stat2AngB0Kstar0*0.062*0.062+syst2AngB0Kstar0*0.013*0.013);
names[360]="P2_B0Kstar0mumu_6-8";

central[361]=0.057;
errors[361]=sqrt(stat2AngB0Kstar0*0.148*0.148+syst2AngB0Kstar0*0.013*0.013);
names[361]="P3_B0Kstar0mumu_6-8";

central[362]=-0.599;
errors[362]=sqrt(stat2AngB0Kstar0*0.135*0.135+syst2AngB0Kstar0*0.010*0.010);
names[362]="P4p_B0Kstar0mumu_6-8";

central[363]=-0.095;
errors[363]=sqrt(stat2AngB0Kstar0*0.135*0.135+syst2AngB0Kstar0*0.011*0.011);
names[363]="P6p_B0Kstar0mumu_6-8";

central[364]=-0.171;
errors[364]=sqrt(stat2AngB0Kstar0*0.143*0.143+syst2AngB0Kstar0*0.006*0.006);
names[364]="P8p_B0Kstar0mumu_6-8";

central[365]=-0.745;
errors[365]=sqrt(stat2AngB0Kstar0*0.230*0.230+syst2AngB0Kstar0*0.015*0.015);
names[365]="P1_B0Kstar0mumu_11-12.5";

central[366]=0.418;
errors[366]=sqrt(stat2AngB0Kstar0*0.053*0.053+syst2AngB0Kstar0*0.005*0.005);
names[366]="P2_B0Kstar0mumu_11-12.5";

central[367]=0.007;
errors[367]=sqrt(stat2AngB0Kstar0*0.141*0.141+syst2AngB0Kstar0*0.010*0.010);
names[367]="P3_B0Kstar0mumu_11-12.5";

central[368]=-0.567;
errors[368]=sqrt(stat2AngB0Kstar0*0.187*0.187+syst2AngB0Kstar0*0.014*0.014);
names[368]="P4p_B0Kstar0mumu_11-12.5";

central[369]=-0.282;
errors[369]=sqrt(stat2AngB0Kstar0*0.151*0.151+syst2AngB0Kstar0*0.007*0.007);
names[369]="P6p_B0Kstar0mumu_11-12.5";

central[370]=-0.015;
errors[370]=sqrt(stat2AngB0Kstar0*0.145*0.145+syst2AngB0Kstar0*0.005*0.005);
names[370]="P8p_B0Kstar0mumu_11-12.5";

central[371]=-0.436;
errors[371]=sqrt(stat2AngB0Kstar0*0.147*0.147+syst2AngB0Kstar0*0.018*0.018);
names[371]="P1_B0Kstar0mumu_15-17";

central[372]=0.421;
errors[372]=sqrt(stat2AngB0Kstar0*0.042*0.042+syst2AngB0Kstar0*0.005*0.005);
names[372]="P2_B0Kstar0mumu_15-17";

central[373]=0.029;
errors[373]=sqrt(stat2AngB0Kstar0*0.084*0.084+syst2AngB0Kstar0*0.006*0.006);
names[373]="P3_B0Kstar0mumu_15-17";

central[374]=-0.672;
errors[374]=sqrt(stat2AngB0Kstar0*0.151*0.151+syst2AngB0Kstar0*0.016*0.016);
names[374]="P4p_B0Kstar0mumu_15-17";

central[375]=0.127;
errors[375]=sqrt(stat2AngB0Kstar0*0.122*0.122+syst2AngB0Kstar0*0.006*0.006);
names[375]="P6p_B0Kstar0mumu_15-17";

central[376]=0.007;
errors[376]=sqrt(stat2AngB0Kstar0*0.129*0.129+syst2AngB0Kstar0*0.005*0.005);
names[376]="P8p_B0Kstar0mumu_15-17";

central[377]=-0.581;
errors[377]=sqrt(stat2AngB0Kstar0*0.263*0.263+syst2AngB0Kstar0*0.037*0.037);
names[377]="P1_B0Kstar0mumu_17-19";

central[378]=0.314;
errors[378]=sqrt(stat2AngB0Kstar0*0.048*0.048+syst2AngB0Kstar0*0.007*0.007);
names[378]="P2_B0Kstar0mumu_17-19";

central[379]=0.145;
errors[379]=sqrt(stat2AngB0Kstar0*0.107*0.107+syst2AngB0Kstar0*0.008*0.008);
names[379]="P3_B0Kstar0mumu_17-19";

central[380]=-0.556;
errors[380]=sqrt(stat2AngB0Kstar0*0.156*0.156+syst2AngB0Kstar0*0.016*0.016);
names[380]="P4p_B0Kstar0mumu_17-19";

central[381]=0.092;
errors[381]=sqrt(stat2AngB0Kstar0*0.152*0.152+syst2AngB0Kstar0*0.025*0.025);
names[381]="P6p_B0Kstar0mumu_17-19";

central[382]=0.027;
errors[382]=sqrt(stat2AngB0Kstar0*0.147*0.147+syst2AngB0Kstar0*0.009*0.009);
names[382]="P8p_B0Kstar0mumu_17-19";

central[383]=-0.162;
errors[383]=sqrt(stat2AngB0Kstar0*0.073*0.073+syst2AngB0Kstar0*0.010*0.010);
names[383]="P2_B0Kstar0mumu_1.1_6";

central[384]=0.205;
errors[384]=sqrt(stat2AngB0Kstar0*0.135*0.135+syst2AngB0Kstar0*0.017*0.017);
names[384]="P3_B0Kstar0mumu_1.1_6";

central[385]=-0.336;
errors[385]=sqrt(stat2AngB0Kstar0*0.124*0.124+syst2AngB0Kstar0*0.012*0.012);
names[385]="P4p_B0Kstar0mumu_1.1_6";

central[386]=-0.166;
errors[386]=sqrt(stat2AngB0Kstar0*0.108*0.108+syst2AngB0Kstar0*0.021*0.021);
names[386]="P6p_B0Kstar0mumu_1.1_6";

central[387]=0.060;
errors[387]=sqrt(stat2AngB0Kstar0*0.124*0.124+syst2AngB0Kstar0*0.009*0.009);
names[387]="P8p_B0Kstar0mumu_1.1_6";

if(likelihood==1)
{
central[388]=-0.497;
errors[388]=sqrt(stat2AngB0Kstar0*0.102*0.102+syst2AngB0Kstar0*0.027*0.027);
names[388]="P1_B0Kstar0mumu_15-19";

central[389]=0.361;
errors[389]=sqrt(stat2AngB0Kstar0*0.026*0.026+syst2AngB0Kstar0*0.010*0.010);
names[389]="P2_B0Kstar0mumu_15-19";

central[390]=0.081;
errors[390]=sqrt(stat2AngB0Kstar0*0.060*0.060+syst2AngB0Kstar0*0.005*0.005);
names[390]="P3_B0Kstar0mumu_15-19";

central[391]=-0.597;
errors[391]=sqrt(stat2AngB0Kstar0*0.085*0.085+syst2AngB0Kstar0*0.015*0.015);
names[391]="P4p_B0Kstar0mumu_15-19";

central[392]=0.101;
errors[392]=sqrt(stat2AngB0Kstar0*0.092*0.092+syst2AngB0Kstar0*0.011*0.011);
names[392]="P6p_B0Kstar0mumu_15-19";

central[393]=0.059;
errors[393]=sqrt(stat2AngB0Kstar0*0.094*0.094+syst2AngB0Kstar0*0.008*0.008);
names[393]="P8p_B0Kstar0mumu_15-19";
}
else
{
central[388]=-0.424;
errors[388]=sqrt(stat2AngB0Kstar0*0.150*0.150+syst2AngB0Kstar0*0.028*0.028);
names[388]="P1_B0Kstar0mumu_15-19";

central[389]=0.385;
errors[389]=sqrt(stat2AngB0Kstar0*0.036*0.036+syst2AngB0Kstar0*0.010*0.010);
names[389]="P2_B0Kstar0mumu_15-19";

central[390]=0.089;
errors[390]=sqrt(stat2AngB0Kstar0*0.072*0.072+syst2AngB0Kstar0*0.019*0.019);
names[390]="P3_B0Kstar0mumu_15-19";

central[391]=-0.663;
errors[391]=sqrt(stat2AngB0Kstar0*0.105*0.105+syst2AngB0Kstar0*0.055*0.055);
names[391]="P4p_B0Kstar0mumu_15-19";

central[392]=0.140;
errors[392]=sqrt(stat2AngB0Kstar0*0.101*0.101+syst2AngB0Kstar0*0.032*0.032);
names[392]="P6p_B0Kstar0mumu_15-19";

central[393]=0.049;
errors[393]=sqrt(stat2AngB0Kstar0*0.106*0.106+syst2AngB0Kstar0*0.021*0.021);
names[393]="P8p_B0Kstar0mumu_15-19";
}

central[394]=0.439;
errors[394]=sqrt(stat2AngB0Kstar0*1.916*1.916+syst2AngB0Kstar0*0.012*0.012);
names[394]="P1_B0Kstar0mumu_1.1-2";

central[395]=-0.667;
errors[395]=sqrt(stat2AngB0Kstar0*1.939*1.939+syst2AngB0Kstar0*0.017*0.017);
names[395]="P2_B0Kstar0mumu_1.1-2";

central[396]=0.363;
errors[396]=sqrt(stat2AngB0Kstar0*1.088*1.088+syst2AngB0Kstar0*0.001*0.001);
names[396]="P3_B0Kstar0mumu_1.1-2";

central[397]=0.266;
errors[397]=sqrt(stat2AngB0Kstar0*0.648*0.648+syst2AngB0Kstar0*0.057*0.057);
names[397]="P4p_B0Kstar0mumu_1.1-2";

central[398]=-0.632;
errors[398]=sqrt(stat2AngB0Kstar0*0.753*0.753+syst2AngB0Kstar0*0.009*0.009);
names[398]="P6p_B0Kstar0mumu_1.1-2";

central[399]=-0.244;
errors[399]=sqrt(stat2AngB0Kstar0*0.645*0.645+syst2AngB0Kstar0*0.012*0.012);
names[399]="P8p_B0Kstar0mumu_1.1-2";

central[400]=0.055;
errors[400]=sqrt(stat2AngB0Kstar0*0.756*0.756+syst2AngB0Kstar0*0.007*0.007);
names[400]="P1_B0Kstar0mumu_2-3";

central[401]=-0.323;
errors[401]=sqrt(stat2AngB0Kstar0*0.316*0.316+syst2AngB0Kstar0*0.033*0.033);
names[401]="P2_B0Kstar0mumu_2-3";

central[402]=0.005;
errors[402]=sqrt(stat2AngB0Kstar0*0.364*0.364+syst2AngB0Kstar0*0.012*0.012);
names[402]="P3_B0Kstar0mumu_2-3";

central[403]=-0.765;
errors[403]=sqrt(stat2AngB0Kstar0*0.359*0.359+syst2AngB0Kstar0*0.099*0.099);
names[403]="P4p_B0Kstar0mumu_2-3";

central[404]=-0.549;
errors[404]=sqrt(stat2AngB0Kstar0*0.393*0.393+syst2AngB0Kstar0*0.005*0.005);
names[404]="P6p_B0Kstar0mumu_2-3";

central[405]=-0.393;
errors[405]=sqrt(stat2AngB0Kstar0*0.388*0.388+syst2AngB0Kstar0*0.002*0.002);
names[405]="P8p_B0Kstar0mumu_2-3";

central[406]=0.421;
errors[406]=sqrt(stat2AngB0Kstar0*18.35*18.25+syst2AngB0Kstar0*0.018*0.018);
names[406]="P1_B0Kstar0mumu_3-4";

central[407]=-0.117;
errors[407]=sqrt(stat2AngB0Kstar0*4.435*4.435+syst2AngB0Kstar0*0.015*0.015);
names[407]="P2_B0Kstar0mumu_3-4";

central[408]=0.905;
errors[408]=sqrt(stat2AngB0Kstar0*17.51*17.51+syst2AngB0Kstar0*0.009*0.009);
names[408]="P3_B0Kstar0mumu_3-4";

central[409]=-0.134;
errors[409]=sqrt(stat2AngB0Kstar0*1.343*1.343+syst2AngB0Kstar0*0.108*0.108);
names[409]="P4p_B0Kstar0mumu_3-4";

central[410]=0.449;
errors[410]=sqrt(stat2AngB0Kstar0*19.04*19.04+syst2AngB0Kstar0*0.007*0.007);
names[410]="P6p_B0Kstar0mumu_3-4";

central[411]=0.303;
errors[411]=sqrt(stat2AngB0Kstar0*1.394*1.394+syst2AngB0Kstar0*0.006*0.006);
names[411]="P8p_B0Kstar0mumu_3-4";

central[412]=2.296;
errors[412]=sqrt(stat2AngB0Kstar0*17.71*17.71+syst2AngB0Kstar0*0.024*0.024);
names[412]="P1_B0Kstar0mumu_4-5";

central[413]=0.174;
errors[413]=sqrt(stat2AngB0Kstar0*3.034*3.034+syst2AngB0Kstar0*0.010*0.010);
names[413]="P2_B0Kstar0mumu_4-5";

central[414]=-0.801;
errors[414]=sqrt(stat2AngB0Kstar0*17.42*17.42+syst2AngB0Kstar0*0.007*0.007);
names[414]="P3_B0Kstar0mumu_4-5";

central[415]=-0.415;
errors[415]=sqrt(stat2AngB0Kstar0*1.911*1.911+syst2AngB0Kstar0*0.104*0.104);
names[415]="P4p_B0Kstar0mumu_4-5";

central[416]=-0.215;
errors[416]=sqrt(stat2AngB0Kstar0*1.243*1.243+syst2AngB0Kstar0*0.006*0.006);
names[416]="P6p_B0Kstar0mumu_4-5";

central[417]=0.293;
errors[417]=sqrt(stat2AngB0Kstar0*1.522*1.522+syst2AngB0Kstar0*0.006*0.006);
names[417]="P8p_B0Kstar0mumu_4-5";

central[418]=-0.540;
errors[418]=sqrt(stat2AngB0Kstar0*1.100*1.100+syst2AngB0Kstar0*0.025*0.025);
names[418]="P1_B0Kstar0mumu_5-6";

central[419]=0.089;
errors[419]=sqrt(stat2AngB0Kstar0*0.227*0.227+syst2AngB0Kstar0*0.012*0.012);
names[419]="P2_B0Kstar0mumu_5-6";

central[420]=0.178;
errors[420]=sqrt(stat2AngB0Kstar0*0.465*0.465+syst2AngB0Kstar0*0.007*0.007);
names[420]="P3_B0Kstar0mumu_5-6";

central[421]=-0.561;
errors[421]=sqrt(stat2AngB0Kstar0*0.465*0.465+syst2AngB0Kstar0*0.101*0.101);
names[421]="P4p_B0Kstar0mumu_5-6";

central[422]=0.074;
errors[422]=sqrt(stat2AngB0Kstar0*0.309*0.309+syst2AngB0Kstar0*0.005*0.005);
names[422]="P6p_B0Kstar0mumu_5-6";

central[423]=-0.068;
errors[423]=sqrt(stat2AngB0Kstar0*0.372*0.372+syst2AngB0Kstar0*0.006*0.006);
names[423]="P8p_B0Kstar0mumu_5-6";

central[424]=-0.353;
errors[424]=sqrt(stat2AngB0Kstar0*0.602*0.602+syst2AngB0Kstar0*0.026*0.026);
names[424]="P1_B0Kstar0mumu_6-7";

central[425]=0.104;
errors[425]=sqrt(stat2AngB0Kstar0*0.136*0.136+syst2AngB0Kstar0*0.013*0.013);
names[425]="P2_B0Kstar0mumu_6-7";

central[426]=-0.161;
errors[426]=sqrt(stat2AngB0Kstar0*0.291*0.291+syst2AngB0Kstar0*0.001*0.001);
names[426]="P3_B0Kstar0mumu_6-7";

central[427]=-0.641;
errors[427]=sqrt(stat2AngB0Kstar0*0.294*0.294+syst2AngB0Kstar0*0.106*0.106);
names[427]="P4p_B0Kstar0mumu_6-7";

central[428]=0.017;
errors[428]=sqrt(stat2AngB0Kstar0*0.267*0.267+syst2AngB0Kstar0*0.007*0.007);
names[428]="P6p_B0Kstar0mumu_6-7";

central[429]=0.162;
errors[429]=sqrt(stat2AngB0Kstar0*0.289*0.289+syst2AngB0Kstar0*0.005*0.005);
names[429]="P8p_B0Kstar0mumu_6-7";

central[430]=-0.284;
errors[430]=sqrt(stat2AngB0Kstar0*0.548*0.548+syst2AngB0Kstar0*0.025*0.025);
names[430]="P1_B0Kstar0mumu_7-8";

central[431]=0.393;
errors[431]=sqrt(stat2AngB0Kstar0*0.231*0.231+syst2AngB0Kstar0*0.013*0.013);
names[431]="P2_B0Kstar0mumu_7-8";

central[432]=-0.063;
errors[432]=sqrt(stat2AngB0Kstar0*0.298*0.298+syst2AngB0Kstar0*0.002*0.002);
names[432]="P3_B0Kstar0mumu_7-8";

central[433]=-0.503;
errors[433]=sqrt(stat2AngB0Kstar0*0.288*0.288+syst2AngB0Kstar0*0.118*0.118);
names[433]="P4p_B0Kstar0mumu_7-8";

central[434]=-0.201;
errors[434]=sqrt(stat2AngB0Kstar0*0.274*0.274+syst2AngB0Kstar0*0.007*0.007);
names[434]="P6p_B0Kstar0mumu_7-8";

central[435]=-0.623;
errors[435]=sqrt(stat2AngB0Kstar0*0.295*0.295+syst2AngB0Kstar0*0.005*0.005);
names[435]="P8p_B0Kstar0mumu_7-8";

central[436]=-0.869;
errors[436]=sqrt(stat2AngB0Kstar0*0.408*0.408+syst2AngB0Kstar0*0.030*0.030);
names[436]="P1_B0Kstar0mumu_11-11.75";

central[437]=0.494;
errors[437]=sqrt(stat2AngB0Kstar0*0.134*0.134+syst2AngB0Kstar0*0.013*0.013);
names[437]="P2_B0Kstar0mumu_11-11.75";

central[438]=0.166;
errors[438]=sqrt(stat2AngB0Kstar0*0.221*0.221+syst2AngB0Kstar0*0.005*0.005);
names[438]="P3_B0Kstar0mumu_11-11.75";

central[439]=-0.522;
errors[439]=sqrt(stat2AngB0Kstar0*0.222*0.222+syst2AngB0Kstar0*0.128*0.128);
names[439]="P4p_B0Kstar0mumu_11-11.75";

central[440]=-0.233;
errors[440]=sqrt(stat2AngB0Kstar0*0.227*0.227+syst2AngB0Kstar0*0.004*0.004);
names[440]="P6p_B0Kstar0mumu_11-11.75";

central[441]=-0.159;
errors[441]=sqrt(stat2AngB0Kstar0*0.250*0.250+syst2AngB0Kstar0*0.007*0.007);
names[441]="P8p_B0Kstar0mumu_11-11.75";

central[442]=-1.002;
errors[442]=sqrt(stat2AngB0Kstar0*1.360*1.360+syst2AngB0Kstar0*0.030*0.030);
names[442]="P1_B0Kstar0mumu_11.75-12.5";

central[443]=0.637;
errors[443]=sqrt(stat2AngB0Kstar0*0.599*0.599+syst2AngB0Kstar0*0.008*0.008);
names[443]="P2_B0Kstar0mumu_11.75-12.5";

central[444]=-0.105;
errors[444]=sqrt(stat2AngB0Kstar0*0.42*0.42+syst2AngB0Kstar0*0.004*0.004);
names[444]="P3_B0Kstar0mumu_11.75-12.5";

central[445]=-0.701;
errors[445]=sqrt(stat2AngB0Kstar0*0.342*0.342+syst2AngB0Kstar0*0.114*0.114);
names[445]="P4p_B0Kstar0mumu_11.75-12.5";

central[446]=-0.473;
errors[446]=sqrt(stat2AngB0Kstar0*0.344*0.344+syst2AngB0Kstar0*0.004*0.004);
names[446]="P6p_B0Kstar0mumu_11.75-12.5";

central[447]=-0.211;
errors[447]=sqrt(stat2AngB0Kstar0*0.274*0.274+syst2AngB0Kstar0*0.007*0.007);
names[447]="P8p_B0Kstar0mumu_11.75-12.5";

central[448]=-0.199;
errors[448]=sqrt(stat2AngB0Kstar0*0.285*0.285+syst2AngB0Kstar0*0.014*0.014);
names[448]="P1_B0Kstar0mumu_15-16";

central[449]=0.433;
errors[449]=sqrt(stat2AngB0Kstar0*0.074*0.074+syst2AngB0Kstar0*0.005*0.005);
names[449]="P2_B0Kstar0mumu_15-16";

central[450]=0.087;
errors[450]=sqrt(stat2AngB0Kstar0*0.144*0.144+syst2AngB0Kstar0*0.007*0.007);
names[450]="P3_B0Kstar0mumu_15-16";

central[451]=-0.673;
errors[451]=sqrt(stat2AngB0Kstar0*0.199*0.199+syst2AngB0Kstar0*0.013*0.013);
names[451]="P4p_B0Kstar0mumu_15-16";

central[452]=0.083;
errors[452]=sqrt(stat2AngB0Kstar0*0.189*0.189+syst2AngB0Kstar0*0.004*0.004);
names[452]="P6p_B0Kstar0mumu_15-16";

central[453]=-0.120;
errors[453]=sqrt(stat2AngB0Kstar0*0.198*0.198+syst2AngB0Kstar0*0.010*0.010);
names[453]="P8p_B0Kstar0mumu_15-16";

central[454]=-0.726;
errors[454]=sqrt(stat2AngB0Kstar0*0.241*0.241+syst2AngB0Kstar0*0.014*0.014);
names[454]="P1_B0Kstar0mumu_16-17";

central[455]=0.430;
errors[455]=sqrt(stat2AngB0Kstar0*0.063*0.063+syst2AngB0Kstar0*0.007*0.007);
names[455]="P2_B0Kstar0mumu_16-17";

central[456]=0.019;
errors[456]=sqrt(stat2AngB0Kstar0*0.122*0.122+syst2AngB0Kstar0*0.006*0.006);
names[456]="P3_B0Kstar0mumu_16-17";

central[457]=-0.552;
errors[457]=sqrt(stat2AngB0Kstar0*0.213*0.213+syst2AngB0Kstar0*0.055*0.055);
names[457]="P4p_B0Kstar0mumu_16-17";

central[458]=0.328;
errors[458]=sqrt(stat2AngB0Kstar0*0.195*0.195+syst2AngB0Kstar0*0.012*0.012);
names[458]="P6p_B0Kstar0mumu_16-17";

central[459]=0.122;
errors[459]=sqrt(stat2AngB0Kstar0*0.199*0.199+syst2AngB0Kstar0*0.010*0.010);
names[459]="P8p_B0Kstar0mumu_16-17";

central[460]=-0.313;
errors[460]=sqrt(stat2AngB0Kstar0*0.293*0.293+syst2AngB0Kstar0*0.019*0.019);
names[460]="P1_B0Kstar0mumu_17-18";

central[461]=0.288;
errors[461]=sqrt(stat2AngB0Kstar0*0.075*0.075+syst2AngB0Kstar0*0.006*0.006);
names[461]="P2_B0Kstar0mumu_17-18";

central[462]=0.144;
errors[462]=sqrt(stat2AngB0Kstar0*0.149*0.149+syst2AngB0Kstar0*0.002*0.002);
names[462]="P3_B0Kstar0mumu_17-18";

central[463]=-0.486;
errors[463]=sqrt(stat2AngB0Kstar0*0.200*0.200+syst2AngB0Kstar0*0.092*0.092);
names[463]="P4p_B0Kstar0mumu_17-18";

central[464]=0.047;
errors[464]=sqrt(stat2AngB0Kstar0*0.198*0.198+syst2AngB0Kstar0*0.023*0.023);
names[464]="P6p_B0Kstar0mumu_17-18";

central[465]=-0.006;
errors[465]=sqrt(stat2AngB0Kstar0*0.215*0.215+syst2AngB0Kstar0*0.001*0.001);
names[465]="P8p_B0Kstar0mumu_17-18";

central[466]=-0.450;
errors[466]=sqrt(stat2AngB0Kstar0*0.447*0.447+syst2AngB0Kstar0*0.022*0.022);
names[466]="P1_B0Kstar0mumu_18-19";

central[467]=0.393;
errors[467]=sqrt(stat2AngB0Kstar0*0.159*0.159+syst2AngB0Kstar0*0.011*0.011);
names[467]="P2_B0Kstar0mumu_18-19";

central[468]=0.134;
errors[468]=sqrt(stat2AngB0Kstar0*0.219*0.219+syst2AngB0Kstar0*0.010*0.010);
names[468]="P3_B0Kstar0mumu_18-19";

central[469]=-1.221;
errors[469]=sqrt(stat2AngB0Kstar0*0.388*0.388+syst2AngB0Kstar0*0.119*0.119);
names[469]="P4p_B0Kstar0mumu_18-19";

central[470]=0.128;
errors[470]=sqrt(stat2AngB0Kstar0*0.265*0.265+syst2AngB0Kstar0*0.012*0.012);
names[470]="P6p_B0Kstar0mumu_18-19";

central[471]=0.300;
errors[471]=sqrt(stat2AngB0Kstar0*0.297*0.297+syst2AngB0Kstar0*0.022*0.022);
names[471]="P8p_B0Kstar0mumu_18-19";

central[472]=0.660-1.;
errors[472]=sqrt(stat2RKstar*0.110*0.110+syst2RKstar*0.024*0.024);
names[472]="RKstar_0.045-1.1-1";

central[473]=0.685-1.;
errors[473]=sqrt(stat2RKstar*0.113*0.113+syst2RKstar*0.047*0.047);
names[473]="RKstar_1.1-6-1";

central[474]=4.33e-5;
errors[474]=0.15e-5;
names[474]="BR_B0Kstar0gamma";

