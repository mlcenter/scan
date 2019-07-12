#include "include.h"

#define ORDER 10

#define CACHE

#ifdef CACHE
	#define NTAB 100
	#define NVAL 100
#endif

/*--------------------------------------------------------------------*/

double complex F17_bsll(double s, double z, double L)
{
#ifdef CACHE
	return F_17re_cache(L,z,s,ORDER)+I*F_17im_cache(L,z,s,ORDER);
#else
	return F_17re(L,z,s,ORDER)+I*F_17im(L,z,s,ORDER);
#endif
}

double complex F27_bsll(double s, double z, double L)
{
#ifdef CACHE
	return F_27re_cache(L,z,s,ORDER)+I*F_27im_cache(L,z,s,ORDER);
#else
	return F_27re(L,z,s,ORDER)+I*F_27im(L,z,s,ORDER);
#endif
}

double complex F19_bsll(double s, double z, double L)
{
#ifdef CACHE
	return F_19re_cache(L,z,s,ORDER)+DeltaF_19re_cache(L,z,s,ORDER)+I*(F_19im_cache(L,z,s,ORDER)+DeltaF_19im_cache(L,z,s,ORDER));
#else
	return F_19re(L,z,s,ORDER)+DeltaF_19re(L,z,s,ORDER)+I*(F_19im(L,z,s,ORDER)+DeltaF_19im(L,z,s,ORDER));
#endif
}

double complex F29_bsll(double s, double z, double L)
{
#ifdef CACHE
	return F_29re_cache(L,z,s,ORDER)+DeltaF_29re_cache(L,z,s,ORDER)+I*(F_29im_cache(L,z,s,ORDER)+DeltaF_29im_cache(L,z,s,ORDER));
#else
	return F_29re(L,z,s,ORDER)+DeltaF_29re(L,z,s,ORDER)+I*(F_29im(L,z,s,ORDER)+DeltaF_29im(L,z,s,ORDER));
#endif
}

double complex F19_bkll(double s, double z, double L)
{
#ifdef CACHE
	return F_19re_cache(L,z,s,ORDER)+I*F_19im_cache(L,z,s,ORDER);
#else
	return F_19re(L,z,s,ORDER)+I*F_19im(L,z,s,ORDER);
#endif
}

double complex F29_bkll(double s, double z, double L)
{
#ifdef CACHE
	return F_29re_cache(L,z,s,ORDER)+I*F_29im_cache(L,z,s,ORDER);
#else
	return F_29re(L,z,s,ORDER)+I*F_29im(L,z,s,ORDER);
#endif
}

/*--------------------------------------------------------------------*/

double F_17re(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double Li3shb = creal(CLi3(1.-sh));
	double Li4sh = creal(CLi4(sh));
	double sqrtsh = sqrt(sh);
	double sqrt4sh = sqrt(4-sh);
	double ash = asin(sqrtsh/2);
	double cl2 = Cl2(2*ash);
	double cl3 = Cl3(2*ash);
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double sh_fourteen = sh*sh_thirteen;
	double sh_fifteen = sh*sh_fourteen;
	double sh_sixteen = sh*sh_fifteen;
	double sh_seventeen = sh*sh_sixteen;
	double sh_eighteen = sh*sh_seventeen;
	double sh_nineteen = sh*sh_eighteen;
	double sh_twenty = sh*sh_nineteen;
	double sh_twenty_one = sh*sh_twenty;
	double lz_sq = pow(lz,2);
	double lz_cube = lz*lz_sq;
	double lz_four = lz*lz_cube;
	double z_sqrt = pow(z,0.5);   
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;      
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);        
	double sh_m4_sh_sqrt = pow(-((-4 + sh)*sh),0.5);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = -1.14266 - 0.8559670781893005*L - 2.07503*sh - 
	0.0164609*Lsh*sh - 0.222222*lz*sh - 2.20356*z - 
	2.608715*lz*z - 25.9259*sh*z + 1.03704*Lsh*sh*z - 
	1.700505*lz*sh*z + 0.2962965*Lsh*lz*sh*z + 
	0.14814825*z*lz_sq - 1.122925*sh*z*lz_sq + 
	0.04938275*z*lz_cube + 0.04938275*sh*z*lz_cube - 
	0.024691375*sh*z*lz_four - 19.4691*sh_sq - 
	0.0164609*Lsh*sh_sq - 5.83895*lz*sh_sq - 
	90.4953*z*sh_sq + 2.37037*Lsh*z*sh_sq + 
	7.46645*lz*z*sh_sq + 0.592595*Lsh*lz*z*sh_sq - 
	0.74074*lz_sq*sh_sq - 
	6.1095*z*lz_sq*sh_sq - 
	0.04938275*lz_cube*sh_sq + 
	0.14814875*z*lz_cube*sh_sq - 
	0.074074375*z*lz_four*sh_sq - 41.9952*sh_cube - 
	0.0164609*Lsh*sh_cube - 15.10455*lz*sh_cube - 
	189.354*z*sh_cube + 3.85185*Lsh*z*sh_cube + 
	21.3283*lz*z*sh_cube + 0.88889*Lsh*lz*z*sh_cube - 
	1.555555*lz_sq*sh_cube - 
	14.44125*z*lz_sq*sh_cube - 
	0.14814875*lz_cube*sh_cube + 
	0.34567875*z*lz_cube*sh_cube - 
	0.148148125*z*lz_four*sh_cube - 
	0.00010778*sh_cube*z_m2 + 
	0.00555556*sh_sq*z_m1 + 
	0.946811*sh_cube*z_m1 + 
	0.2444445*lz*sh_cube*z_m1 + 
	0.02469135*lz_sq*sh_cube*z_m1 + 
	11.6973*sh*z_sqrt + 70.1839*sh_sq*z_sqrt + 
	sh_cube*(-3.8991*pow(z,-0.5) + 159.863*z_sqrt) + 
	1.94955*pow(z,1.5) + 1.86366*z_sq - 
	2.331735*lz*z_sq + 11.4229*sh*z_sq - 
	4.66347*Lsh*sh*z_sq - 17.0403*lz*sh*z_sq + 
	2.5926*sh*lz_sq*z_sq + 
	0.5925925*Lsh*sh*lz_sq*z_sq + 
	0.04938275*lz_cube*z_sq + 
	0.29629625*sh*lz_cube*z_sq + 
	23.8816*sh_sq*z_sq - 
	13.9904*Lsh*sh_sq*z_sq - 
	41.39575*lz*sh_sq*z_sq + 
	1.185185*Lsh*lz*sh_sq*z_sq + 
	8.074075*lz_sq*sh_sq*z_sq + 
	1.7777775*Lsh*lz_sq*sh_sq*z_sq + 
	0.74074125*lz_cube*sh_sq*z_sq + 
	45.1784*sh_cube*z_sq - 
	27.3882*Lsh*sh_cube*z_sq - 
	72.5905*lz*sh_cube*z_sq + 
	4.14815*Lsh*lz*sh_cube*z_sq + 
	17.7284*lz_sq*sh_cube*z_sq + 
	3.55555*Lsh*lz_sq*sh_cube*z_sq + 
	1.3827125*lz_cube*sh_cube*z_sq - 
	1.21131*z_cube + 1.49794*lz*z_cube + 
	11.7509*sh*z_cube + 6.73754*Lsh*sh*z_cube + 
	9.4782*lz*sh*z_cube + 0.592595*Lsh*lz*sh*z_cube - 
	1.0370375*lz_sq*z_cube - 
	3.654325*sh*lz_sq*z_cube - 
	0.5925925*Lsh*sh*lz_sq*z_cube + 
	38.1415*sh_sq*z_cube + 
	27.5428*Lsh*sh_sq*z_cube + 
	19.3218*lz*sh_sq*z_cube + 
	1.185185*Lsh*lz*sh_sq*z_cube - 
	10.39505*lz_sq*sh_sq*z_cube - 
	2.37037*Lsh*lz_sq*sh_sq*z_cube + 
	77.3602*sh_cube*z_cube + 
	69.4495*Lsh*sh_cube*z_cube + 
	29.22455*lz*sh_cube*z_cube + 
	0.592595*Lsh*lz*sh_cube*z_cube - 
	24.0247*lz_sq*sh_cube*z_cube - 
	5.925925*Lsh*lz_sq*sh_cube*z_cube;
	else{if(0<=maxpow){
	if(sh<.900001)
	res += (-0.8559670781893004*L - 
	0.2962962962962963*sh*pow(ash,2)*sh_m1_pow_m4 - 
	0.00823045267489712*Lsh*
	(29. - 18.*Lshb*(-1. + sh) - 47.*sh)*sh*
	sh_m1_pow_m2 + 0.14814814814814814*Li2sh*sh*
	sh_m1_pow_m1 - 0.07407407407407407*pow(Lsh,2)*
	sh_m1_pow_m3*sh_cube - 
	0.01646090534979424*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m1*(-4. + 9.*sh - 15.*sh_sq + 
	4.*sh_cube) - 0.0013717421124828531*
	sh_m1_pow_m4*(785. - 
	2.*sh*(1585. + 12.*pi_sq) + 
	6.*(803. + 9.*pi_sq)*sh_sq - 
	2.*(1633. + 27.*pi_sq)*
	sh_cube + (833. + 18.*pi_sq)*
	sh_four));
	else
	res += (-0.7509911973035993 - 0.21726138504882622*(1. - sh) - 
	0.8559670781893004*L - 
	0.0659657518981201*msh_p1_pow_2 - 
	0.03129640614981805*msh_p1_pow_3 - 
	0.017719628843137502*msh_p1_pow_4);}

	if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}

	if(2<=maxpow){
	if(sh<.900001)
	res += (0.04938271604938271*sh*pow(Lsh,4)*sh_m1_pow_m3 - 
	0.2962962962962963*Li2sh*(-1. + 4.*sh)*
	sh_m1_pow_m2 + 0.09876543209876543*(-3. + sh)*
	pow(Lsh,3)*sh_m1_pow_m2 - 
	1.*lz*(-0.2962962962962963*Lsh*
	(2. + 2.*Lshb*(-1. + sh) + sh)*sh_m1_pow_m2 - 
	0.5925925925925926*Li2sh*sh_m1_pow_m1 + 
	0.09876543209876543*(9. + pi_sq)*
	sh_m1_pow_m1 + 
	0.2962962962962963*pow(Lsh,2)*sh_m1_pow_m1) - 
	0.35616500834358344*sh_m1_pow_m3*
	(-5. + 5.*sh + 2.*sh_sq) - 
	0.2962962962962963*Li3sh*sh_m1_pow_m3*
	(-1. - 5.*sh + 4.*sh_sq) - 
	1.*pow(Lsh,2)*(-0.5925925925925926*Li2sh*sh*
	sh_m1_pow_m3 - 
	0.04938271604938271*sh_m1_pow_m3*
	(12. + 3.*Lshb*(-5. + 3.*sh) + 
	sh*(-15. + 2.*pi_sq) + 
	6.*sh_sq)) - 
	1.*Lsh*(-2.8493200667486676*sh*sh_m1_pow_m3 + 
	2.3703703703703702*Li3sh*sh*sh_m1_pow_m3 + 
	0.2962962962962963*Li2sh*sh_m1_pow_m3*
	(3. + sh - 2.*sh_sq) + 
	0.04938271604938271*sh_m1_pow_m3*
	(12. + 21.*sh + 
	(-5. + 3.*sh)*pi_sq - 
	33.*sh_sq + 
	6.*Lshb*(1. - 5.*sh + 4.*sh_sq))) - 
	1.7777777777777777*pow(ash,2)*sh_m1_pow_m4*
	(-1. - 3.*sh_sq + sh_cube) + 
	0.009876543209876543*sh_m1_pow_m4*
	(-4.*(-1. + sh)*sh*pow(pi,4.) - 
	15.*(-1. + sh)*(10. - (17. + 24.*Li4sh)*sh + 
	7.*sh_sq) + 
	5.*pi_sq*
	(-6. + 16.*sh - 20.*sh_sq + 7.*sh_cube)) - 
	0.8888888888888888*ash*(1. + sh)*sh_m1_pow_m3*
	sh_m4_sh_sqrt)*z;
	else
	res += (-1.5156467116109078 - 0.06864318910496076*(1. - sh) + 
	0.19606872711735043*msh_p1_pow_2 + 
	0.3270513130738385*msh_p1_pow_3 + 
	log(z)*(0.4444444444444444 + 
	0.2962962962962963*(1. - sh) + 
	0.23868312757201646*msh_p1_pow_2 + 
	0.20493827160493827*msh_p1_pow_3 + 
	0.1817283950617284*msh_p1_pow_4) + 
	0.38673058027546225*msh_p1_pow_4)*z;}

	if(3<=maxpow){
	if(sh<.900001)
	res += 0.19753086419753085*(2. + sh)*pi_sq*
	sh_m1*pow(z,1.5);
	else
	res += (5.848654459904805 + 3.8991029732698697*(1. - sh) + 
	3.8991029732698697*msh_p1_pow_2 + 
	3.8991029732698697*msh_p1_pow_3 + 
	3.8991029732698697*msh_p1_pow_4)*pow(z,1.5);}

	if(4<=maxpow){
	if(sh<.900001)
	res += (-0.2962962962962963*Li3sh*(5. + 3.*sh)*sh_m1_pow_m3 + 
	0.2962962962962963*Li3shb*(1. + 6.*sh)*
	sh_m1_pow_m3 + 0.14814814814814814*pow(Lsh,3)*
	sh_m1_pow_m3 - 0.8888888888888888*ash*(-3. + sh)*
	sqrt4sh*sqrtsh*sh_m1_pow_m3*sh_sq + 
	0.35616500834358344*sh_m1_pow_m4*
	(-5. + 8.*sh + 3.*sh_sq) + 
	0.2962962962962963*Li2sh*sh_m1_pow_m3*sh_m1*
	(3. - 10.*sh + 8.*sh_sq) - 
	0.8888888888888888*pow(ash,2)*sh_m1_pow_m4*
	sh_sq*(2. + 9.*sh - 6.*sh_sq + sh_cube) - 
	0.07407407407407407*pow(Lsh,2)*sh_m1_pow_m4*
	sh_m2*(-4. + 7.*sh + (5. + 6.*Lshb)*sh_sq - 
	1.*(25. + 8.*Lshb)*sh_cube + 
	(11. + 2.*Lshb)*sh_four) - 
	1.*lz_sq*(0.14814814814814814*Lsh*(1. + 6.*sh)*
	sh_m1_pow_m3 - 
	0.07407407407407407*sh_m1_pow_m2*sh_m2*
	(1. - 2.*sh + 8.*sh_sq + 10.*sh_cube - 
	4.*sh_four + sh_five)) - 
	1.*lz*(-10.666666666666666*sh*pow(ash,2)*
	sh_m1_pow_m4 - 
	0.2962962962962963*Li2sh*(1. + 6.*sh)*
	sh_m1_pow_m3 - 
	0.14814814814814814*pow(Lsh,2)*sh_m1_pow_m2 - 
	0.14814814814814814*Lsh*sh_m1_pow_m3*sh_m2*
	(2. - 3.*sh + 2.*(1. + Lshb)*sh_sq + 
	3.*(-1. + 4.*Lshb)*sh_cube) + 
	0.024691358024691357*sh_m1_pow_m4*sh_m2*
	(9. - 42.*sh + (81. - 2.*pi_sq)*
	sh_sq + 2.*(-6. + pi_sq)*
	sh_cube + 6.*
	(-23. + 2.*pi_sq)*sh_four\
	+ 153.*sh_five - 60.*sh_six + 9.*sh_seven)
	 - 0.5925925925925926*ash*sh_m1_pow_m3*
	(-4. - 3.*sh_sq + sh_cube)*
	sh_m4_sh_sqrt) - 
	1.*Lsh*(-5.333333333333333*sh*pow(ash,2)*
	sh_m1_pow_m4 - 
	0.2962962962962963*Li2sh*(4. + sh)*
	sh_m1_pow_m3 - 
	0.024691358024691357*sh_m1_pow_m3*sh_m2*
	(-18. + (57. + 36.*Lshb)*sh - 
	6.*(5. + 20.*Lshb + pi_sq)*
	sh_sq + (-93. + 96.*Lshb + 
	50.*pi_sq)*sh_cube) - 
	0.2962962962962963*ash*sh_m1_pow_m3*
	(-4. - 3.*sh_sq + sh_cube)*
	sh_m4_sh_sqrt) + 
	0.012345679012345678*sh_m1_pow_m4*sh_m2*
	(2.*pi_sq*
	(-6. + 45.*sh - 133.*sh_sq + 157.*sh_cube - 
	59.*sh_four + 12.*sh_five - 12.*sh_six + 
	2.*sh_seven) + 
	3.*(7. - 78.*sh + 7.*sh_seven + 
	sh_five*(155. - 
	32.*cl2*sh_m4_sh_sqrt) + 
	8.*sh_six*
	(-7. + cl2*sh_m4_sh_sqrt) - 
	8.*sh_cube*(27. + 18.*cl3 + 
	4.*cl2*sh_m4_sh_sqrt) + 
	sh_four*(-62. + 
	24.*cl2*sh_m4_sh_sqrt) + 
	sh_sq*(243. + 
	32.*cl2*sh_m4_sh_sqrt))))*z_sq;
	else
	res += z_sq*(-0.647637155328409 - 
	0.9413088193612421*(1. - sh) - 
	0.8385255835044603*msh_p1_pow_2 - 
	1.1019705279460212*msh_p1_pow_3 - 
	1.6181578946132955*msh_p1_pow_4 + 
	log(z)*(1.7526351083839322 + 
	0.5609956120680696*(1. - sh) + 
	1.1754785561935277*msh_p1_pow_2 + 
	1.6549606263554564*msh_p1_pow_3 + 
	2.1141011282028064*msh_p1_pow_4) + 
	(0.09876543209876543 + 
	0.1111111111111111*(1. - sh) + 
	0.23703703703703705*msh_p1_pow_2 + 
	0.3012345679012346*msh_p1_pow_3 + 
	0.37037037037037035*msh_p1_pow_4)*pow(log(z),2)
	);}

	if(5<=maxpow){
	if(sh<.900001)
	res += 0.02633744855967078*pi_sq*sh_m3*
	(3. + 14.*sh + 3.*sh_sq)*pow(z,2.5);
	else
	res += (5.1988039643598265 + 10.397607928719653*(1. - sh) + 
	16.376232487733454*msh_p1_pow_2 + 
	23.13467764140123*msh_p1_pow_3 + 
	30.672943389722974*msh_p1_pow_4)*pow(z,2.5);}

	if(6<=maxpow){
	if(sh<.900001)
	res += (1.1851851851851851*Li3shb*sh*sh_m1_pow_m4 - 
	0.5925925925925926*Li3sh*sh_m1_pow_m3*sh_m1 + 
	0.19753086419753085*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m1 - 0.7123300166871669*sh_m1_pow_m4*
	sh_m1*(1. - sh + 2.*sh_sq) - 
	0.2962962962962963*Li2sh*sh_m1_pow_m4*sh_m2*
	(4. - 14.*sh + 18.*sh_sq - 11.*sh_cube + 
	sh_four) - 0.19753086419753085*cl2*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-9. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) - 0.06584362139917696*ash*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(27. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) + 0.5925925925925926*sh*pow(ash,2)*
	sh_m1_pow_m4*(-6. + 27.*sh - 30.*sh_cube + 
	27.*sh_four - 9.*sh_five + sh_six) - 
	1.*Lsh*(3.5555555555555554*sh*pow(ash,2)*
	sh_m1_pow_m4 - 
	0.5925925925925926*Li2sh*sh_m1_pow_m3*
	sh_m1 + 0.19753086419753085*ash*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-9. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) + 0.00823045267489712*sh_m1_pow_m4*
	sh_m4*(-10. + 10.*sh + 
	6.*(23. + 24.*Lshb)*sh_sq - 
	2.*(67. + 252.*Lshb + 
	30.*pi_sq)*sh_cube + 
	(-163. + 648.*Lshb + 
	60.*pi_sq)*sh_four - 
	12.*(-5. + 33.*Lshb + 
	6.*pi_sq)*sh_five + 
	9.*(-5. + 4.*Lshb)*sh_six)) - 
	0.012345679012345678*pow(Lsh,2)*sh_m1_pow_m5*
	sh_m4*(8. + 8.*sh - 94.*sh_sq + 
	(94. - 24.*Lshb)*sh_cube + 
	(52. + 48.*Lshb)*sh_four - 
	8.*(19. + 3.*Lshb)*sh_five + 72.*sh_six - 
	12.*sh_seven) - 
	1.*lz_sq*(0.5925925925925926*Lsh*sh*
	sh_m1_pow_m4 + 
	0.024691358024691357*sh_m1_pow_m3*sh_m4*
	(1. + 3.*sh - 18.*sh_sq + 25.*sh_cube - 
	36.*sh_four + 26.*sh_five - 27.*sh_six - 
	22.*sh_seven + 38.*sh_eight - 16.*sh_nine + 
	2.*sh_ten)) - 
	1.*lz*(-1.1851851851851851*Li2sh*sh*sh_m1_pow_m4 + 
	7.111111111111111*sh*pow(ash,2)*sh_m1_pow_m4 + 
	0.2962962962962963*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m1 + 0.3950617283950617*ash*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-9. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) + 0.04938271604938271*Lsh*
	sh_m1_pow_m4*sh_m4*
	(2. + 4.*sh - 27.*sh_sq + 22.*sh_cube + 
	14.*sh_four - 12.*(1. + 2.*Lshb)*sh_five + 
	9.*sh_six) + 
	0.00823045267489712*sh_m1_pow_m3*sh_m4*
	(-5. + 57.*sh_sq - 46.*sh_cube - 
	85.*sh_four + 164.*sh_five - 73.*sh_six + 
	102.*sh_seven - 40.*sh_eight - 4.*sh_nine + 
	2.*sh_ten)) - 
	0.0013717421124828531*sh_m1_pow_m4*sh_m4*
	(-19. + 46.*sh + 1365.*sh_sq - 4435.*sh_cube + 
	3757.*sh_four + (610. - 2592.*cl3)*sh_five - 
	3103.*sh_six + 2113.*sh_seven + 
	594.*sh_eight - 1548.*sh_nine + 
	714.*sh_ten - 94.*sh_eleven + 
	6.*pi_sq*
	(6. + 12.*sh - 195.*sh_sq + 442.*sh_cube - 
	434.*sh_four + 124.*sh_five + 
	137.*sh_six - 16.*sh_seven - 
	120.*sh_eight + 108.*sh_nine - 
	36.*sh_ten + 4.*sh_eleven)))*z_cube;
	else
	res += pow(z,3.)*(2.0558846293203854 - 
	5.023381570340977*(1. - sh) - 
	7.761640024983791*msh_p1_pow_2 - 
	13.35266138630302*msh_p1_pow_3 - 
	24.65891981493292*msh_p1_pow_4 + 
	log(z)*(2.772759928702021 - 
	0.24604278849693273*(1. - sh) + 
	3.842007139630303*msh_p1_pow_2 + 
	8.253815580345757*msh_p1_pow_3 + 
	13.162562518769866*msh_p1_pow_4) + 
	(-0.2222222222222222 + 
	0.044444444444444446*(1. - sh) + 
	0.9679012345679012*msh_p1_pow_2 + 
	1.689594356261023*msh_p1_pow_3 + 
	2.6807760141093473*msh_p1_pow_4)*pow(log(z),2));}

	if(7<=maxpow){
	if(sh<.900001)
	res += 0.002257495590828924*pi_sq*sh_m5*
	(15. + 108.*sh + 314.*sh_sq + 108.*sh_cube + 
	15.*sh_four)*pow(z,3.5);
	else
	res += (12.477129514463583 + 37.43138854339075*(1. - sh) + 
	78.60591594112057*msh_p1_pow_2 + 
	139.74385056199213*msh_p1_pow_3 + 
	224.92254008662476*msh_p1_pow_4)*pow(z,3.5);}

	if(8<=maxpow){
	if(sh<.900001)
	res += (1.1872166944786116*sh_m1_pow_m3*sh_m2 - 
	0.9876543209876543*Li3sh*sh_m1_pow_m3*sh_m2 + 
	0.3292181069958848*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m2 + 0.2962962962962963*cl2*sh*sqrt4sh*sqrtsh*
	(8. + 6.*sh - 6.*sh_sq + sh_cube) - 
	0.04938271604938271*Li2sh*sh_m1_pow_m5*sh_m3*
	(-50. + 230.*sh - 420.*sh_sq + 380.*sh_cube - 
	146.*sh_four + 15.*sh_five) + 
	0.04938271604938271*ash*sh*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*(16. + 72.*sh + 18.*sh_sq - 
	203.*sh_cube + 189.*sh_four - 63.*sh_five + 
	7.*sh_six) - 0.8888888888888888*pow(ash,2)*
	sh_m1_pow_m4*sh_sq*
	(36. - 60.*sh + 99.*sh_cube - 112.*sh_four + 
	54.*sh_five - 12.*sh_six + sh_seven) - 
	1.*Lsh*(-0.9876543209876543*Li2sh*sh_m1_pow_m3*
	sh_m2 - 0.2962962962962963*ash*sh*sqrt4sh*
	sqrtsh*(8. + 6.*sh - 6.*sh_sq + sh_cube) + 
	0.0013717421124828531*sh_m1_pow_m5*sh_m6*
	(21. - 10.*sh - 307.*sh_sq - 
	1.*(283. + 1800.*Lshb)*sh_cube + 
	5.*(437. + 1656.*Lshb + 
	120.*pi_sq)*sh_four - 
	2.*(1081. + 7560.*Lshb + 
	600.*pi_sq)*sh_five + 
	(253. + 13680.*Lshb + 
	600.*pi_sq)*sh_six + 
	(861. - 5256.*Lshb)*sh_seven + 
	6.*(49. + 90.*Lshb)*sh_eight - 78.*sh_nine))\
	+ 0.00823045267489712*pow(Lsh,2)*sh_m1_pow_m6*
	sh_m6*(6. + 14.*sh - 30.*sh_sq - 
	255.*sh_cube + (729. - 60.*Lshb)*sh_four + 
	45.*(-13. + 4.*Lshb)*sh_five - 
	15.*(11. + 12.*Lshb)*sh_six + 
	(541. + 60.*Lshb)*sh_seven - 270.*sh_eight + 
	45.*sh_nine) + 
	0.00411522633744856*lz_sq*sh_m6*
	(3. + 25.*sh + 90.*sh_sq - 60.*sh_cube + 
	5.*sh_four - 3.*sh_five - 144.*sh_six - 
	252.*sh_seven + 72.*sh_eight + 288.*sh_nine - 
	144.*sh_ten + 18.*sh_eleven) - 
	1.*lz*(0.49382716049382713*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 - 0.5925925925925926*ash*sh*sqrt4sh*
	sqrtsh*(8. + 6.*sh - 6.*sh_sq + sh_cube) + 
	0.01646090534979424*Lsh*sh_m1_pow_m5*sh_m6*
	(-3. - 10.*sh + 5.*sh_sq + 185.*sh_cube - 
	425.*sh_four + 343.*sh_five - 
	125.*sh_six + 3.*sh_seven - 15.*sh_eight + 
	15.*sh_nine) + 
	0.00013717421124828533*sh_m1_pow_m4*sh_m6*
	(105. + 55.*sh - 1480.*sh_sq + 90.*sh_cube + 
	10470.*sh_four - 25782.*sh_five + 
	43242.*sh_six - 52678.*sh_seven + 
	52987.*sh_eight - 14283.*sh_nine - 
	23976.*sh_ten + 3240.*sh_eleven + 
	28080.*sh_twelve - 22680.*sh_thirteen + 
	6480.*sh_fourteen - 630.*sh_fifteen)) + 
	2.2862368541380887e-6*sh_m1_pow_m5*sh_m6*
	(-2775. + 5550.*sh + 34425.*sh_sq + 
	1.48905e6*sh_cube - 7.4002e6*sh_four + 
	1.4083748e7*sh_five - 1.4200216e7*sh_six + 
	9.34426e6*sh_seven - 5.568515e6*sh_eight + 
	141850.*sh_nine + 4.948647e6*sh_ten - 
	27324.*sh_eleven - 8.9613e6*sh_twelve + 
	9.8388e6*sh_thirteen - 4.6764e6*sh_fourteen + 
	1.03725e6*sh_fifteen - 86850.*sh_sixteen + 
	1200.*pi_sq*
	(9. + 30.*sh - 15.*sh_sq - 1245.*sh_cube + 
	4325.*sh_four - 6463.*sh_five + 
	5357.*sh_six - 2657.*sh_seven + 
	295.*sh_eight + 1405.*sh_nine - 
	1014.*sh_ten - 1782.*sh_eleven + 
	3798.*sh_twelve - 2988.*sh_thirteen + 
	1188.*sh_fourteen - 234.*sh_fifteen + 
	18.*sh_sixteen)))*pow(z,4.);
	else
	res += pow(z,4.)*(4.166076007573638 - 
	29.482916540765096*(1. - sh) - 
	51.76439651455935*msh_p1_pow_2 - 
	100.58418925538838*msh_p1_pow_3 - 
	217.9981110758939*msh_p1_pow_4 + 
	log(z)*(6.036523956202085 - 
	6.188787330687681*(1. - sh) + 
	14.248013358351397*msh_p1_pow_2 + 
	42.62670356928015*msh_p1_pow_3 + 
	77.13858012891833*msh_p1_pow_4) + 
	(-0.41975308641975306 + 
	0.24691358024691357*(1. - sh) + 
	5.111111111111111*msh_p1_pow_2 + 
	9.74485596707819*msh_p1_pow_3 + 
	17.88477366255144*msh_p1_pow_4)*pow(log(z),2));}

	if(9<=maxpow){
	if(sh<.900001)
	res += 0.0005374989501973629*pi_sq*sh_m7*
	(35. + 330.*sh + 1389.*sh_sq + 3212.*sh_cube + 
	1389.*sh_four + 330.*sh_five + 35.*sh_six)*
	pow(z,4.5);
	else
	res += (35.64894146989595 + 142.5957658795838*(1. - sh) + 
	372.5314383604127*msh_p1_pow_2 + 
	793.1889477051849*msh_p1_pow_3 + 
	1491.2079534505585*msh_p1_pow_4)*pow(z,4.5);}

	if(10<=maxpow){
	if(sh<.900001)
	res += (2.493155058405084*sh_m1_pow_m3*sh_m3 - 
	2.074074074074074*Li3sh*sh_m1_pow_m3*sh_m3 + 
	0.691358024691358*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m3 - 0.5925925925925926*cl2*sqrt4sh*sqrtsh*
	sh_sq*(-20. - 11.*sh + 24.*sh_sq - 
	9.*sh_cube + sh_four) - 
	0.019753086419753086*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_sq*
	(630. - 1853.*sh - 99.*sh_sq + 4458.*sh_cube - 
	5217.*sh_four + 2538.*sh_five - 
	564.*sh_six + 47.*sh_seven) + 
	0.09876543209876543*Li2sh*sh_m1_pow_m6*sh_m4*
	(-63. + 355.*sh - 826.*sh_sq + 1008.*sh_cube - 
	665.*sh_four + 213.*sh_five - 7.*sh_seven + 
	sh_eight) + 1.7777777777777777*pow(ash,2)*
	sh_m1_pow_m4*sh_sq*
	(10. - 90.*sh + 165.*sh_sq - 333.*sh_four + 
	445.*sh_five - 275.*sh_six + 90.*sh_seven - 
	15.*sh_eight + sh_nine) - 
	1.*Lsh*(-2.074074074074074*Li2sh*sh_m1_pow_m3*
	sh_m3 + 0.5925925925925926*ash*sqrt4sh*sqrtsh*
	sh_sq*(-20. - 11.*sh + 24.*sh_sq - 
	9.*sh_cube + sh_four) + 
	0.0000823045267489712*sh_m1_pow_m6*sh_m8*
	(-162. + 36.*sh + 2914.*sh_sq + 
	5862.*sh_cube + 
	30.*(-767. + 2520.*Lshb)*sh_four - 
	8.*(2131. + 53250.*Lshb + 
	2625.*pi_sq)*sh_five + 
	8.*(12113. + 123900.*Lshb + 
	7875.*pi_sq)*sh_six - 
	18.*(5693. + 67200.*Lshb + 
	3500.*pi_sq)*sh_seven + 
	3.*(15101. + 266000.*Lshb + 
	7000.*pi_sq)*sh_eight - 
	50.*(259. + 5112.*Lshb)*sh_nine + 
	9600.*sh_ten + 
	150.*(11. + 56.*Lshb)*sh_eleven - 
	25.*(25. + 48.*Lshb)*sh_twelve)) - 
	0.0024691358024691358*pow(Lsh,2)*sh_m1_pow_m7*
	sh_m8*(12. + 42.*sh - 28.*sh_sq - 
	308.*sh_cube - 1533.*sh_four + 
	(8318. - 420.*Lshb)*sh_five + 
	6.*(-2237. + 280.*Lshb)*sh_six - 
	56.*(-143. + 45.*Lshb)*sh_seven + 
	12.*(149. + 140.*Lshb)*sh_eight - 
	1.*(4507. + 420.*Lshb)*sh_nine + 
	1540.*sh_ten + 140.*sh_eleven - 
	160.*sh_twelve + 20.*sh_thirteen) - 
	0.0012345679012345679*lz_sq*sh_m8*
	(-6. - 63.*sh - 301.*sh_sq - 840.*sh_cube + 
	630.*sh_four - 99.*sh_five + 153.*sh_six + 
	6.*sh_seven + 1200.*sh_eight + 1440.*sh_nine + 
	3360.*sh_ten - 2400.*sh_eleven - 
	5400.*sh_twelve + 4800.*sh_thirteen - 
	1320.*sh_fourteen + 120.*sh_fifteen) - 
	1.*lz*(1.037037037037037*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 + 1.1851851851851851*ash*sqrt4sh*sqrtsh*
	sh_sq*(-20. - 11.*sh + 24.*sh_sq - 
	9.*sh_cube + sh_four) + 
	0.0049382716049382715*Lsh*sh_m1_pow_m6*
	sh_m8*(6. + 27.*sh + 13.*sh_sq - 
	141.*sh_cube - 1380.*sh_four + 
	5549.*sh_five - 8207.*sh_six + 
	6207.*sh_seven - 2809.*sh_eight + 
	1080.*sh_nine - 30.*sh_ten - 
	30.*sh_eleven + 35.*sh_twelve) + 
	0.000011757789535567314*sh_m1_pow_m5*sh_m8*
	(-567. - 441.*sh + 9758.*sh_sq + 
	30275.*sh_cube - 166495.*sh_four + 
	20923.*sh_five + 866054.*sh_six - 
	1.788427e6*sh_seven + 1.255755e6*sh_eight + 
	458890.*sh_nine - 1.258076e6*sh_ten + 
	55832.*sh_eleven - 524481.*sh_twelve + 
	2.826e6*sh_thirteen - 901320.*sh_fourteen - 
	5.25168e6*sh_fifteen + 7.8771e6*sh_sixteen - 
	5.0799e6*sh_seventeen + 1.7073e6*sh_eighteen - 
	290640.*sh_nineteen + 19740.*sh_twenty)) - 
	1.3997368494722993e-8*sh_m1_pow_m6*sh_m8*
	(-161406. + 318843.*sh + 1.910657e6*sh_sq - 
	738969.*sh_cube + 6.33135615e8*sh_four - 
	3.654733992e9*sh_five + 8.580730746e9*sh_six - 
	1.0387902666e10*sh_seven + 
	6.272402572e9*sh_eight - 4.35513965e8*sh_nine - 
	2.310469467e9*sh_ten + 
	2.213537811e9*sh_eleven - 
	4.191847021e9*sh_twelve + 
	8.659496442e9*sh_thirteen - 
	5.00515056e9*sh_fourteen - 
	1.057709688e10*sh_fifteen + 
	2.294859924e10*sh_sixteen - 
	2.03364504e10*sh_seventeen + 
	1.00250472e10*sh_eighteen - 
	2.83732932e9*sh_nineteen + 4.2896364e8*sh_twenty - 
	2.674812e7*sh_twenty_one + 
	58800.*pi_sq*
	(18. + 81.*sh + 39.*sh_sq - 423.*sh_cube - 
	9810.*sh_four + 47729.*sh_five - 
	95467.*sh_six + 103887.*sh_seven - 
	63289.*sh_eight + 15730.*sh_nine + 
	8388.*sh_ten - 21404.*sh_eleven + 
	46029.*sh_twelve - 50988.*sh_thirteen - 
	20160.*sh_fourteen + 133320.*sh_fifteen - 
	179760.*sh_sixteen + 130200.*sh_seventeen - 
	56400.*sh_eighteen + 14520.*sh_nineteen - 
	2040.*sh_twenty + 120.*sh_twenty_one)))*pow(z,5.);
	else
	res += pow(z,5.)*(9.495992115022151 - 
	144.07780670229167*(1. - sh) - 
	287.1615196453904*msh_p1_pow_2 - 
	629.2386817194127*msh_p1_pow_3 - 
	1577.85840284407*msh_p1_pow_4 + 
	log(z)*(18.273866188842664 - 
	36.50758667966486*(1. - sh) + 
	57.239414092843454*msh_p1_pow_2 + 
	217.34199886507017*msh_p1_pow_3 + 
	420.0317629604588*msh_p1_pow_4) + 
	(-1.5802469135802468 + 
	0.2962962962962963*(1. - sh) + 
	25.*msh_p1_pow_2 + 
	51.65432098765432*msh_p1_pow_3 + 
	108.14814814814815*msh_p1_pow_4)*pow(log(z),2));}

	}	
	return res;
}

#ifdef CACHE
double F_17re_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
			
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else {currentmax=0;loop++;}
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_17re(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_17im(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double lz_sq = pow(lz,2);
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;     
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = -0.517135 + 1.39626*sh + 0.465421*lz*sh + 1.59186*z + 
	0.93084*lz*z + 5.78065*sh*z + 0.930842*Lsh*sh*z + 
	6.5159*lz*sh*z + 0.9308425*z*lz_sq + 
	0.9308425*sh*z*lz_sq + 1.59019*sh_sq + 
	0.465421*lz*sh_sq + 14.7788*z*sh_sq + 
	1.86168*Lsh*z*sh_sq + 11.1701*lz*z*sh_sq + 
	0.9308425*z*lz_sq*sh_sq + 1.63673*sh_cube + 
	0.465421*lz*sh_cube + 25.8196*z*sh_cube + 
	2.79253*Lsh*z*sh_cube + 15.51405*lz*z*sh_cube + 
	0.9308425*z*lz_sq*sh_cube + 
	0.00258567*sh_cube*z_m2 - 
	0.0258567*sh_cube*z_m1 - 3.06235*z_sq - 
	15.2375*sh*z_sq + 5.58505*lz*sh*z_sq + 
	3.72337*Lsh*lz*sh*z_sq + 
	0.9308425*lz_sq*z_sq + 
	4.6542*sh*lz_sq*z_sq - 
	32.8021*sh_sq*z_sq + 
	3.72337*Lsh*sh_sq*z_sq + 
	19.5477*lz*sh_sq*z_sq + 
	11.1701*Lsh*lz*sh_sq*z_sq + 
	11.1701*lz_sq*sh_sq*z_sq - 
	52.4207*sh_cube*z_sq + 
	13.0318*Lsh*sh_cube*z_sq + 
	44.37015*lz*sh_cube*z_sq + 
	22.3402*Lsh*lz*sh_cube*z_sq + 
	20.478525*lz_sq*sh_cube*z_sq + 
	2.89595*z_cube - 1.241125*lz*z_cube + 
	15.6984*sh*z_cube + 1.86168*Lsh*sh*z_cube - 
	12.41125*lz*sh*z_cube - 3.72337*Lsh*lz*sh*z_cube + 
	34.8683*sh_sq*z_cube + 
	3.72337*Lsh*sh_sq*z_cube - 
	40.3365*lz*sh_sq*z_cube - 
	14.8935*Lsh*lz*sh_sq*z_cube + 
	54.2499*sh_cube*z_cube + 
	1.86168*Lsh*sh_cube*z_cube - 
	92.4635*lz*sh_cube*z_cube - 
	37.2337*Lsh*lz*sh_cube*z_cube;
	else{
		
	if(0<=maxpow){
	if(sh<.900001)
	res += (0.4654211338651545*Lsh*sh*(-1. + 2.*sh)*sh_m1_pow_m3 - 
	0.025856729659175254*sh_m1_pow_m2*
	(20. - 49.*sh + 47.*sh_sq));
	else
	res += (-0.8274153490936081 + 0.11635528346628862*(1. - sh) + 
	0.054299132284268026*msh_p1_pow_2 + 
	0.031028075591010302*msh_p1_pow_3 + 
	0.019946620022792336*msh_p1_pow_4);}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (-8.951402989423594*sh*sh_m1_pow_m3 + 
	7.446738141842472*Li3sh*sh*sh_m1_pow_m3 - 
	0.930842267730309*Li2sh*(-5. + 3.*sh)*
	sh_m1_pow_m3 - 0.620561511820206*sh*pow(Lsh,3)*
	sh_m1_pow_m3 + 1.861684535460618*pow(Lsh,2)*
	sh_m1_pow_m2 - Lsh*
	(3.723369070921236*Li2sh*sh*sh_m1_pow_m3 + 
	0.930842267730309*sh_m1_pow_m3*
	(3. + Lshb*(-5. + 3.*sh) + 
	2.*sh*pi_sq - 2.*sh_sq))\
	+ 0.1551403779550515*sh_m1_pow_m3*
	(-6.*(2. - 5.*sh + 3.*sh_sq) + 
	pi_sq*
	(-9. + 3.*sh + 4.*sh_sq)))*z;
	else
	res += (4.1876741044473365 + 2.7402552937030498*(1. - sh) + 
	2.1232783768063475*msh_p1_pow_2 + 
	1.7630333644902134*msh_p1_pow_3 + 
	1.5217395214542553*msh_p1_pow_4)*z;}

	/*if(3<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(4<=maxpow){
	if(sh<.900001)
	res += (-0.930842267730309*Li2sh*(2. - 7.*sh)*sh_m1_pow_m3 + 
	0.930842267730309*pow(Lsh,2)*sh_m1_pow_m2 - 
	1.*lz*(6.5158958741121635*Lsh*sh*sh_m1_pow_m3 - 
	0.4654211338651545*sh_m1_pow_m2*sh_m2*
	(2. - 7.*sh + 15.*sh_sq + 4.*sh_cube)) - 
	0.07757018897752575*sh_m1_pow_m3*sh_m2*
	(-18. + 63.*sh - 2.*
	(63. + 2.*pi_sq)*sh_sq + 
	(9. + 14.*pi_sq)*sh_cube + 
	36.*sh_four) + 
	0.4654211338651545*Lsh*sh_m1_pow_m4*sh_m2*
	(-4. + 13.*sh + (-21. + 4.*Lshb)*sh_sq + 
	(11. - 18.*Lshb)*sh_cube + 
	(-5. + 14.*Lshb)*sh_four))*z_sq;
	else
	res += (1.52554704989134 + 3.0730723199929786*(1. - sh) + 
	4.914071471726257*msh_p1_pow_2 + 
	7.194598089480287*msh_p1_pow_3 + 
	9.866611424925042*msh_p1_pow_4 + 
	log(z)*(0.620561511820206 + 
	1.0084124567078347*(1. - sh) + 
	1.7220581953010718*msh_p1_pow_2 + 
	2.544302198462845*msh_p1_pow_3 + 
	3.413088315011133*msh_p1_pow_4))*z_sq;}

	if(5<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}

	if(6<=maxpow){
	if(sh<.900001)
	res += (-1.861684535460618*pow(Lsh,2)*sh_m1_pow_m3*sh_m1 + 
	1.861684535460618*Li2sh*sh_m1_pow_m4*sh_m1*
	(1. - sh + 2.*sh_sq) + 
	0.1551403779550515*Lsh*sh_m1_pow_m5*sh_m4*
	(4. + 4.*sh - 71.*sh_sq + 
	(155. - 12.*Lshb)*sh_cube + 
	2.*(-83. + 12.*Lshb)*sh_four + 
	(98. - 36.*Lshb)*sh_five + 
	12.*(-3. + 2.*Lshb)*sh_six) - 
	1.*lz*(1.861684535460618*Lsh*sh_m1_pow_m4*
	sh_m1*(1. - sh + 2.*sh_sq) + 
	0.1551403779550515*sh_m1_pow_m3*sh_m4*
	(2. + 6.*sh - 45.*sh_sq + 61.*sh_cube - 
	41.*sh_four - 11.*sh_five + 4.*sh_six)) - 
	0.025856729659175254*sh_m1_pow_m4*sh_m4*
	(10. - 10.*sh - 54.*sh_sq + 
	2.*(-23. + 6.*pi_sq)*sh_cube - 
	3.*(-53. + 4.*pi_sq)*sh_four + 
	24.*(-10. + pi_sq)*sh_five + 
	141.*sh_six - 32.*sh_seven))*z_cube;
	else
	res += (2.32969134229169 + 7.64686922940449*(1. - sh) + 
	18.94249239557088*msh_p1_pow_2 + 
	38.204310125548965*msh_p1_pow_3 + 
	67.3328659256586*msh_p1_pow_4 + 
	log(z)*(0.1551403779550515 + 
	3.07177948351002*(1. - sh) + 
	8.625805014300864*msh_p1_pow_2 + 
	17.193986459646997*msh_p1_pow_3 + 
	29.117632651049522*msh_p1_pow_4))*z_cube;}

	/*if(7<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(8<=maxpow){
	if(sh<.900001)
	res += (-3.1028075591010302*Li2sh*sh_m1_pow_m3*sh_m2 - 
	3.1028075591010302*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 - lz*
	(-3.1028075591010302*Lsh*sh_m1_pow_m3*sh_m2 - 
	0.05171345931835051*sh_m1_pow_m2*sh_m6*
	(3. + 19.*sh + 43.*sh_sq - 260.*sh_cube + 
	225.*sh_four - 105.*sh_five + 15.*sh_six))
	 - 0.05171345931835051*Lsh*sh_m1_pow_m6*sh_m6*
	(6. + 14.*sh - 30.*sh_sq - 405.*sh_cube + 
	(1569. - 60.*Lshb)*sh_four + 
	15.*(-169. + 12.*Lshb)*sh_five + 
	(2235. - 180.*Lshb)*sh_six + 
	(-1037. + 60.*Lshb)*sh_seven + 213.*sh_eight) + 
	0.004309454943195875*sh_m1_pow_m5*sh_m6*
	(21. - 10.*sh - 307.*sh_sq + 1247.*sh_cube + 
	(-3143. + 120.*pi_sq)*
	sh_four + (3526. - 
	240.*pi_sq)*sh_five + 
	(397. + 120.*pi_sq)*sh_six - 
	2577.*sh_seven + 1140.*sh_eight + 
	138.*sh_nine - 72.*sh_ten))*pow(z,4.);
	else
	res += (2.999380640464329 + 22.55925786392586*(1. - sh) + 
	78.72635417514317*msh_p1_pow_2 + 
	200.15639021889723*msh_p1_pow_3 + 
	424.28359587096105*msh_p1_pow_4 + 
	log(z)*(2.482246047280824 + 
	14.11777439390969*(1. - sh) + 
	44.215007717189685*msh_p1_pow_2 + 
	104.04748014852122*msh_p1_pow_3 + 
	206.9868147402206*msh_p1_pow_4))*pow(z,4.);}

	/*if(9<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(10<=maxpow){
	if(sh<.900001)
	res += (-6.5158958741121635*Li2sh*sh_m1_pow_m3*sh_m3 - 
	6.5158958741121635*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 - lz*
	(-6.5158958741121635*Lsh*sh_m1_pow_m3*sh_m3 - 
	0.015514037795505151*sh_m1_pow_m2*sh_m8*
	(6. + 51.*sh + 181.*sh_sq + 301.*sh_cube - 
	2324.*sh_four + 2220.*sh_five - 
	880.*sh_six - 70.*sh_seven + 95.*sh_eight))\
	- 0.015514037795505151*Lsh*sh_m1_pow_m7*sh_m8*
	(-12. - 42.*sh + 28.*sh_sq + 308.*sh_cube + 
	2793.*sh_four + 
	2.*(-8339. + 210.*Lshb)*sh_five + 
	(37042. - 1680.*Lshb)*sh_six + 
	168.*(-266. + 15.*Lshb)*sh_seven - 
	8.*(-3959. + 210.*Lshb)*sh_eight + 
	3.*(-4351. + 140.*Lshb)*sh_nine + 2720.*sh_ten
	) - 0.00025856729659175253*sh_m1_pow_m6*
	sh_m8*(162. - 36.*sh - 2914.*sh_sq - 
	5862.*sh_cube + 97530.*sh_four - 
	329152.*sh_five - 
	4200.*pi_sq*pow(-1. + sh,3)*
	sh_five + 497696.*sh_six - 
	286926.*sh_seven - 114303.*sh_eight + 
	259230.*sh_nine - 141000.*sh_ten + 
	22550.*sh_eleven - 2975.*sh_twelve + 
	600.*sh_thirteen))*pow(z,5.);
	else
	res += (7.5730667353201575 + 76.74295831935828*(1. - sh) + 
	333.4828305424486*msh_p1_pow_2 + 
	1010.6544814238534*msh_p1_pow_3 + 
	2487.030661606718*msh_p1_pow_4 + 
	log(z)*(9.618703433213193 + 
	58.798203244964526*(1. - sh) + 
	212.46474760944307*msh_p1_pow_2 + 
	578.673609772342*msh_p1_pow_3 + 
	1319.624054885668*msh_p1_pow_4))*pow(z,5.);}

	}
	return res;
}

#ifdef CACHE
double F_17im_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_17im(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_19re(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double Li3shb = creal(CLi3(1.-sh));
	double Li4sh = creal(CLi4(sh));
	double sqrtsh = sqrt(sh);
	double sqrt4sh = sqrt(4-sh);
	double ash = asin(sqrtsh/2);
	double cl2 = Cl2(2*ash);
	double cl3 = Cl3(2*ash);
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double sh_fourteen = sh*sh_thirteen;
	double sh_fifteen = sh*sh_fourteen;
	double sh_sixteen = sh*sh_fifteen;
	double sh_seventeen = sh*sh_sixteen;
	double sh_eighteen = sh*sh_seventeen;
	double sh_nineteen = sh*sh_eighteen;
	double sh_twenty = sh*sh_nineteen;
	double sh_twenty_one = sh*sh_twenty;
	double lz_sq = pow(lz,2);
	double lz_cube = lz*lz_sq;   
	double lz_four = lz*lz_cube;     
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;      
	double z_m3 = z_m2*z_m1;    
	double z_sqrt = pow(z,0.5);
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m4_sh_sqrt = pow(-((-4 + sh)*sh),0.5);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = -4.61812 - 0.0493827*Lsh - 0.06584362139917696*L*Lsh + 
	2.814815*lz + L*(-1.953360768175583 + 
	1.1851851851851851*lz) + 4.47441*sh + 0.74074*lz*sh + 
	14.4621*z - 0.592593*Lsh*z + 4.796605*lz*z + 
	71.3855*sh*z - 2.66667*Lsh*sh*z + 4.238385*lz*sh*z - 
	0.592595*Lsh*lz*sh*z - 1.0534979423868314*pow(L,2) - 
	0.2962975*z*lz_sq + 3.134725*sh*z*lz_sq - 
	0.098765375*z*lz_cube - 0.098765375*sh*z*lz_cube + 
	0.0493826875*sh*z*lz_four + 37.1282*sh_sq + 
	0.0164609*Lsh*sh_sq + 11.03105*lz*sh_sq + 
	212.74*z*sh_sq - 5.33333*Lsh*z*sh_sq - 
	10.96075*lz*z*sh_sq - 1.185185*Lsh*lz*z*sh_sq + 
	1.3333325*lz_sq*sh_sq + 
	14.2931*z*lz_sq*sh_sq + 
	0.098765375*lz_cube*sh_sq - 
	0.29629625*z*lz_cube*sh_sq + 
	0.148148125*z*lz_four*sh_sq + 79.7475*sh_cube + 
	0.0219479*Lsh*sh_cube + 28.65855*lz*sh_cube + 
	425.579*z*sh_cube - 8.2963*Lsh*z*sh_cube - 
	34.4008*lz*z*sh_cube - 1.77778*Lsh*lz*z*sh_cube + 
	2.814825*lz_sq*sh_cube + 
	32.33925*z*lz_sq*sh_cube + 
	0.29629625*lz_cube*sh_cube - 
	0.6913575*z*lz_cube*sh_cube + 
	0.29629625*z*lz_four*sh_cube + 
	L*sh_cube*(0.00020902736952119668 - 
	0.0037624926513815404*z_m3) - 
	0.0759415*sh_cube*z_m3 + 
	L*sh_sq*(0.0014109347442680777 - 
	0.025396825396825397*z_m2) - 
	0.403158*sh_sq*z_m2 - 
	0.00480894*sh_cube*z_m2 + 
	L*sh*(0.01316872427983539 - 
	0.23703703703703705*z_m1) - 
	2.48507*sh*z_m1 - 0.0613169*sh_sq*z_m1 - 
	1.81002*sh_cube*z_m1 - 
	0.4597295*lz*sh_cube*z_m1 - 
	0.04938275*lz_sq*sh_cube*z_m1 + 
	sh_cube*(7.79821*pow(z,-0.5) - 319.726*z_sqrt) - 
	23.3946*sh*z_sqrt - 140.368*sh_sq*z_sqrt + 
	3.8991*pow(z,1.5) - 16.0864*z_sq + 
	4.95977*Lsh*z_sq + 27.12195*lz*z_sq - 
	0.592595*Lsh*lz*z_sq - 18.1301*sh*z_sq + 
	18.6539*Lsh*sh*z_sq + 74.798*lz*sh*z_sq - 
	2.37037*Lsh*lz*sh*z_sq - 
	3.85185*lz_sq*z_sq - 
	0.5925925*Lsh*lz_sq*z_sq - 
	12.2963*sh*lz_sq*z_sq - 
	2.37037*Lsh*sh*lz_sq*z_sq - 
	0.4938275*lz_cube*z_sq - 
	1.3827125*sh*lz_cube*z_sq - 
	44.6829*sh_sq*z_sq + 
	40.786*Lsh*sh_sq*z_sq + 
	136.0075*lz*sh_sq*z_sq - 
	7.1111*Lsh*lz*sh_sq*z_sq - 
	29.77775*lz_sq*sh_sq*z_sq - 
	5.333325*Lsh*lz_sq*sh_sq*z_sq - 
	2.6666625*lz_cube*sh_sq*z_sq - 
	87.8946*sh_cube*z_sq + 
	70.2698*Lsh*sh_cube*z_sq + 
	208.806*lz*sh_cube*z_sq - 
	15.80245*Lsh*lz*sh_cube*z_sq - 
	56.79*lz_sq*sh_cube*z_sq - 
	9.481475*Lsh*lz_sq*sh_cube*z_sq - 
	4.345675*lz_cube*sh_cube*z_sq - 14.73*z_cube - 
	9.20287*Lsh*z_cube - 14.28805*lz*z_cube - 
	0.52675*Lsh*lz*z_cube - 72.89*sh*z_cube - 
	41.6104*Lsh*sh*z_cube - 34.0675*lz*sh*z_cube - 
	1.185185*Lsh*lz*sh*z_cube + 
	5.037025*lz_sq*z_cube + 
	0.7901225*Lsh*lz_sq*z_cube + 
	15.901225*sh*lz_sq*z_cube + 
	3.55555*Lsh*sh*lz_sq*z_cube - 
	137.203*sh_sq*z_cube - 
	111.356*Lsh*sh_sq*z_cube - 
	49.7185*lz*sh_sq*z_cube + 
	42.22225*lz_sq*sh_sq*z_cube + 
	9.481475*Lsh*lz_sq*sh_sq*z_cube - 
	279.268*sh_cube*z_cube - 
	231.893*Lsh*sh_cube*z_cube - 
	73.4265*lz*sh_cube*z_cube + 
	5.92595*Lsh*lz*sh_cube*z_cube + 
	82.81475*lz_sq*sh_cube*z_cube + 
	19.753075*Lsh*lz_sq*sh_cube*z_cube;
	else{
	
	if(0<=maxpow){
	if(sh<.900001)
	res += (11.39728026699467 - 1.0534979423868314*pow(L,2) - 
	0.03292181069958848*Li2sh*(-8. + 17.*sh)*
	sh_m1_pow_m1 - 0.03292181069958848*cl2*(2. + sh)*
	sqrt4sh*sqrtsh*sh_m2 - 
	1.*L*(-1.1193415637860082*Lsh + 
	0.13168724279835392*ash*(2. + sh)*sqrt4sh*sqrtsh*
	sh_m2 + 0.13168724279835392*(-2. + 29.*sh)*
	sh_m1) + 0.14814814814814814*pow(Lsh,2)*
	sh_m1_pow_m3*sh_sq - 
	1.*Lsh*(0.03292181069958848*ash*(2. + sh)*sqrt4sh*
	sqrtsh*sh_m2 + 
	0.01646090534979424*sh_m1_pow_m2*
	(-168. + 307.*sh - 121.*sh_sq + 
	2.*Lshb*(8. - 25.*sh + 17.*sh_sq))) - 
	0.03292181069958848*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m2*(-4. + 3.*sh + 18.*sh_sq - 
	16.*sh_cube + 5.*sh_four) - 
	0.09876543209876543*pow(ash,2)*sh_m1_pow_m4*
	sh_m1*(12. - 33.*sh + 18.*sh_sq - 
	4.*sh_four + sh_five) + 
	0.001828989483310471*sh_m1_pow_m4*sh_m1*
	(144. + 2.*sh*(-2068. + 
	75.*pi_sq) + 
	(15077. - 618.*pi_sq)*
	sh_sq + 102.*
	(-214. + 9.*pi_sq)*sh_cube + 
	(14249. - 615.*pi_sq)*
	sh_four + 2.*(-1753. + 
	78.*pi_sq)*sh_five));
	else
	res += 1.*(5.605483340182181 - 
	1.9301201034078073*(1. - sh) - 
	1.1708211541924651*msh_p1_pow_2 - 
	0.8303057591438543*msh_p1_pow_3 + 
	L*(-3.9138369114536733 - 
	1.136144688522592*(1. - sh) - 
	0.5573265536192934*msh_p1_pow_2 - 
	0.37357014711105585*msh_p1_pow_3 - 
	0.27973271946953165*msh_p1_pow_4) - 
	0.6410250317311895*msh_p1_pow_4 - 
	1.0534979423868314*pow(L,2));}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (-0.09876543209876543*sh*pow(Lsh,4)*sh_m1_pow_m3 + 
	0.5925925925925926*Li2sh*(1. + 2.*sh)*
	sh_m1_pow_m2 - 0.19753086419753085*(-3. + sh)*
	pow(Lsh,3)*sh_m1_pow_m2 - 
	7.111111111111111*L*sh_m1 - 
	1.*lz*(0.5925925925925926*Lsh*
	(2. + 2.*Lshb*(-1. + sh) + sh)*sh_m1_pow_m2 + 
	1.1851851851851851*Li2sh*sh_m1_pow_m1 - 
	0.5925925925925926*pow(Lsh,2)*sh_m1_pow_m1 - 
	0.19753086419753085*
	(144. + sh*(-135. + pi_sq))*
	sh_m1_pow_m1*sh_m1) + 
	0.7123300166871669*sh_m1_pow_m3*
	(-3. + 3.*sh + 2.*sh_sq) - 
	0.5925925925925926*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m1*(2. - 11.*sh + 3.*sh_sq) + 
	0.5925925925925926*Li3sh*sh_m1_pow_m3*
	(-3. - 3.*sh + 4.*sh_sq) - 
	1.*pow(Lsh,2)*(1.1851851851851851*Li2sh*sh*
	sh_m1_pow_m3 + 
	0.09876543209876543*sh_m1_pow_m3*
	(12. + 3.*Lshb*(-7. + 5.*sh) + 
	2.*sh*(-6. + pi_sq) + 
	3.*sh_sq)) - 
	3.5555555555555554*pow(ash,2)*sh_m1_pow_m4*
	(-1. + 8.*sh - 5.*sh_sq + sh_cube) - 
	1.*Lsh*(5.698640133497335*sh*sh_m1_pow_m3 - 
	4.7407407407407405*Li3sh*sh*sh_m1_pow_m3 + 
	0.5925925925925926*Li2sh*sh_m1_pow_m3*
	(-5. + sh + 2.*sh_sq) - 
	0.09876543209876543*sh_m1_pow_m3*sh_m1*
	(-288. + sh*(906. - 6.*Lshb - 
	7.*pi_sq) + 
	(-903. - 6.*Lshb + 5.*pi_sq)*
	sh_sq + 3.*(95. + 4.*Lshb)*sh_cube)) - 
	0.019753086419753086*sh_m1_pow_m4*sh_m1*
	(-4.*(-1. + sh)*pow(pi,4.)*sh_sq + 
	15.*(-1. + sh)*(-20. + 48.*sh + 
	3.*(-13. + 8.*Li4sh)*sh_sq + 11.*sh_cube)\
	+ 5.*pi_sq*
	(-3. + 8.*sh - 10.*sh_sq + 2.*sh_four)))*
	z;
	else
	res += (-0.6372090343808106 - 30.841194974865694*(1. - sh) - 
	45.8723194923496*msh_p1_pow_2 - 
	55.64666652217136*msh_p1_pow_3 + 
	log(z)*(-29.333333333333332 - 
	29.037037037037038*(1. - sh) - 
	28.921810699588477*msh_p1_pow_2 - 
	28.854320987654322*msh_p1_pow_3 - 
	28.8079012345679*msh_p1_pow_4) + 
	L*(-7.111111111111111 - 
	7.111111111111111*(1. - sh) - 
	7.111111111111111*msh_p1_pow_2 - 
	7.111111111111111*msh_p1_pow_3 - 
	7.111111111111111*msh_p1_pow_4) - 
	62.91148562558235*msh_p1_pow_4)*z;}

	if(3<=maxpow){
	if(sh<.900001)
	res += 0.3950617283950617*(-2. + sh)*pi_sq*
	sh_m1*pow(z,1.5);
	else
	res += (-3.8991029732698697 - 7.798205946539739*(1. - sh) - 
	7.798205946539739*msh_p1_pow_2 - 
	7.798205946539739*msh_p1_pow_3 - 
	7.798205946539739*msh_p1_pow_4)*pow(z,1.5);}

	if(4<=maxpow){
	if(sh<.900001)
	res += (-0.5925925925925926*Li3shb*(6. + sh)*sh_m1_pow_m3 - 
	0.5925925925925926*Li3sh*(-11. + 3.*sh)*
	sh_m1_pow_m3 - 0.09876543209876543*(2. + sh)*
	pow(Lsh,3)*sh_m1_pow_m3 - 
	1.*L*(3.5555555555555554*sh_m2 + 
	7.111111111111111*Lsh*sh_m2) + 
	0.5925925925925926*Li2sh*sh_m1_pow_m3*sh_m2*
	(-3. + 10.*sh - 9.*sh_sq + sh_cube) + 
	0.5925925925925926*cl2*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m1*(2. + 8.*sh_sq - 5.*sh_cube + 
	sh_four) - 0.5925925925925926*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(-2. - 4.*sh + 24.*sh_sq - 15.*sh_cube + 
	3.*sh_four) - 0.7123300166871669*sh_m1_pow_m4*
	sh_m2*(64. - 256.*sh + 375.*sh_sq - 
	236.*sh_cube + 59.*sh_four) + 
	0.14814814814814814*pow(Lsh,2)*sh_m1_pow_m4*
	sh_m2*(185. - 751.*sh + 
	3.*(379. + 2.*Lshb)*sh_sq - 
	1.*(777. + 8.*Lshb)*sh_cube + 
	2.*(100. + Lshb)*sh_four) - 
	1.*Lsh*(-0.5925925925925926*Li2sh*(-7. + 2.*sh)*
	sh_m1_pow_m3 - 
	3.5555555555555554*pow(ash,2)*sh_m1_pow_m4*
	(-1. - 3.*sh + sh_sq) - 
	0.04938271604938271*sh_m1_pow_m3*sh_m2*
	(-423. + 1236.*sh - 
	1.*(1179. + 46.*pi_sq)*
	sh_sq + 2.*
	(225. + pi_sq)*sh_cube + 
	12.*Lshb*(-3. + 10.*sh - 9.*sh_sq + 
	sh_cube)) - 
	0.5925925925925926*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(2. + 8.*sh_sq - 5.*sh_cube + sh_four)) - 
	1.7777777777777777*pow(ash,2)*sh_m1_pow_m4*
	(4. - 2.*sh - 24.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five) - 
	1.*lz_sq*(-0.2962962962962963*Lsh*(6. + sh)*
	sh_m1_pow_m3 - 
	0.14814814814814814*sh_m1_pow_m2*sh_m2*
	(188. - 378.*sh + 169.*sh_sq + 12.*sh_cube - 
	6.*sh_four + sh_five)) - 
	0.024691358024691357*sh_m1_pow_m4*sh_m2*
	(88. - 448.*sh - 3.*(-415. + 48.*cl3)*sh_sq - 
	2.*(917. + 216.*cl3)*sh_cube + 
	(1453. + 144.*cl3)*sh_four - 693.*sh_five + 
	210.*sh_six - 21.*sh_seven - 
	2.*pi_sq*
	(-603. + 2441.*sh - 3677.*sh_sq + 
	2457.*sh_cube - 662.*sh_four + 
	52.*sh_five - 16.*sh_six + 2.*sh_seven)) - 
	1.*lz*(0.5925925925925926*Li2sh*(6. + sh)*
	sh_m1_pow_m3 - 
	0.2962962962962963*pow(Lsh,2)*sh_m1_pow_m2 - 
	7.111111111111111*L*sh_m2 - 
	7.111111111111111*pow(ash,2)*sh_m1_pow_m4*
	(-1. - 3.*sh + sh_sq) + 
	0.2962962962962963*Lsh*sh_m1_pow_m3*sh_m2*
	(-187. + 566.*sh + 
	3.*(-191. + 4.*Lshb)*sh_sq + 
	2.*(96. + Lshb)*sh_cube) - 
	1.1851851851851851*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(2. + 8.*sh_sq - 5.*sh_cube + sh_four) + 
	0.04938271604938271*sh_m1_pow_m4*sh_m2*
	(432. - 1728.*sh + 
	(2511. + 8.*pi_sq)*
	sh_sq - 2.*
	(735. + 11.*pi_sq)*sh_cube\
	+ (99. + 2.*pi_sq)*
	sh_four + 225.*sh_five - 78.*sh_six + 
	9.*sh_seven)))*z_sq;
	else
	res += z_sq*(-336.0382986459167 - 
	695.8939576145283*(1. - sh) - 
	1040.3241576582361*msh_p1_pow_2 - 
	1365.7546667417848*msh_p1_pow_3 - 
	1671.7824637501126*msh_p1_pow_4 + 
	L*(-3.5555555555555554 + 
	7.111111111111111*msh_p1_pow_2 + 
	16.59259259259259*msh_p1_pow_3 + 
	27.85185185185185*msh_p1_pow_4) + 
	log(z)*(-21.887216272530303 + 
	9.991951790611958*(1. - sh) + 
	71.95973563934483*msh_p1_pow_2 + 
	151.94549827445954*msh_p1_pow_3 + 
	245.72125954180052*msh_p1_pow_4 + 
	L*(7.111111111111111 + 
	14.222222222222221*(1. - sh) + 
	21.333333333333332*msh_p1_pow_2 + 
	28.444444444444443*msh_p1_pow_3 + 
	35.55555555555556*msh_p1_pow_4)) + 
	(27.65432098765432 + 55.67901234567901*(1. - sh) + 
	83.6*msh_p1_pow_2 + 
	111.39753086419753*msh_p1_pow_3 + 
	139.20987654320987*msh_p1_pow_4)*pow(log(z),2));}

	if(5<=maxpow){
	if(sh<.900001)
	res += -0.05267489711934156*pi_sq*sh_m3*
	(3. + 14.*sh + 3.*sh_sq)*pow(z,2.5);
	else
	res += (-10.397607928719653 - 20.795215857439306*(1. - sh) - 
	32.75246497546691*msh_p1_pow_2 - 
	46.26935528280246*msh_p1_pow_3 - 
	61.34588677944595*msh_p1_pow_4)*pow(z,2.5);}

	if(6<=maxpow){
	if(sh<.900001)
	res += (-0.7901234567901234*Li3shb*(2. + sh)*sh_m1_pow_m4 - 
	1.*L*(-6.320987654320987*sh_m3 + 
	9.481481481481481*Lsh*sh_m3) + 
	1.1851851851851851*Li3sh*sh_m1_pow_m3*sh_m1 - 
	0.3950617283950617*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m1 + 0.4748866777914446*sh_m1_pow_m4*
	sh_m1*(3. + sh + 2.*sh_sq) - 
	0.06584362139917696*Li2sh*sh_m1_pow_m4*sh_m3*
	(-36. + 144.*sh - 234.*sh_sq + 208.*sh_cube - 
	73.*sh_four + 9.*sh_five) - 
	0.01646090534979424*pow(Lsh,2)*sh_m1_pow_m5*
	sh_m4*(-12. + 3682.*sh - 18347.*sh_sq + 
	(36889. + 36.*Lshb)*sh_cube - 
	2.*(18563. + 36.*Lshb)*sh_four + 
	4.*(4679. + 9.*Lshb)*sh_five - 3748.*sh_six - 
	18.*sh_seven) + 
	0.3950617283950617*pow(ash,2)*sh_m1_pow_m4*
	(16. - 22.*sh - 171.*sh_sq + 336.*sh_cube - 
	300.*sh_four + 141.*sh_five - 33.*sh_six + 
	3.*sh_seven) - lz_sq*
	(-0.3950617283950617*Lsh*(2. + sh)*sh_m1_pow_m4 + 
	0.01646090534979424*sh_m1_pow_m3*sh_m4*
	(-3. + 3685.*sh - 11046.*sh_sq + 
	11061.*sh_cube - 3604.*sh_four - 
	132.*sh_five + 321.*sh_six - 
	378.*sh_seven + 222.*sh_eight - 
	60.*sh_nine + 6.*sh_ten)) - 
	0.13168724279835392*ash*sh_m1_pow_m3*
	(-8. - 21.*sh - 17.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five)*sh_m4_sh_sqrt\
	- Lsh*(-2.3703703703703702*(2. + sh)*pow(ash,2)*
	sh_m1_pow_m4 + 
	1.1851851851851851*Li2sh*sh_m1_pow_m3*
	sh_m1 + 0.0054869684499314125*sh_m1_pow_m4*
	sh_m4*(30. + (12452. - 432.*Lshb)*sh + 
	2.*(-25315. + 864.*Lshb)*sh_sq + 
	36.*(-78.*Lshb + 
	5.*(424. + pi_sq))*
	sh_cube + (-50317. + 2496.*Lshb - 
	36.*pi_sq)*sh_four + 
	(12496. - 876.*Lshb + 
	72.*pi_sq)*sh_five + 
	27.*(3. + 4.*Lshb)*sh_six) + 
	0.3950617283950617*ash*sh_m1_pow_m3*
	(4. + 15.*sh - 29.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five)*
	sh_m4_sh_sqrt) - 
	1.*lz*(0.7901234567901234*Li2sh*(2. + sh)*
	sh_m1_pow_m4 - 
	4.7407407407407405*(2. + sh)*pow(ash,2)*
	sh_m1_pow_m4 - 
	9.481481481481481*L*sh_m3 - 
	0.5925925925925926*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m1 + 0.03292181069958848*Lsh*
	sh_m1_pow_m4*sh_m4*
	(-6. + 3682.*sh - 14713.*sh_sq + 
	22170.*sh_cube + 
	6.*(-2491. + 8.*Lshb)*sh_four + 
	6.*(631. + 4.*Lshb)*sh_five - 9.*sh_six) + 
	0.0054869684499314125*sh_m1_pow_m3*sh_m4*
	(15. + 12482.*sh - 37473.*sh_sq + 
	37278.*sh_cube - 12305.*sh_four + 
	498.*sh_five - 1215.*sh_six + 
	606.*sh_seven - 84.*sh_eight - 24.*sh_nine + 
	6.*sh_ten) + 
	0.7901234567901234*ash*sh_m1_pow_m3*
	(4. + 15.*sh - 29.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five)*
	sh_m4_sh_sqrt) - 
	0.0027434842249657062*sh_m1_pow_m4*sh_m4*
	(19. + 31332.*sh - 125971.*sh_sq + 
	189649.*sh_cube - 94.*sh_eleven + 
	6.*pi_sq*
	(-6. + 3838.*sh - 15355.*sh_sq + 
	23258.*sh_cube - 15790.*sh_four + 
	4134.*sh_five - 299.*sh_six + 
	440.*sh_seven - 400.*sh_eight + 
	188.*sh_nine - 44.*sh_ten + 4.*sh_eleven)\
	+ sh_four*(-125525. + 1728.*cl3 - 
	576.*cl2*sh_m4_sh_sqrt) + 
	12.*sh_five*(2563. + 72.*cl3 - 
	132.*cl2*sh_m4_sh_sqrt) + 
	2.*sh_ten*(451. + 
	72.*cl2*sh_m4_sh_sqrt) - 
	4.*sh_nine*(791. + 
	324.*cl2*sh_m4_sh_sqrt) + 
	sh_eight*(4810. + 
	4464.*cl2*sh_m4_sh_sqrt) + 
	sh_six*(-343. + 
	6336.*cl2*sh_m4_sh_sqrt) - 
	1.*sh_seven*(2371. + 
	7488.*cl2*sh_m4_sh_sqrt)))*z_cube;
	else
	res += pow(z,3.)*(-691.9182670329545 - 
	2017.613655726338*(1. - sh) - 
	3881.550454689942*msh_p1_pow_2 - 
	6209.7474017909*msh_p1_pow_3 - 
	8925.087722919276*msh_p1_pow_4 + 
	L*(6.320987654320987 + 
	28.444444444444443*(1. - sh) + 
	71.11111111111111*msh_p1_pow_2 + 
	137.4814814814815*msh_p1_pow_3 + 
	229.92592592592592*msh_p1_pow_4) + 
	log(z)*(65.87937511351353 + 
	322.6369826409675*(1. - sh) + 
	831.0074785910778*msh_p1_pow_2 + 
	1624.0261584651576*msh_p1_pow_3 + 
	2734.9197954441984*msh_p1_pow_4 + 
	L*(9.481481481481481 + 
	28.444444444444443*(1. - sh) + 
	56.888888888888886*msh_p1_pow_2 + 
	94.81481481481481*msh_p1_pow_3 + 
	142.22222222222223*msh_p1_pow_4)) + 
	(59.7037037037037 + 181.14567901234568*(1. - sh) + 
	362.4493827160494*msh_p1_pow_2 + 
	603.5673133450912*msh_p1_pow_3 + 
	905.1428571428571*msh_p1_pow_4)*pow(log(z),2));}

	if(7<=maxpow){
	if(sh<.900001)
	res += -0.004514991181657848*pi_sq*sh_m5*
	(15. + 108.*sh + 314.*sh_sq + 108.*sh_cube + 
	15.*sh_four)*pow(z,3.5);
	else
	res += (-24.954259028927165 - 74.8627770867815*(1. - sh) - 
	157.21183188224114*msh_p1_pow_2 - 
	279.48770112398427*msh_p1_pow_3 - 
	449.8450801732495*msh_p1_pow_4)*pow(z,3.5);}

	if(8<=maxpow){
	if(sh<.900001)
	res += (-1.*L*(-19.555555555555557*sh_m4 + 
	21.333333333333332*Lsh*sh_m4) - 
	2.374433388957223*sh_m1_pow_m3*sh_m2 + 
	1.9753086419753085*Li3sh*sh_m1_pow_m3*sh_m2 - 
	0.6584362139917695*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m2 - 0.09876543209876543*Li2sh*
	sh_m1_pow_m5*sh_m4*
	(54. - 268.*sh + 544.*sh_sq - 570.*sh_cube + 
	292.*sh_four - 58.*sh_five - 6.*sh_six + 
	3.*sh_seven) - 1.7777777777777777*sh*pow(ash,2)*
	sh_m1_pow_m4*(-24. - 3.*sh + 174.*sh_sq - 
	390.*sh_cube + 411.*sh_four - 241.*sh_five + 
	80.*sh_six - 14.*sh_seven + sh_eight) - 
	0.01646090534979424*pow(Lsh,2)*sh_m1_pow_m6*
	sh_m6*(6. + 14.*sh - 9747.*sh_sq + 
	58077.*sh_cube - 
	3.*(48391. + 20.*Lshb)*sh_four + 
	3.*(64669. + 60.*Lshb)*sh_five - 
	15.*(9737. + 12.*Lshb)*sh_six + 
	(58789. + 60.*Lshb)*sh_seven - 9906.*sh_eight + 
	27.*sh_nine - 9.*sh_ten) + 
	0.00823045267489712*lz_sq*sh_m6*
	(-3. - 25.*sh + 19344.*sh_sq + 73.*sh_four + 
	3.*sh_five + 180.*sh_six + 36.*sh_seven - 
	738.*sh_eight + 612.*sh_nine - 180.*sh_ten + 
	18.*sh_eleven) + 
	0.09876543209876543*ash*sh_m1_pow_m3*
	(-6. - 41.*sh - 239.*sh_sq + 671.*sh_cube - 
	680.*sh_four + 329.*sh_five - 77.*sh_six + 
	7.*sh_seven)*sh_m4_sh_sqrt - 
	1.*lz*(-21.333333333333332*L*sh_m4 - 
	0.9876543209876543*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 - 0.03292181069958848*Lsh*
	sh_m1_pow_m5*sh_m6*
	(-3. - 10.*sh + 9722.*sh_sq - 
	48430.*sh_cube + 96898.*sh_four - 
	97142.*sh_five + 48736.*sh_six - 
	9804.*sh_seven + 3.*sh_eight + 3.*sh_nine) - 
	0.00027434842249657066*sh_m1_pow_m4*sh_m6*
	(105. + 55.*sh + 963530.*sh_sq - 
	3.85233e6*sh_cube + 5.76189e6*sh_four - 
	3.806202e6*sh_five + 923844.*sh_six + 
	1154.*sh_seven + 54349.*sh_eight - 
	76131.*sh_nine + 4716.*sh_ten + 
	82440.*sh_eleven - 83970.*sh_twelve + 
	36900.*sh_thirteen - 7740.*sh_fourteen + 
	630.*sh_fifteen) - 
	1.1851851851851851*ash*
	(-6. - 13.*sh + 20.*sh_sq - 8.*sh_cube + 
	sh_four)*sh_m4_sh_sqrt) - 
	1.*Lsh*(1.9753086419753085*Li2sh*sh_m1_pow_m3*
	sh_m2 + 0.0027434842249657062*sh_m1_pow_m5*
	sh_m6*(-21. + 10.*sh + 
	2.*(-48097. + 972.*Lshb)*sh_sq + 
	(483970. - 9648.*Lshb)*sh_cube + 
	(-972667. + 19584.*Lshb - 
	600.*pi_sq)*sh_four + 
	8.*(122107. - 2565.*Lshb + 
	150.*pi_sq)*sh_five + 
	2.*(-245429. + 5256.*Lshb - 
	300.*pi_sq)*sh_six + 
	(98142. - 2088.*Lshb)*sh_seven - 
	6.*(7. + 36.*Lshb)*sh_eight + 
	6.*(5. + 18.*Lshb)*sh_nine) - 
	0.5925925925925926*ash*
	(-6. - 13.*sh + 20.*sh_sq - 8.*sh_cube + 
	sh_four)*sh_m4_sh_sqrt) - 
	4.572473708276177e-6*sh_m1_pow_m5*sh_m6*
	(-1200.*pi_sq*
	(-9. - 30.*sh + 30192.*sh_sq - 
	150534.*sh_cube + 301774.*sh_four - 
	303542.*sh_five + 152986.*sh_six - 
	30076.*sh_seven - 1033.*sh_eight - 
	2857.*sh_nine + 10122.*sh_ten - 
	14418.*sh_eleven + 11736.*sh_twelve - 
	5778.*sh_thirteen + 1692.*sh_fourteen - 
	270.*sh_fifteen + 18.*sh_sixteen) + 
	(-1. + sh)*(2775. - 2775.*sh + 
	1.836255e7*sh_sq - 7.415055e7*sh_cube + 
	1.1249425e8*sh_four - 7.6337698e7*sh_five + 
	86850.*sh_fifteen + 
	sh_eight*(4.883011e6 - 
	4.6656e6*cl2*sh_m4_sh_sqrt) - 
	900.*sh_fourteen*
	(1249. + 144.*cl2*sh_m4_sh_sqrt) + 
	900.*sh_thirteen*
	(6445. + 1728.*cl2*sh_m4_sh_sqrt) - 
	4050.*sh_twelve*
	(3755. + 1856.*cl2*sh_m4_sh_sqrt) + 
	4500.*sh_eleven*
	(4681. + 4176.*cl2*sh_m4_sh_sqrt) + 
	12.*sh_six*
	(1.746143e6 + 
	64800.*cl2*sh_m4_sh_sqrt) - 
	2.*sh_seven*(1.393817e6 + 
	712800.*cl2*sh_m4_sh_sqrt) - 
	36.*sh_ten*
	(374591. + 716400.*cl2*sh_m4_sh_sqrt)
	 + 3.*sh_nine*
	(-523. + 6.1344e6*cl2*sh_m4_sh_sqrt))
	))*pow(z,4.);
	else
	res += pow(z,4.)*(-1672.4687901566804 - 
	6425.515582048289*(1. - sh) - 
	15411.783043279445*msh_p1_pow_2 - 
	29568.2502704632*msh_p1_pow_3 - 
	49534.42028818623*msh_p1_pow_4 + 
	L*(19.555555555555557 + 
	99.55555555555556*(1. - sh) + 
	291.55555555555554*msh_p1_pow_2 + 
	654.2222222222222*msh_p1_pow_3 + 
	1251.5555555555557*msh_p1_pow_4) + 
	log(z)*(253.5498382001439 + 
	1369.6655601827902*(1. - sh) + 
	4072.304539664025*msh_p1_pow_2 + 
	9174.433775454903*msh_p1_pow_3 + 
	17618.760992657026*msh_p1_pow_4 + 
	L*(21.333333333333332 + 
	85.33333333333333*(1. - sh) + 
	213.33333333333334*msh_p1_pow_2 + 
	426.6666666666667*msh_p1_pow_3 + 
	746.6666666666666*msh_p1_pow_4)) + 
	(159.01234567901236 + 638.8148148148148*(1. - sh) + 
	1591.9506172839506*msh_p1_pow_2 + 
	3177.448559670782*msh_p1_pow_3 + 
	5557.119341563786*msh_p1_pow_4)*pow(log(z),2));}

	if(9<=maxpow){
	if(sh<.900001)
	res += -0.0010749979003947259*pi_sq*sh_m7*
	(35. + 330.*sh + 1389.*sh_sq + 3212.*sh_cube + 
	1389.*sh_four + 330.*sh_five + 35.*sh_six)*
	pow(z,4.5);
	else
	res += (-71.2978829397919 - 285.1915317591676*(1. - sh) - 
	745.0628767208254*msh_p1_pow_2 - 
	1586.3778954103698*msh_p1_pow_3 - 
	2982.415906901117*msh_p1_pow_4)*pow(z,4.5);}

	if(10<=maxpow){
	if(sh<.900001)
	res += (-1.*L*(-58.785185185185185*sh_m5 + 
	56.888888888888886*Lsh*sh_m5) - 
	4.986310116810168*sh_m1_pow_m3*sh_m3 + 
	4.148148148148148*Li3sh*sh_m1_pow_m3*sh_m3 - 
	1.382716049382716*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m3 - 0.3950617283950617*cl2*sh*sqrt4sh*sqrtsh*
	(44. + 108.*sh - 221.*sh_sq + 132.*sh_cube - 
	33.*sh_four + 3.*sh_five) - 
	0.03950617283950617*ash*sh*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*(-496. + 326.*sh + 6147.*sh_sq - 
	16495.*sh_cube + 18680.*sh_four - 
	11251.*sh_five + 3760.*sh_six - 
	658.*sh_seven + 47.*sh_eight) - 
	0.19753086419753085*Li2sh*sh_m1_pow_m6*sh_m5*
	(-72. + 429.*sh - 1073.*sh_sq + 1442.*sh_cube - 
	1092.*sh_four + 435.*sh_five - 47.*sh_six - 
	6.*sh_seven - sh_eight + sh_nine) + 
	1.1851851851851851*sh*pow(ash,2)*sh_m1_pow_m4*
	(-24. + 162.*sh + 114.*sh_sq - 1809.*sh_cube + 
	4224.*sh_four - 4957.*sh_five + 
	3423.*sh_six - 1445.*sh_seven + 
	366.*sh_eight - 51.*sh_nine + 3.*sh_ten) - 
	1.*Lsh*(4.148148148148148*Li2sh*sh_m1_pow_m3*
	sh_m3 + 0.3950617283950617*ash*sh*sqrt4sh*
	sqrtsh*(44. + 108.*sh - 221.*sh_sq + 
	132.*sh_cube - 33.*sh_four + 3.*sh_five) + 
	0.00005486968449931413*sh_m1_pow_m6*sh_m8*
	(486. - 108.*sh - 8742.*sh_sq + 
	(1.6072738e7 - 259200.*Lshb)*sh_cube + 
	6.*(-1.6103209e7 + 257400.*Lshb)*sh_four + 
	24.*(1.0092821e7 - 160950.*Lshb + 
	2625.*pi_sq)*sh_five + 
	4.*(-8.0991623e7 + 1.2978e6*Lshb - 
	47250.*pi_sq)*sh_six + 
	6.*(4.0636147e7 - 655200.*Lshb + 
	31500.*pi_sq)*sh_seven + 
	3.*(-3.2644121e7 + 522000.*Lshb - 
	21000.*pi_sq)*sh_eight - 
	2.*(-8.202737e6 + 84600.*Lshb)*sh_nine - 
	900.*(17. + 24.*Lshb)*sh_ten - 
	450.*(-1. + 8.*Lshb)*sh_eleven + 
	75.*(7. + 48.*Lshb)*sh_twelve)) - 
	0.0016460905349794238*pow(Lsh,2)*sh_m1_pow_m7*
	sh_m8*(-36. - 126.*sh + 84.*sh_sq + 
	299648.*sh_cube - 2.087099e6*sh_four + 
	60.*(104197. + 21.*Lshb)*sh_five - 
	2.*(5.211467e6 + 2520.*Lshb)*sh_six + 
	56.*(186416. + 135.*Lshb)*sh_seven - 
	12.*(523459. + 420.*Lshb)*sh_eight + 
	(2.104169e6 + 1260.*Lshb)*sh_nine - 
	303404.*sh_ten + 300.*sh_eleven + 
	120.*sh_twelve - 60.*sh_thirteen) + 
	0.0008230452674897119*lz_sq*sh_m8*
	(-18. - 189.*sh - 903.*sh_sq + 
	594928.*sh_cube + 630.*sh_four + 
	1383.*sh_five + 1839.*sh_six + 18.*sh_seven + 
	3120.*sh_eight + 8160.*sh_nine + 
	1440.*sh_ten - 46560.*sh_eleven + 
	51720.*sh_twelve - 23040.*sh_thirteen + 
	4680.*sh_fourteen - 360.*sh_fifteen) - 
	1.*lz*(-56.888888888888886*L*sh_m5 - 
	2.074074074074074*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 + 0.7901234567901234*ash*sh*sqrt4sh*
	sqrtsh*(44. + 108.*sh - 221.*sh_sq + 
	132.*sh_cube - 33.*sh_four + 3.*sh_five) - 
	0.0032921810699588477*Lsh*sh_m1_pow_m6*
	sh_m8*(18. + 81.*sh + 39.*sh_sq - 
	299147.*sh_cube + 1.788834e6*sh_four - 
	4.467993e6*sh_five + 5.959219e6*sh_six - 
	4.474299e6*sh_seven + 1.792497e6*sh_eight - 
	297944.*sh_nine - 450.*sh_ten + 
	90.*sh_eleven + 15.*sh_twelve) + 
	7.838526357044876e-6*sh_m1_pow_m5*sh_m8*
	(1701. + 1323.*sh - 29274.*sh_sq + 
	1.12541443e8*sh_cube - 
	5.61871835e8*sh_four + 
	1.121615551e9*sh_five - 
	1.117743382e9*sh_six + 
	5.54683001e8*sh_seven - 
	1.07408773e8*sh_eight - 4.03983e6*sh_nine + 
	4.371768e6*sh_ten + 4.551324e6*sh_eleven - 
	1.1897697e7*sh_twelve - 1.323264e7*sh_thirteen + 
	6.277908e7*sh_fourteen - 8.508528e7*sh_fifteen + 
	6.212682e7*sh_sixteen - 2.701314e7*sh_seventeen + 
	6.98418e6*sh_eighteen - 990360.*sh_nineteen + 
	59220.*sh_twenty)) - 
	2.7994736989445985e-8*sh_m1_pow_m6*sh_m8*
	(161406. - 318843.*sh - 1.910657e6*sh_sq + 
	2.6587841e7*sh_cube - 4.35567027e8*sh_four + 
	2.014488392e9*sh_five - 4.239785466e9*sh_six + 
	4.472946066e9*sh_seven - 1.480821324e9*sh_eight - 
	2.199797883e9*sh_nine + 
	2.944201467e9*sh_ten + 
	1.246759469e9*sh_eleven - 
	2.466843039e9*sh_twelve - 
	1.8171607842e10*sh_thirteen + 
	6.053910752e10*sh_fourteen - 
	8.996385496e10*sh_fifteen + 
	7.96292826e10*sh_sixteen - 
	4.509057616e10*sh_seventeen + 
	1.647058952e10*sh_eighteen - 
	3.74875284e9*sh_nineteen + 4.8245988e8*sh_twenty - 
	2.674812e7*sh_twenty_one + 
	58800.*pi_sq*
	(-18. - 81.*sh - 39.*sh_sq + 
	308147.*sh_cube - 1.843944e6*sh_four + 
	4.610871e6*sh_five - 6.160013e6*sh_six + 
	4.635393e6*sh_seven - 1.862435e6*sh_eight + 
	312614.*sh_nine + 3132.*sh_ten - 
	3736.*sh_eleven - 77479.*sh_twelve + 
	318588.*sh_thirteen - 608560.*sh_fourteen + 
	702440.*sh_fifteen - 529920.*sh_sixteen + 
	267160.*sh_seventeen - 89120.*sh_eighteen + 
	18840.*sh_nineteen - 2280.*sh_twenty + 
	120.*sh_twenty_one)))*pow(z,5.);
	else
	res += pow(z,5.)*(-4875.039604802492 - 
	23394.026843086856*(1. - sh) - 
	67666.42652379921*msh_p1_pow_2 - 
	152266.3652960926*msh_p1_pow_3 - 
	292874.8144929396*msh_p1_pow_4 + 
	L*(58.785185185185185 + 
	350.81481481481484*(1. - sh) + 
	1194.6666666666667*msh_p1_pow_2 + 
	3072.*msh_p1_pow_3 + 
	6641.777777777777*msh_p1_pow_4) + 
	log(z)*(847.677812246958 + 
	5397.215321438122*(1. - sh) + 
	18627.480751808722*msh_p1_pow_2 + 
	48121.18079162945*msh_p1_pow_3 + 
	104649.22541796633*msh_p1_pow_4 + 
	L*(56.888888888888886 + 
	284.44444444444446*(1. - sh) + 
	853.3333333333334*msh_p1_pow_2 + 
	1991.111111111111*msh_p1_pow_3 + 
	3982.222222222222*msh_p1_pow_4)) + 
	(491.23292181069957 + 
	2460.5102880658437*(1. - sh) + 
	7344.419753086419*msh_p1_pow_2 + 
	17100.213991769546*msh_p1_pow_3 + 
	34182.45267489712*msh_p1_pow_4)*pow(log(z),2));}

	}
	return res;
}

#ifdef CACHE
double F_19re_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_19re(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_19im(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double lz_sq = pow(lz,2);
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;      
	double z_m3 = z_m2*z_m1;            
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = 3.67166 + 0.20685383727340204*L - 0.103427*Lsh + 
	0.93084*lz - 0.310281*sh - 0.93084*lz*sh - 16.2155*z - 
	5.58505*lz*z - 30.7987*sh*z - 1.86168*Lsh*sh*z - 
	16.75515*lz*sh*z - 1.861685*z*lz_sq - 
	1.861685*sh*z*lz_sq - 1.36524*sh_sq - 
	0.93084*lz*sh_sq - 52.2081*z*sh_sq - 
	3.72337*Lsh*z*sh_sq - 26.0636*lz*z*sh_sq - 
	1.861685*z*lz_sq*sh_sq - 1.72206*sh_cube - 
	0.93084*lz*sh_cube - 76.6479*z*sh_cube - 
	5.58505*Lsh*z*sh_cube - 34.75145*lz*z*sh_cube - 
	1.861685*z*lz_sq*sh_cube - 
	0.00295505*sh_cube*z_m3 - 
	0.0199466*sh_sq*z_m2 + 
	0.00369382*sh_cube*z_m2 - 0.186168*sh*z_m1 + 
	0.0620562*sh_sq*z_m1 + 
	0.0871741*sh_cube*z_m1 + 26.7517*z_sq - 
	1.86168*Lsh*z_sq - 7.44675*lz*z_sq - 
	3.72337*Lsh*lz*z_sq + 66.1439*sh*z_sq - 
	7.44674*Lsh*sh*z_sq - 33.5103*lz*sh*z_sq - 
	14.8935*Lsh*lz*sh*z_sq - 
	7.44675*lz_sq*z_sq - 
	20.478525*sh*lz_sq*z_sq + 
	108.713*sh_sq*z_sq - 
	22.3402*Lsh*sh_sq*z_sq - 
	81.914*lz*sh_sq*z_sq - 
	33.5103*Lsh*lz*sh_sq*z_sq - 
	39.0955*lz_sq*sh_sq*z_sq + 
	148.481*sh_cube*z_sq - 
	49.6449*Lsh*sh_cube*z_sq - 
	155.761*lz*sh_cube*z_sq - 
	59.574*Lsh*lz*sh_cube*z_sq - 
	63.29725*lz_sq*sh_cube*z_sq - 
	23.6892*z_cube - 1.65483*Lsh*z_cube + 
	17.3757*lz*z_cube + 4.96449*Lsh*lz*z_cube - 
	63.7828*sh*z_cube - 3.72337*Lsh*sh*z_cube + 
	67.0205*lz*sh*z_cube + 22.3402*Lsh*lz*sh*z_cube - 
	106.832*sh_sq*z_cube + 
	165.0695*lz*sh_sq*z_cube + 
	59.574*Lsh*lz*sh_sq*z_cube - 
	135.118*sh_cube*z_cube + 
	18.6168*Lsh*sh_cube*z_cube + 
	326.4155*lz*sh_cube*z_cube + 
	124.1125*Lsh*lz*sh_cube*z_cube;
	else{
		
	if(0<=maxpow){
	if(sh<.900001)
	res += (-3.5165152336478345*L - 
	0.05171345931835051*sh_m1_pow_m2*
	(130. - 269.*sh + 121.*sh_sq) + 
	0.10342691863670102*Lsh*sh_m1_pow_m3*
	(-8. + 33.*sh - 51.*sh_sq + 17.*sh_cube));
	else
	res += (-5.636767065700205 - 1.3704066719362884*(1. - sh) - 
	3.5165152336478345*L - 
	0.5998761280928658*msh_p1_pow_2 - 
	0.36888934313756694*msh_p1_pow_3 - 
	0.26226111511449185*msh_p1_pow_4);}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (17.902805978847187*sh*sh_m1_pow_m3 - 
	14.893476283684944*Li3sh*sh*sh_m1_pow_m3 + 
	1.861684535460618*Li2sh*(-7. + 5.*sh)*
	sh_m1_pow_m3 + 1.241123023640412*sh*pow(Lsh,3)*
	sh_m1_pow_m3 - 3.723369070921236*pow(Lsh,2)*
	sh_m1_pow_m2 - Lsh*
	(-7.446738141842472*Li2sh*sh*sh_m1_pow_m3 + 
	1.861684535460618*sh_m1_pow_m3*
	(-5. + Lshb*(7. - 5.*sh) + 
	sh*(3. - 2.*pi_sq) + sh_sq
	)) - 0.310280755910103*sh_m1_pow_m3*
	sh_m1*(-306. + 
	sh*(912. - 11.*pi_sq) + 
	5.*(-180. + pi_sq)*sh_sq + 
	(294. + 4.*pi_sq)*sh_cube))*
	z;
	else
	res += (-101.4595749819256 - 99.28872579089392*(1. - sh) - 
	98.33919598335144*msh_p1_pow_2 - 
	97.78315475935152*msh_p1_pow_3 - 
	97.41019960703451*msh_p1_pow_4)*z;}

	/*if(3<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(4<=maxpow){
	if(sh<.900001)
	res += (-1.861684535460618*Li2sh*(3. + 2.*sh)*sh_m1_pow_m3 + 
	1.861684535460618*pow(Lsh,2)*sh_m1_pow_m2 + 
	22.340214425527417*L*sh_m2 - 
	1.*lz*(-13.031791748224327*Lsh*sh_m1_pow_m3 - 
	0.930842267730309*sh_m1_pow_m2*sh_m2*
	(193. - 393.*sh + 186.*sh_sq)) + 
	0.1551403779550515*sh_m1_pow_m3*sh_m2*
	(441. - 1368.*sh + 
	3.*(437. + 2.*pi_sq)*sh_sq + 
	4.*(-111. + pi_sq)*sh_cube + 
	24.*sh_four) - 
	0.930842267730309*Lsh*sh_m1_pow_m4*sh_m2*
	(191. - 777.*sh + (1175. - 6.*Lshb)*sh_sq + 
	(-797. + 2.*Lshb)*sh_cube + 
	(202. + 4.*Lshb)*sh_four))*z_sq;
	else
	res += (-74.72594871501649 + 29.277574993084137*(1. - sh) + 
	223.8210232757528*msh_p1_pow_2 + 
	478.4917845842381*msh_p1_pow_3 + 
	778.1470920996388*msh_p1_pow_4 + 
	L*(22.340214425527417 + 
	44.680428851054835*(1. - sh) + 
	67.02064327658225*msh_p1_pow_2 + 
	89.36085770210967*msh_p1_pow_3 + 
	111.70107212763709*msh_p1_pow_4) + 
	log(z)*(177.48059238057894 + 
	356.0471674068432*(1. - sh) + 
	535.0481354913816*msh_p1_pow_2 + 
	714.2663001050572*msh_p1_pow_3 + 
	893.6085770210967*msh_p1_pow_4))*z_sq;}

	/*if(5<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(6<=maxpow){
	if(sh<.900001)
	res += (29.786952567369887*L*sh_m3 + 
	3.723369070921236*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m1 - 1.241123023640412*Li2sh*sh_m1_pow_m4*
	sh_m1*(3. + sh + 2.*sh_sq) - 
	0.10342691863670102*Lsh*sh_m1_pow_m5*sh_m4*
	(12. - 3754.*sh + 18707.*sh_sq - 
	1.*(37645. + 36.*Lshb)*sh_cube + 
	6.*(6335. + 4.*Lshb)*sh_four - 
	6.*(3213. + 2.*Lshb)*sh_five + 
	24.*(163. + Lshb)*sh_six) - 
	1.*lz*(-1.241123023640412*Lsh*sh_m1_pow_m4*
	sh_m1*(3. + sh + 2.*sh_sq) + 
	0.10342691863670102*sh_m1_pow_m3*sh_m4*
	(-6. + 3748.*sh - 11253.*sh_sq + 
	11385.*sh_cube - 3841.*sh_four + 
	27.*sh_five + 12.*sh_six)) + 
	0.0172378197727835*sh_m1_pow_m4*sh_m4*
	(30. + 12740.*sh - 51314.*sh_sq + 
	36.*(2123. + pi_sq)*sh_cube + 
	(-50297. + 12.*pi_sq)*
	sh_four + 8.*(1474. + 
	3.*pi_sq)*sh_five + 
	381.*sh_six + 24.*sh_seven))*z_cube;
	else
	res += (211.9165849406685 + 1019.1512736698835*(1. - sh) + 
	2614.8727490925203*msh_p1_pow_2 + 
	5126.049113728145*msh_p1_pow_3 + 
	8646.890723525932*msh_p1_pow_4 + 
	L*(29.786952567369887 + 
	89.36085770210967*(1. - sh) + 
	178.72171540421934*msh_p1_pow_2 + 
	297.8695256736989*msh_p1_pow_3 + 
	446.80428851054836*msh_p1_pow_4) + 
	log(z)*(378.2322414544156 + 
	1142.5157994121814*(1. - sh) + 
	2290.864877035473*msh_p1_pow_2 + 
	3822.6086768805603*msh_p1_pow_3 + 
	5737.100041942259*msh_p1_pow_4))*z_cube;}

	/*if(7<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(8<=maxpow){
	if(sh<.900001)
	res += (67.02064327658225*L*sh_m4 + 
	6.2056151182020605*Li2sh*sh_m1_pow_m3*sh_m2 + 
	6.2056151182020605*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 - lz*
	(6.2056151182020605*Lsh*sh_m1_pow_m3*sh_m2 + 
	0.10342691863670102*sh_m1_pow_m2*sh_m6*
	(3. + 19.*sh - 9836.*sh_sq + 
	19672.*sh_cube - 9987.*sh_four + 
	39.*sh_five + 30.*sh_six)) + 
	0.10342691863670102*Lsh*sh_m1_pow_m6*sh_m6*
	(6. + 14.*sh - 9909.*sh_sq + 59043.*sh_cube - 
	3.*(49203. + 20.*Lshb)*sh_four + 
	3.*(65783. + 60.*Lshb)*sh_five - 
	3.*(49547. + 60.*Lshb)*sh_six + 
	(59839. + 60.*Lshb)*sh_seven - 10062.*sh_eight) - 
	0.00861890988639175*sh_m1_pow_m5*sh_m6*
	(21. - 10.*sh + 97976.*sh_sq - 
	490846.*sh_cube + 
	(981523. + 120.*pi_sq)*
	sh_four - 8.*(122341. + 
	30.*pi_sq)*sh_five + 
	2.*(242981. + 60.*pi_sq)*
	sh_six - 94596.*sh_seven - 678.*sh_eight - 
	408.*sh_nine + 144.*sh_ten))*pow(z,4.);
	else
	res += (818.5206340908517 + 4293.785518767847*(1. - sh) + 
	12756.23214766427*msh_p1_pow_2 + 
	28863.589912766278*msh_p1_pow_3 + 
	55512.63989848175*msh_p1_pow_4 + 
	L*(67.02064327658225 + 
	268.082573106329*(1. - sh) + 
	670.2064327658226*msh_p1_pow_2 + 
	1340.4128655316451*msh_p1_pow_3 + 
	2345.722514680379*msh_p1_pow_4) + 
	log(z)*(997.2423494950712 + 
	4001.6909089725987*(1. - sh) + 
	10016.483362289946*msh_p1_pow_2 + 
	20040.827170396282*msh_p1_pow_3 + 
	35069.73361639894*msh_p1_pow_4))*pow(z,4.);}

	/*if(9<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(10<=maxpow){
	if(sh<.900001)
	res += (178.72171540421934*L*sh_m5 + 
	13.031791748224327*Li2sh*sh_m1_pow_m3*sh_m3 + 
	13.031791748224327*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 - lz*
	(13.031791748224327*Lsh*sh_m1_pow_m3*sh_m3 + 
	0.0103426918636701*sh_m1_pow_m2*sh_m8*
	(18. + 153.*sh + 543.*sh_sq - 
	302141.*sh_cube + 603346.*sh_four - 
	304124.*sh_five + 60.*sh_six + 
	690.*sh_seven + 195.*sh_eight)) - 
	0.0103426918636701*Lsh*sh_m1_pow_m7*sh_m8*
	(36. + 126.*sh - 84.*sh_sq - 303968.*sh_cube + 
	2.117159e6*sh_four - 
	180.*(35233. + 7.*Lshb)*sh_five + 
	2.*(5.286917e6 + 2520.*Lshb)*sh_six - 
	56.*(189131. + 135.*Lshb)*sh_seven + 
	24.*(265547. + 210.*Lshb)*sh_eight - 
	7.*(304727. + 180.*Lshb)*sh_nine + 
	305864.*sh_ten) - 
	0.00017237819772783502*sh_m1_pow_m6*sh_m8*
	(-486. + 108.*sh + 8742.*sh_sq - 
	1.6340578e7*sh_cube + 9.7955934e7*sh_four - 
	2.44785504e8*sh_five + 
	12600.*pi_sq*pow(-1. + sh,3)*
	sh_five + 3.26063492e8*sh_six - 
	2.44009482e8*sh_seven + 9.7070643e7*sh_eight - 
	1.5960394e7*sh_nine + 92700.*sh_ten - 
	88650.*sh_eleven + 4275.*sh_twelve + 
	5400.*sh_thirteen))*pow(z,5.);
	else
	res += (2736.4139567610805 + 16812.602775356263*(1. - sh) + 
	58181.13183757553*msh_p1_pow_2 + 
	151135.2971056933*msh_p1_pow_3 + 
	329016.24546199996*msh_p1_pow_4 + 
	L*(178.72171540421934 + 
	893.6085770210967*(1. - sh) + 
	2680.8257310632903*msh_p1_pow_2 + 
	6255.260039147677*msh_p1_pow_3 + 
	12510.520078295354*msh_p1_pow_4) + 
	log(z)*(3064.787823810179 + 
	15358.380282956916*(1. - sh) + 
	46110.04743391064*msh_p1_pow_2 + 
	107597.09199613279*msh_p1_pow_3 + 
	215122.81941840626*msh_p1_pow_4))*pow(z,5.);}

	}
	return res;
}

#ifdef CACHE
double F_19im_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_19im(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double DeltaF_19re(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double lz_sq = pow(lz,2);
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;      
	double z_m3 = z_m2*z_m1;          
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	
	double res=0;

	if(sh<0.4)
	res = 0.09029982363315697*(0.6666666666666666 + L - 0.5*lz)*
	(105. + sh_cube*z_m3 + 4.5*sh_sq*z_m2 + 
	21.*sh*z_m1);
	else{
		
	/*if(0<=maxpow)
	res += 0.;*/

	/*if(1<=maxpow)
	res += 0.;*/

	if(2<=maxpow)
	res += (-37.925925925925924*sh_m1 - 
	56.888888888888886*L*sh_m1 + 
	28.444444444444443*lz*sh_m1)*z;

	if(3<=maxpow)
	res += 0.;

	if(4<=maxpow)
	res += (-75.85185185185185*Lsh*sh_m2 - 
	113.77777777777777*L*Lsh*sh_m2 - 
	56.888888888888886*lz_sq*sh_m2 - 
	1.*lz*(-75.85185185185185*sh_m2 - 
	113.77777777777777*L*sh_m2 - 
	56.888888888888886*Lsh*sh_m2))*z_sq;

	/*if(5<=maxpow)
	res += 0.;*/

	if(6<=maxpow)
	res += (151.7037037037037*sh_m3 - 
	151.7037037037037*Lsh*sh_m3 - 
	113.77777777777777*lz_sq*sh_m3 - 
	1.*lz*(-37.925925925925924*sh_m3 - 
	227.55555555555554*L*sh_m3 - 
	113.77777777777777*Lsh*sh_m3) - 
	1.*L*(-227.55555555555554*sh_m3 + 
	227.55555555555554*Lsh*sh_m3))*z_cube;

	/*if(7<=maxpow)
	res += 0.;*/

	if(8<=maxpow)
	res += (530.9629629629629*sh_m4 - 
	455.1111111111111*Lsh*sh_m4 - 
	341.3333333333333*lz_sq*sh_m4 - 
	1.*lz*(-56.888888888888886*sh_m4 - 
	682.6666666666666*L*sh_m4 - 
	341.3333333333333*Lsh*sh_m4) - 
	1.*L*(-796.4444444444445*sh_m4 + 
	682.6666666666666*Lsh*sh_m4))*pow(z,4.);

	/*if(9<=maxpow)
	res += 0.;*/

	if(10<=maxpow)
	res += (1871.0123456790122*sh_m5 - 
	1517.037037037037*Lsh*sh_m5 - 
	1137.7777777777778*lz_sq*sh_m5 - 
	1.*lz*(-113.77777777777777*sh_m5 - 
	2275.5555555555557*L*sh_m5 - 
	1137.7777777777778*Lsh*sh_m5) - 
	1.*L*(-2806.5185185185187*sh_m5 + 
	2275.5555555555557*Lsh*sh_m5))*pow(z,5.);

	}
	return res;
}

#ifdef CACHE
double DeltaF_19re_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=DeltaF_19re(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double DeltaF_19im(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;    
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;

	double res=0;

	if(sh<0.4)
	res = 0.;
	else{
	
	/*if(0<=maxpow)
	res += 0.;

	if(1<=maxpow)
	res += 0.;

	if(2<=maxpow)
	res += 0.;

	if(3<=maxpow)
	res += 0.;*/

	if(4<=maxpow)
	res += (238.2956205389591*sh_m2 + 
	357.4434308084387*L*sh_m2 - 
	178.72171540421934*lz*sh_m2)*z_sq;

	/*if(5<=maxpow)
	res += 0.;*/

	if(6<=maxpow)
	res += (476.5912410779182*sh_m3 + 
	714.8868616168774*L*sh_m3 - 
	357.4434308084387*lz*sh_m3)*z_cube;

	/*if(7<=maxpow)
	res += 0.;*/

	if(8<=maxpow)
	res += (1429.7737232337547*sh_m4 + 
	2144.660584850632*L*sh_m4 - 
	1072.330292425316*lz*sh_m4)*pow(z,4.);

	/*if(9<=maxpow)
	res += 0.;*/

	if(10<=maxpow)
	res += (4765.912410779182*sh_m5 + 
	7148.868616168774*L*sh_m5 - 
	3574.434308084387*lz*sh_m5)*pow(z,5.);

	}
	return res;
}

#ifdef CACHE
double DeltaF_19im_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=DeltaF_19im(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_27re(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double Li3shb = creal(CLi3(1.-sh));
	double Li4sh = creal(CLi4(sh));
	double sqrtsh = sqrt(sh);
	double sqrt4sh = sqrt(4-sh);
	double ash = asin(sqrtsh/2);
	double cl2 = Cl2(2*ash);
	double cl3 = Cl3(2*ash);   
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double sh_fourteen = sh*sh_thirteen;
	double sh_fifteen = sh*sh_fourteen;
	double sh_sixteen = sh*sh_fifteen;
	double sh_seventeen = sh*sh_sixteen;
	double sh_eighteen = sh*sh_seventeen;
	double sh_nineteen = sh*sh_eighteen;
	double sh_twenty = sh*sh_nineteen;
	double sh_twenty_one = sh*sh_twenty;
	double lz_sq = pow(lz,2);
	double lz_cube = lz*lz_sq;  
	double lz_four = lz*lz_cube;   
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;        
	double z_sqrt = pow(z,0.5);
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = 6.85597 + 5.135802469135802*L + 12.4502*sh + 
	0.0987654*Lsh*sh + 1.333335*lz*sh + 13.2214*z + 
	15.6523*lz*z + 155.555*sh*z - 6.22222*Lsh*sh*z + 
	10.20305*lz*sh*z - 1.77778*Lsh*lz*sh*z - 
	0.88889*z*lz_sq + 6.73755*sh*z*lz_sq - 
	0.29629625*z*lz_cube - 0.29629625*sh*z*lz_cube + 
	0.148148125*sh*z*lz_four + 116.815*sh_sq + 
	0.0987654*Lsh*sh_sq + 35.03385*lz*sh_sq + 
	542.972*z*sh_sq - 14.2222*Lsh*z*sh_sq - 
	44.79855*lz*z*sh_sq - 3.555555*Lsh*lz*z*sh_sq + 
	4.44445*lz_sq*sh_sq + 
	36.657*z*lz_sq*sh_sq + 
	0.29629625*lz_cube*sh_sq - 
	0.88888875*z*lz_cube*sh_sq + 
	0.444444375*z*lz_four*sh_sq + 251.971*sh_cube + 
	0.0987654*Lsh*sh_cube + 90.6275*lz*sh_cube + 
	1136.13*z*sh_cube - 23.1111*Lsh*z*sh_cube - 
	127.97*lz*z*sh_cube - 5.33335*Lsh*lz*z*sh_cube + 
	9.333325*lz_sq*sh_cube + 
	86.6475*z*lz_sq*sh_cube + 
	0.88888875*lz_cube*sh_cube - 
	2.074075*z*lz_cube*sh_cube + 
	0.8888875*z*lz_four*sh_cube + 
	0.000646678*sh_cube*z_m2 - 
	0.0333333*sh_sq*z_m1 - 
	5.68087*sh_cube*z_m1 - 
	1.466665*lz*sh_cube*z_m1 - 
	0.14814825*lz_sq*sh_cube*z_m1 + 
	sh_cube*(23.3946*pow(z,-0.5) - 959.179*z_sqrt) - 
	70.1839*sh*z_sqrt - 421.103*sh_sq*z_sqrt - 
	11.6973*pow(z,1.5) - 11.182*z_sq + 
	13.9904*lz*z_sq - 68.5374*sh*z_sq + 
	27.9808*Lsh*sh*z_sq + 102.242*lz*sh*z_sq - 
	15.55555*sh*lz_sq*z_sq - 
	3.55555*Lsh*sh*lz_sq*z_sq - 
	0.29629625*lz_cube*z_sq - 
	1.777775*sh*lz_cube*z_sq - 
	143.29*sh_sq*z_sq + 
	83.9424*Lsh*sh_sq*z_sq + 
	248.3745*lz*sh_sq*z_sq - 
	7.1111*Lsh*lz*sh_sq*z_sq - 
	48.4445*lz_sq*sh_sq*z_sq - 
	10.666675*Lsh*lz_sq*sh_sq*z_sq - 
	4.44445*lz_cube*sh_sq*z_sq - 
	271.07*sh_cube*z_sq + 
	164.329*Lsh*sh_cube*z_sq + 
	435.5445*lz*sh_cube*z_sq - 
	24.8889*Lsh*lz*sh_cube*z_sq - 
	106.37025*lz_sq*sh_cube*z_sq - 
	21.333325*Lsh*lz_sq*sh_cube*z_sq - 
	8.2963*lz_cube*sh_cube*z_sq + 7.26787*z_cube - 
	8.98765*lz*z_cube - 70.5057*sh*z_cube - 
	40.4253*Lsh*sh*z_cube - 56.869*lz*sh*z_cube - 
	3.555555*Lsh*lz*sh*z_cube + 
	6.222225*lz_sq*z_cube + 
	21.925925*sh*lz_sq*z_cube + 
	3.55555*Lsh*sh*lz_sq*z_cube - 
	228.849*sh_sq*z_cube - 
	165.257*Lsh*sh_sq*z_cube - 
	115.931*lz*sh_sq*z_cube - 
	7.1111*Lsh*lz*sh_sq*z_cube + 
	62.37025*lz_sq*sh_sq*z_cube + 
	14.222225*Lsh*lz_sq*sh_sq*z_cube - 
	464.161*sh_cube*z_cube - 
	416.697*Lsh*sh_cube*z_cube - 
	175.3475*lz*sh_cube*z_cube - 
	3.555555*Lsh*lz*sh_cube*z_cube + 
	144.14825*lz_sq*sh_cube*z_cube + 
	35.5555*Lsh*lz_sq*sh_cube*z_cube;
	else{
		
	if(0<=maxpow){
	if(sh<.900001)
	res += (5.135802469135802*L + 
	1.7777777777777777*sh*pow(ash,2)*sh_m1_pow_m4 - 
	0.04938271604938271*Lsh*sh*
	(-29. + 18.*Lshb*(-1. + sh) + 47.*sh)*sh_m1_pow_m2\
	- 0.8888888888888888*Li2sh*sh*sh_m1_pow_m1 + 
	0.4444444444444444*pow(Lsh,2)*sh_m1_pow_m3*
	sh_cube + 0.09876543209876543*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(-4. + 9.*sh - 15.*sh_sq + 4.*sh_cube) + 
	0.00823045267489712*sh_m1_pow_m4*
	(pow(-1. + sh,2)*(785. - 1600.*sh + 833.*sh_sq) + 
	6.*sh*pi_sq*
	(-4. + 9.*sh - 9.*sh_sq + 3.*sh_cube)));
	else
	res += (4.505947183821596 + 1.3035683102929574*(1. - sh) + 
	5.135802469135802*L + 
	0.39579451138872057*msh_p1_pow_2 + 
	0.18777843689890833*msh_p1_pow_3 + 
	0.10631777305882502*msh_p1_pow_4);}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (5.333333333333333*ash*(1. + sh)*sqrt4sh*sqrtsh*
	sh_m1_pow_m3 - 0.2962962962962963*sh*pow(Lsh,4)*
	sh_m1_pow_m3 + 1.7777777777777777*Li2sh*
	(-1. + 4.*sh)*sh_m1_pow_m2 - 
	0.5925925925925926*(-3. + sh)*pow(Lsh,3)*
	sh_m1_pow_m2 - lz*
	(16.*pow(Lsh,2)*pow(9 - 9*sh,-1) + 
	1.7777777777777777*Lsh*
	(2. + 2.*Lshb*(-1. + sh) + sh)*sh_m1_pow_m2 + 
	3.5555555555555554*Li2sh*sh_m1_pow_m1 - 
	0.5925925925925926*(9. + pi_sq)*
	sh_m1_pow_m1) - 
	1.7777777777777777*Li3sh*sh_m1_pow_m3*
	(1. + 5.*sh - 4.*sh_sq) + 
	2.1369900500615007*sh_m1_pow_m3*
	(-5. + 5.*sh + 2.*sh_sq) - 
	1.*pow(Lsh,2)*(3.5555555555555554*Li2sh*sh*
	sh_m1_pow_m3 + 
	0.2962962962962963*sh_m1_pow_m3*
	(12. + 3.*Lshb*(-5. + 3.*sh) + 
	sh*(-15. + 2.*pi_sq) + 
	6.*sh_sq)) - 
	1.*Lsh*(17.095920400492005*sh*sh_m1_pow_m3 - 
	14.222222222222221*Li3sh*sh*sh_m1_pow_m3 + 
	1.7777777777777777*Li2sh*sh_m1_pow_m3*
	(-3. - sh + 2.*sh_sq) - 
	0.2962962962962963*sh_m1_pow_m3*
	(12. + 21.*sh + 
	(-5. + 3.*sh)*pi_sq - 
	33.*sh_sq + 
	6.*Lshb*(1. - 5.*sh + 4.*sh_sq))) + 
	10.666666666666666*pow(ash,2)*sh_m1_pow_m4*
	(-1. - 3.*sh_sq + sh_cube) - 
	0.05925925925925926*sh_m1_pow_m4*
	(-4.*(-1. + sh)*sh*pow(pi,4.) - 
	15.*(-1. + sh)*(10. - (17. + 24.*Li4sh)*sh + 
	7.*sh_sq) + 
	5.*pi_sq*
	(-6. + 16.*sh - 20.*sh_sq + 7.*sh_cube)))*
	z;
	else
	res += (9.093880269665448 + 0.4118591346297644*(1. - sh) - 
	1.1764123627041028*msh_p1_pow_2 - 
	1.9623078784430308*msh_p1_pow_3 + 
	log(z)*(-2.6666666666666665 - 
	1.7777777777777777*(1. - sh) - 
	1.4320987654320987*msh_p1_pow_2 - 
	1.2296296296296296*msh_p1_pow_3 - 
	1.0903703703703704*msh_p1_pow_4) - 
	2.3203834816527733*msh_p1_pow_4)*z;}

	if(3<=maxpow){
	if(sh<.900001)
	res += -1.1851851851851851*(2. + sh)*pi_sq*
	sh_m1*pow(z,1.5);
	else
	res += (-35.09192675942883 - 23.39461783961922*(1. - sh) - 
	23.39461783961922*msh_p1_pow_2 - 
	23.39461783961922*msh_p1_pow_3 - 
	23.39461783961922*msh_p1_pow_4)*pow(z,1.5);}

	if(4<=maxpow){
	if(sh<.900001)
	res += (1.7777777777777777*Li3sh*(5. + 3.*sh)*sh_m1_pow_m3 - 
	1.7777777777777777*Li3shb*(1. + 6.*sh)*
	sh_m1_pow_m3 - 0.8888888888888888*pow(Lsh,3)*
	sh_m1_pow_m3 + 5.333333333333333*ash*(-3. + sh)*
	sqrt4sh*sqrtsh*sh_m1_pow_m3*sh_sq - 
	2.1369900500615007*sh_m1_pow_m4*
	(-5. + 8.*sh + 3.*sh_sq) - 
	1.7777777777777777*Li2sh*sh_m1_pow_m3*sh_m1*
	(3. - 10.*sh + 8.*sh_sq) + 
	5.333333333333333*pow(ash,2)*sh_m1_pow_m4*sh_sq*
	(2. + 9.*sh - 6.*sh_sq + sh_cube) - 
	1.7777777777777777*cl2*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	(-4. - 3.*sh_sq + sh_cube) - 
	1.*Lsh*(32.*sh*pow(ash,2)*sh_m1_pow_m4 + 
	1.7777777777777777*Li2sh*(4. + sh)*
	sh_m1_pow_m3 + 
	1.7777777777777777*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*(-4. - 3.*sh_sq + sh_cube) + 
	0.14814814814814814*sh_m1_pow_m3*sh_m2*
	(-18. + (57. + 36.*Lshb)*sh - 
	6.*(5. + 20.*Lshb + pi_sq)*
	sh_sq + (-93. + 96.*Lshb + 
	50.*pi_sq)*sh_cube)) + 
	0.4444444444444444*pow(Lsh,2)*sh_m1_pow_m4*
	sh_m2*(-4. + 7.*sh + (5. + 6.*Lshb)*sh_sq - 
	1.*(25. + 8.*Lshb)*sh_cube + 
	(11. + 2.*Lshb)*sh_four) - 
	1.*lz_sq*(-0.8888888888888888*Lsh*(1. + 6.*sh)*
	sh_m1_pow_m3 + 
	0.4444444444444444*sh_m1_pow_m2*sh_m2*
	(1. - 2.*sh + 8.*sh_sq + 10.*sh_cube - 
	4.*sh_four + sh_five)) - 
	1.*lz*(64.*sh*pow(ash,2)*sh_m1_pow_m4 + 
	1.7777777777777777*Li2sh*(1. + 6.*sh)*
	sh_m1_pow_m3 + 
	0.8888888888888888*pow(Lsh,2)*sh_m1_pow_m2 + 
	3.5555555555555554*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*(-4. - 3.*sh_sq + sh_cube) + 
	0.8888888888888888*Lsh*sh_m1_pow_m3*sh_m2*
	(2. - 3.*sh + 2.*(1. + Lshb)*sh_sq + 
	3.*(-1. + 4.*Lshb)*sh_cube) - 
	0.14814814814814814*sh_m1_pow_m4*sh_m2*
	(2.*pi_sq*sh_sq*
	(-1. + sh + 6.*sh_sq) + 
	3.*pow(-1. + sh,2)*
	(3. - 8.*sh + 8.*sh_sq + 20.*sh_cube - 
	14.*sh_four + 3.*sh_five))) - 
	0.07407407407407407*sh_m1_pow_m4*sh_m2*
	(-432.*cl3*sh_cube + 
	3.*pow(-1. + sh,2)*
	(7. - 64.*sh + 108.*sh_sq + 64.*sh_cube - 
	42.*sh_four + 7.*sh_five) + 
	2.*pi_sq*
	(-6. + 45.*sh - 133.*sh_sq + 157.*sh_cube - 
	59.*sh_four + 12.*sh_five - 12.*sh_six + 
	2.*sh_seven)))*z_sq;
	else
	res += z_sq*(3.8858229319704556 + 
	5.647852916167452*(1. - sh) + 
	5.03115350102676*msh_p1_pow_2 + 
	6.6118231676761265*msh_p1_pow_3 + 
	log(z)*(-10.515810650303592 - 
	3.3659736724084177*(1. - sh) - 
	7.052871337161167*msh_p1_pow_2 - 
	9.929763758132738*msh_p1_pow_3 - 
	12.684606769216838*msh_p1_pow_4) + 
	9.708947367679775*msh_p1_pow_4 + 
	(-0.5925925925925926 - 
	0.6666666666666666*(1. - sh) - 
	1.4222222222222223*msh_p1_pow_2 - 
	1.8074074074074074*msh_p1_pow_3 - 
	2.2222222222222223*msh_p1_pow_4)*pow(log(z),2));}

	if(5<=maxpow){
	if(sh<.900001)
	res += -0.1580246913580247*pi_sq*sh_m3*
	(3. + 14.*sh + 3.*sh_sq)*pow(z,2.5);
	else
	res += (-31.192823786158957 - 62.385647572317914*(1. - sh) - 
	98.25739492640072*msh_p1_pow_2 - 
	138.80806584840735*msh_p1_pow_3 - 
	184.03766033833784*msh_p1_pow_4)*pow(z,2.5);}

	if(6<=maxpow){
	if(sh<.900001)
	res += (-7.111111111111111*Li3shb*sh*sh_m1_pow_m4 + 
	3.5555555555555554*Li3sh*sh_m1_pow_m3*sh_m1 - 
	1.1851851851851851*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m1 + 4.273980100123001*sh_m1_pow_m4*
	sh_m1*(1. - sh + 2.*sh_sq) + 
	1.7777777777777777*Li2sh*sh_m1_pow_m4*sh_m2*
	(4. - 14.*sh + 18.*sh_sq - 11.*sh_cube + 
	sh_four) + 1.1851851851851851*cl2*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-9. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) + 0.3950617283950617*ash*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(27. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) - 3.5555555555555554*sh*pow(ash,2)*
	sh_m1_pow_m4*(-6. + 27.*sh - 30.*sh_cube + 
	27.*sh_four - 9.*sh_five + sh_six) - 
	1.*Lsh*(-21.333333333333332*sh*pow(ash,2)*
	sh_m1_pow_m4 + 
	3.5555555555555554*Li2sh*sh_m1_pow_m3*
	sh_m1 - 1.1851851851851851*ash*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-9. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) - 0.04938271604938271*sh_m1_pow_m4*
	sh_m4*(-10. + 10.*sh + 
	6.*(23. + 24.*Lshb)*sh_sq - 
	2.*(67. + 252.*Lshb + 
	30.*pi_sq)*sh_cube + 
	(-163. + 648.*Lshb + 
	60.*pi_sq)*sh_four - 
	12.*(-5. + 33.*Lshb + 
	6.*pi_sq)*sh_five + 
	9.*(-5. + 4.*Lshb)*sh_six)) - 
	0.14814814814814814*pow(Lsh,2)*sh_m1_pow_m5*
	sh_m4*(-4. - 4.*sh + 47.*sh_sq + 
	(-47. + 12.*Lshb)*sh_cube - 
	2.*(13. + 12.*Lshb)*sh_four + 
	4.*(19. + 3.*Lshb)*sh_five - 36.*sh_six + 
	6.*sh_seven) - lz_sq*
	(-3.5555555555555554*Lsh*sh*sh_m1_pow_m4 - 
	0.14814814814814814*sh_m1_pow_m3*sh_m4*
	(1. + 3.*sh - 18.*sh_sq + 25.*sh_cube - 
	36.*sh_four + 26.*sh_five - 27.*sh_six - 
	22.*sh_seven + 38.*sh_eight - 16.*sh_nine + 
	2.*sh_ten)) - 
	1.*lz*(7.111111111111111*Li2sh*sh*sh_m1_pow_m4 - 
	42.666666666666664*sh*pow(ash,2)*sh_m1_pow_m4 - 
	1.7777777777777777*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m1 - 2.3703703703703702*ash*sh*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-9. - sh + 9.*sh_sq - 6.*sh_cube + 
	sh_four) - 0.2962962962962963*Lsh*
	sh_m1_pow_m4*sh_m4*
	(2. + 4.*sh - 27.*sh_sq + 22.*sh_cube + 
	14.*sh_four - 12.*(1. + 2.*Lshb)*sh_five + 
	9.*sh_six) - 
	0.04938271604938271*sh_m1_pow_m3*sh_m4*
	(-5. + 57.*sh_sq - 46.*sh_cube - 
	85.*sh_four + 164.*sh_five - 73.*sh_six + 
	102.*sh_seven - 40.*sh_eight - 4.*sh_nine + 
	2.*sh_ten)) + 
	0.00823045267489712*sh_m1_pow_m4*sh_m4*
	(-19. + 46.*sh + 1365.*sh_sq - 4435.*sh_cube + 
	3757.*sh_four + (610. - 2592.*cl3)*sh_five - 
	3103.*sh_six + 2113.*sh_seven + 
	594.*sh_eight - 1548.*sh_nine + 
	714.*sh_ten - 94.*sh_eleven + 
	6.*pi_sq*
	(6. + 12.*sh - 195.*sh_sq + 442.*sh_cube - 
	434.*sh_four + 124.*sh_five + 
	137.*sh_six - 16.*sh_seven - 
	120.*sh_eight + 108.*sh_nine - 
	36.*sh_ten + 4.*sh_eleven)))*z_cube;
	else
	res += pow(z,3.)*(-12.335307775922313 + 
	30.140289422045853*(1. - sh) + 
	46.569840149902745*msh_p1_pow_2 + 
	80.11596831781813*msh_p1_pow_3 + 
	log(z)*(-16.636559572212125 + 
	1.4762567309815953*(1. - sh) - 
	23.052042837781816*msh_p1_pow_2 - 
	49.522893482074544*msh_p1_pow_3 - 
	78.97537511261919*msh_p1_pow_4) + 
	147.9535188895975*msh_p1_pow_4 + 
	(1.3333333333333333 - 
	0.26666666666666666*(1. - sh) - 
	5.807407407407408*msh_p1_pow_2 - 
	10.137566137566138*msh_p1_pow_3 - 
	16.084656084656086*msh_p1_pow_4)*pow(log(z),2));}

	if(7<=maxpow){
	if(sh<.900001)
	res += -0.013544973544973546*pi_sq*sh_m5*
	(15. + 108.*sh + 314.*sh_sq + 108.*sh_cube + 
	15.*sh_four)*pow(z,3.5);
	else
	res += (-74.8627770867815 - 224.58833126034452*(1. - sh) - 
	471.6354956467235*msh_p1_pow_2 - 
	838.4631033719528*msh_p1_pow_3 - 
	1349.5352405197486*msh_p1_pow_4)*pow(z,3.5);}

	if(8<=maxpow){
	if(sh<.900001)
	res += (-7.123300166871669*sh_m1_pow_m3*sh_m2 + 
	5.925925925925926*Li3sh*sh_m1_pow_m3*sh_m2 - 
	1.9753086419753085*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m2 - 1.7777777777777777*cl2*sh*sqrt4sh*sqrtsh*
	(8. + 6.*sh - 6.*sh_sq + sh_cube) + 
	0.2962962962962963*Li2sh*sh_m1_pow_m5*sh_m3*
	(-50. + 230.*sh - 420.*sh_sq + 380.*sh_cube - 
	146.*sh_four + 15.*sh_five) - 
	0.2962962962962963*ash*sh*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*(16. + 72.*sh + 18.*sh_sq - 
	203.*sh_cube + 189.*sh_four - 63.*sh_five + 
	7.*sh_six) + 5.333333333333333*pow(ash,2)*
	sh_m1_pow_m4*sh_sq*
	(36. - 60.*sh + 99.*sh_cube - 112.*sh_four + 
	54.*sh_five - 12.*sh_six + sh_seven) - 
	1.*Lsh*(5.925925925925926*Li2sh*sh_m1_pow_m3*
	sh_m2 + 1.7777777777777777*ash*sh*sqrt4sh*
	sqrtsh*(8. + 6.*sh - 6.*sh_sq + sh_cube) - 
	0.00823045267489712*sh_m1_pow_m5*sh_m6*
	(21. - 10.*sh - 307.*sh_sq - 
	1.*(283. + 1800.*Lshb)*sh_cube + 
	5.*(437. + 1656.*Lshb + 
	120.*pi_sq)*sh_four - 
	2.*(1081. + 7560.*Lshb + 
	600.*pi_sq)*sh_five + 
	(253. + 13680.*Lshb + 
	600.*pi_sq)*sh_six + 
	(861. - 5256.*Lshb)*sh_seven + 
	6.*(49. + 90.*Lshb)*sh_eight - 78.*sh_nine))\
	- 0.04938271604938271*pow(Lsh,2)*sh_m1_pow_m6*
	sh_m6*(6. + 14.*sh - 30.*sh_sq - 
	255.*sh_cube + 729.*sh_four + 
	60.*Lshb*pow(-1. + sh,3)*sh_four - 
	585.*sh_five - 165.*sh_six + 541.*sh_seven - 
	270.*sh_eight + 45.*sh_nine) - 
	0.024691358024691357*lz_sq*sh_m6*
	(3. + 25.*sh + 90.*sh_sq - 60.*sh_cube + 
	5.*sh_four - 3.*sh_five - 144.*sh_six - 
	252.*sh_seven + 72.*sh_eight + 288.*sh_nine - 
	144.*sh_ten + 18.*sh_eleven) - 
	1.*lz*(-2.962962962962963*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 + 3.5555555555555554*ash*sh*sqrt4sh*
	sqrtsh*(8. + 6.*sh - 6.*sh_sq + sh_cube) - 
	0.09876543209876543*Lsh*sh_m1_pow_m5*sh_m6*
	(-3. - 10.*sh + 5.*sh_sq + 185.*sh_cube - 
	425.*sh_four + 343.*sh_five - 
	125.*sh_six + 3.*sh_seven - 15.*sh_eight + 
	15.*sh_nine) - 
	0.0008230452674897119*sh_m1_pow_m4*sh_m6*
	(105. + 55.*sh - 1480.*sh_sq + 90.*sh_cube + 
	10470.*sh_four - 25782.*sh_five + 
	43242.*sh_six - 52678.*sh_seven + 
	52987.*sh_eight - 14283.*sh_nine - 
	23976.*sh_ten + 3240.*sh_eleven + 
	28080.*sh_twelve - 22680.*sh_thirteen + 
	6480.*sh_fourteen - 630.*sh_fifteen)) - 
	0.000013717421124828532*sh_m1_pow_m5*sh_m6*
	(-2775. + 5550.*sh + 34425.*sh_sq + 
	1.48905e6*sh_cube - 7.4002e6*sh_four + 
	1.4083748e7*sh_five - 1.4200216e7*sh_six + 
	9.34426e6*sh_seven - 5.568515e6*sh_eight + 
	141850.*sh_nine + 4.948647e6*sh_ten - 
	27324.*sh_eleven - 8.9613e6*sh_twelve + 
	9.8388e6*sh_thirteen - 4.6764e6*sh_fourteen + 
	1.03725e6*sh_fifteen - 86850.*sh_sixteen + 
	1200.*pi_sq*
	(9. + 30.*sh - 15.*sh_sq - 1245.*sh_cube + 
	4325.*sh_four - 6463.*sh_five + 
	5357.*sh_six - 2657.*sh_seven + 
	295.*sh_eight + 1405.*sh_nine - 
	1014.*sh_ten - 1782.*sh_eleven + 
	3798.*sh_twelve - 2988.*sh_thirteen + 
	1188.*sh_fourteen - 234.*sh_fifteen + 
	18.*sh_sixteen)))*pow(z,4.);
	else
	res += pow(z,4.)*(-24.99645604544183 + 
	176.89749924459056*(1. - sh) + 
	310.586379087356*msh_p1_pow_2 + 
	603.5051355323303*msh_p1_pow_3 + 
	log(z)*(-36.21914373721251 + 
	37.13272398412609*(1. - sh) - 
	85.48808015010837*msh_p1_pow_2 - 
	255.76022141568095*msh_p1_pow_3 - 
	462.83148077351*msh_p1_pow_4) + 
	1307.9886664553635*msh_p1_pow_4 + 
	(2.5185185185185186 - 
	1.4814814814814814*(1. - sh) - 
	30.666666666666668*msh_p1_pow_2 - 
	58.46913580246913*msh_p1_pow_3 - 
	107.30864197530865*msh_p1_pow_4)*pow(log(z),2));}

	if(9<=maxpow){
	if(sh<.900001)
	res += -0.0032249937011841773*pi_sq*sh_m7*
	(35. + 330.*sh + 1389.*sh_sq + 3212.*sh_cube + 
	1389.*sh_four + 330.*sh_five + 35.*sh_six)*
	pow(z,4.5);
	else
	res += (-213.8936488193757 - 855.5745952775028*(1. - sh) - 
	2235.1886301624763*msh_p1_pow_2 - 
	4759.133686231109*msh_p1_pow_3 - 
	8947.24772070335*msh_p1_pow_4)*pow(z,4.5);}

	if(10<=maxpow){
	if(sh<.900001)
	res += (-14.958930350430506*sh_m1_pow_m3*sh_m3 + 
	12.444444444444445*Li3sh*sh_m1_pow_m3*sh_m3 - 
	4.148148148148148*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m3 + 3.5555555555555554*cl2*sqrt4sh*sqrtsh*
	sh_sq*(-20. - 11.*sh + 24.*sh_sq - 
	9.*sh_cube + sh_four) + 
	0.11851851851851852*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_sq*(630. - 1853.*sh - 99.*sh_sq + 
	4458.*sh_cube - 5217.*sh_four + 
	2538.*sh_five - 564.*sh_six + 47.*sh_seven) - 
	0.5925925925925926*Li2sh*sh_m1_pow_m6*sh_m4*
	(-63. + 355.*sh - 826.*sh_sq + 1008.*sh_cube - 
	665.*sh_four + 213.*sh_five - 7.*sh_seven + 
	sh_eight) - 10.666666666666666*pow(ash,2)*
	sh_m1_pow_m4*sh_sq*
	(10. - 90.*sh + 165.*sh_sq - 333.*sh_four + 
	445.*sh_five - 275.*sh_six + 90.*sh_seven - 
	15.*sh_eight + sh_nine) - 
	1.*Lsh*(12.444444444444445*Li2sh*sh_m1_pow_m3*
	sh_m3 - 3.5555555555555554*ash*sqrt4sh*sqrtsh*
	sh_sq*(-20. - 11.*sh + 24.*sh_sq - 
	9.*sh_cube + sh_four) - 
	0.0004938271604938272*sh_m1_pow_m6*sh_m8*
	(-162. + 36.*sh + 2914.*sh_sq + 
	5862.*sh_cube + 
	30.*(-767. + 2520.*Lshb)*sh_four - 
	8.*(2131. + 53250.*Lshb + 
	2625.*pi_sq)*sh_five + 
	8.*(12113. + 123900.*Lshb + 
	7875.*pi_sq)*sh_six - 
	18.*(5693. + 67200.*Lshb + 
	3500.*pi_sq)*sh_seven + 
	3.*(15101. + 266000.*Lshb + 
	7000.*pi_sq)*sh_eight - 
	50.*(259. + 5112.*Lshb)*sh_nine + 
	9600.*sh_ten + 
	150.*(11. + 56.*Lshb)*sh_eleven - 
	25.*(25. + 48.*Lshb)*sh_twelve)) - 
	0.014814814814814815*pow(Lsh,2)*sh_m1_pow_m7*
	sh_m8*(-12. - 42.*sh + 28.*sh_sq + 
	308.*sh_cube + 1533.*sh_four - 
	8318.*sh_five + 
	420.*Lshb*pow(-1. + sh,4)*sh_five + 
	13422.*sh_six - 8008.*sh_seven - 
	1788.*sh_eight + 4507.*sh_nine - 
	1540.*sh_ten - 140.*sh_eleven + 
	160.*sh_twelve - 20.*sh_thirteen) - 
	0.007407407407407408*lz_sq*sh_m8*
	(6. + 63.*sh + 301.*sh_sq + 840.*sh_cube - 
	630.*sh_four + 99.*sh_five - 153.*sh_six - 
	6.*sh_seven - 1200.*sh_eight - 1440.*sh_nine - 
	3360.*sh_ten + 2400.*sh_eleven + 
	5400.*sh_twelve - 4800.*sh_thirteen + 
	1320.*sh_fourteen - 120.*sh_fifteen) - 
	1.*lz*(-6.222222222222222*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 - 7.111111111111111*ash*sqrt4sh*sqrtsh*
	sh_sq*(-20. - 11.*sh + 24.*sh_sq - 
	9.*sh_cube + sh_four) - 
	0.02962962962962963*Lsh*sh_m1_pow_m6*sh_m8*
	(6. + 27.*sh + 13.*sh_sq - 141.*sh_cube - 
	1380.*sh_four + 5549.*sh_five - 
	8207.*sh_six + 6207.*sh_seven - 
	2809.*sh_eight + 1080.*sh_nine - 
	30.*sh_ten - 30.*sh_eleven + 35.*sh_twelve)
	 + 0.00007054673721340387*sh_m1_pow_m5*
	sh_m8*(567. + 441.*sh - 9758.*sh_sq - 
	30275.*sh_cube + 166495.*sh_four - 
	20923.*sh_five - 866054.*sh_six + 
	1.788427e6*sh_seven - 1.255755e6*sh_eight - 
	458890.*sh_nine + 1.258076e6*sh_ten - 
	55832.*sh_eleven + 524481.*sh_twelve - 
	2.826e6*sh_thirteen + 901320.*sh_fourteen + 
	5.25168e6*sh_fifteen - 7.8771e6*sh_sixteen + 
	5.0799e6*sh_seventeen - 1.7073e6*sh_eighteen + 
	290640.*sh_nineteen - 19740.*sh_twenty)) - 
	8.398421096833795e-8*sh_m1_pow_m6*sh_m8*
	(161406. - 318843.*sh - 1.910657e6*sh_sq + 
	738969.*sh_cube - 6.33135615e8*sh_four + 
	3.654733992e9*sh_five - 8.580730746e9*sh_six + 
	1.0387902666e10*sh_seven - 
	6.272402572e9*sh_eight + 4.35513965e8*sh_nine + 
	2.310469467e9*sh_ten - 
	2.213537811e9*sh_eleven + 
	4.191847021e9*sh_twelve - 
	8.659496442e9*sh_thirteen + 
	5.00515056e9*sh_fourteen + 
	1.057709688e10*sh_fifteen - 
	2.294859924e10*sh_sixteen + 
	2.03364504e10*sh_seventeen - 
	1.00250472e10*sh_eighteen + 
	2.83732932e9*sh_nineteen - 4.2896364e8*sh_twenty + 
	2.674812e7*sh_twenty_one - 
	58800.*pi_sq*
	(18. + 81.*sh + 39.*sh_sq - 423.*sh_cube - 
	9810.*sh_four + 47729.*sh_five - 
	95467.*sh_six + 103887.*sh_seven - 
	63289.*sh_eight + 15730.*sh_nine + 
	8388.*sh_ten - 21404.*sh_eleven + 
	46029.*sh_twelve - 50988.*sh_thirteen - 
	20160.*sh_fourteen + 133320.*sh_fifteen - 
	179760.*sh_sixteen + 130200.*sh_seventeen - 
	56400.*sh_eighteen + 14520.*sh_nineteen - 
	2040.*sh_twenty + 120.*sh_twenty_one)))*pow(z,5.);
	else
	res += pow(z,5.)*(-56.97595269013286 + 
	864.4668402137501*(1. - sh) + 
	1722.9691178723424*msh_p1_pow_2 + 
	3775.432090316476*msh_p1_pow_3 + 
	log(z)*(-109.64319713305599 + 
	219.04552007798918*(1. - sh) - 
	343.43648455706074*msh_p1_pow_2 - 
	1304.051993190421*msh_p1_pow_3 - 
	2520.1905777627526*msh_p1_pow_4) + 
	9467.150417064422*msh_p1_pow_4 + 
	(9.481481481481481 - 1.7777777777777777*(1. - sh) - 
	150.*msh_p1_pow_2 - 
	309.9259259259259*msh_p1_pow_3 - 
	648.8888888888889*msh_p1_pow_4)*pow(log(z),2));}

	}
	return res;
}

#ifdef CACHE
double F_27re_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_27re(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_27im(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double lz_sq = pow(lz,2);
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;    
	double z_sq = pow(z,2.);
	double z_cube = pow(z,3.);
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = 3.10281 - 8.37758*sh - 2.792525*lz*sh - 9.55118*z - 
	5.58505*lz*z - 34.6839*sh*z - 5.58505*Lsh*sh*z - 
	39.0954*lz*sh*z - 5.58505*z*lz_sq - 
	5.58505*sh*z*lz_sq - 9.54113*sh_sq - 
	2.792525*lz*sh_sq - 88.6728*z*sh_sq - 
	11.1701*Lsh*z*sh_sq - 67.0205*lz*z*sh_sq - 
	5.58505*z*lz_sq*sh_sq - 9.82039*sh_cube - 
	2.792525*lz*sh_cube - 154.918*z*sh_cube - 
	16.7552*Lsh*z*sh_cube - 93.084*lz*z*sh_cube - 
	5.58505*z*lz_sq*sh_cube - 
	0.015514*sh_cube*z_m2 + 
	0.15514*sh_cube*z_m1 + 18.3741*z_sq + 
	91.4251*sh*z_sq - 33.5103*lz*sh*z_sq - 
	22.3402*Lsh*lz*sh*z_sq - 
	5.58505*lz_sq*z_sq - 
	27.92525*sh*lz_sq*z_sq + 
	196.813*sh_sq*z_sq - 
	22.3402*Lsh*sh_sq*z_sq - 
	117.286*lz*sh_sq*z_sq - 
	67.0205*Lsh*lz*sh_sq*z_sq - 
	67.02075*lz_sq*sh_sq*z_sq + 
	314.524*sh_cube*z_sq - 
	78.1908*Lsh*sh_cube*z_sq - 
	266.221*lz*sh_cube*z_sq - 
	134.0415*Lsh*lz*sh_cube*z_sq - 
	122.87125*lz_sq*sh_cube*z_sq - 
	17.3757*z_cube + 7.44675*lz*z_cube - 
	94.1903*sh*z_cube - 11.1701*Lsh*sh*z_cube + 
	74.4675*lz*sh*z_cube + 22.3402*Lsh*lz*sh*z_cube - 
	209.21*sh_sq*z_cube - 
	22.3402*Lsh*sh_sq*z_cube + 
	242.019*lz*sh_sq*z_cube + 
	89.361*Lsh*lz*sh_sq*z_cube - 
	325.499*sh_cube*z_cube - 
	11.1701*Lsh*sh_cube*z_cube + 
	554.78*lz*sh_cube*z_cube + 
	223.402*Lsh*lz*sh_cube*z_cube;
	else{
		
	if(0<=maxpow){
	if(sh<.900001)
	res += (-2.792526803190927*Lsh*sh*(-1. + 2.*sh)*sh_m1_pow_m3 + 
	0.1551403779550515*sh_m1_pow_m2*
	(20. - 49.*sh + 47.*sh_sq));
	else
	res += (4.964492094561648 - 0.6981317007977318*(1. - sh) - 
	0.32579479370560815*msh_p1_pow_2 - 
	0.18616845354606182*msh_p1_pow_3 - 
	0.11967972013675403*msh_p1_pow_4);}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (53.70841793654156*sh*sh_m1_pow_m3 - 
	44.680428851054835*Li3sh*sh*sh_m1_pow_m3 + 
	5.585053606381854*Li2sh*(-5. + 3.*sh)*
	sh_m1_pow_m3 + 3.723369070921236*sh*pow(Lsh,3)*
	sh_m1_pow_m3 - 11.170107212763709*pow(Lsh,2)*
	sh_m1_pow_m2 - Lsh*
	(-22.340214425527417*Li2sh*sh*sh_m1_pow_m3 - 
	5.585053606381854*sh_m1_pow_m3*
	(3. + Lshb*(-5. + 3.*sh) + 
	2.*sh*pi_sq - 2.*sh_sq))\
	- 0.930842267730309*sh_m1_pow_m3*
	(-6.*(2. - 5.*sh + 3.*sh_sq) + 
	pi_sq*
	(-9. + 3.*sh + 4.*sh_sq)))*z;
	else
	res += (-25.126044626684017 - 16.441531762218297*(1. - sh) - 
	12.739670260838084*msh_p1_pow_2 - 
	10.578200186941283*msh_p1_pow_3 - 
	9.130437128725532*msh_p1_pow_4)*z;}

	/*if(3<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(4<=maxpow){
	if(sh<.900001)
	res += (-5.585053606381854*Li2sh*(-2. + 7.*sh)*sh_m1_pow_m3 - 
	5.585053606381854*pow(Lsh,2)*sh_m1_pow_m2 - 
	1.*lz*(-39.09537524467298*Lsh*sh*sh_m1_pow_m3 + 
	2.792526803190927*sh_m1_pow_m2*sh_m2*
	(2. - 7.*sh + 15.*sh_sq + 4.*sh_cube)) + 
	0.4654211338651545*sh_m1_pow_m3*sh_m2*
	(-18. + 63.*sh - 2.*
	(63. + 2.*pi_sq)*sh_sq + 
	(9. + 14.*pi_sq)*sh_cube + 
	36.*sh_four) - 
	2.792526803190927*Lsh*sh_m1_pow_m4*sh_m2*
	(-4. + 13.*sh + (-21. + 4.*Lshb)*sh_sq + 
	(11. - 18.*Lshb)*sh_cube + 
	(-5. + 14.*Lshb)*sh_four))*z_sq;
	else
	res += (-9.153282299348039 - 18.43843391995787*(1. - sh) - 
	29.484428830357537*msh_p1_pow_2 - 
	43.167588536881716*msh_p1_pow_3 + 
	log(z)*(-3.723369070921236 - 
	6.050474740247009*(1. - sh) - 
	10.332349171806431*msh_p1_pow_2 - 
	15.265813190777068*msh_p1_pow_3 - 
	20.4785298900668*msh_p1_pow_4) - 
	59.199668549550246*msh_p1_pow_4)*z_sq;}

	/*if(5<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(6<=maxpow){
	if(sh<.900001)
	res += (11.170107212763709*pow(Lsh,2)*sh_m1_pow_m3*sh_m1 - 
	11.170107212763709*Li2sh*sh_m1_pow_m4*sh_m1*
	(1. - sh + 2.*sh_sq) - 
	0.930842267730309*Lsh*sh_m1_pow_m5*sh_m4*
	(4. + 4.*sh - 71.*sh_sq + 
	(155. - 12.*Lshb)*sh_cube + 
	2.*(-83. + 12.*Lshb)*sh_four + 
	(98. - 36.*Lshb)*sh_five + 
	12.*(-3. + 2.*Lshb)*sh_six) - 
	1.*lz*(-11.170107212763709*Lsh*sh_m1_pow_m4*
	sh_m1*(1. - sh + 2.*sh_sq) - 
	0.930842267730309*sh_m1_pow_m3*sh_m4*
	(2. + 6.*sh - 45.*sh_sq + 61.*sh_cube - 
	41.*sh_four - 11.*sh_five + 4.*sh_six)) - 
	0.1551403779550515*sh_m1_pow_m4*sh_m4*
	(-10. + 10.*sh + 54.*sh_sq + 
	(46. - 12.*pi_sq)*sh_cube + 
	3.*(-53. + 4.*pi_sq)*sh_four - 
	24.*(-10. + pi_sq)*sh_five - 
	141.*sh_six + 32.*sh_seven))*z_cube;
	else
	res += (-13.97814805375014 - 45.88121537642694*(1. - sh) - 
	113.65495437342528*msh_p1_pow_2 - 
	229.22586075329383*msh_p1_pow_3 + 
	log(z)*(-0.930842267730309 - 
	18.430676901060117*(1. - sh) - 
	51.75483008580519*msh_p1_pow_2 - 
	103.16391875788196*msh_p1_pow_3 - 
	174.70579590629714*msh_p1_pow_4) - 
	403.9971955539516*msh_p1_pow_4)*z_cube;}

	/*if(7<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(8<=maxpow){
	if(sh<.900001)
	res += (18.61684535460618*Li2sh*sh_m1_pow_m3*sh_m2 + 
	18.61684535460618*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 - lz*
	(18.61684535460618*Lsh*sh_m1_pow_m3*sh_m2 + 
	0.310280755910103*sh_m1_pow_m2*sh_m6*
	(3. + 19.*sh + 43.*sh_sq - 260.*sh_cube + 
	225.*sh_four - 105.*sh_five + 15.*sh_six))
	 + 0.310280755910103*Lsh*sh_m1_pow_m6*sh_m6*
	(6. + 14.*sh - 30.*sh_sq - 405.*sh_cube + 
	(1569. - 60.*Lshb)*sh_four + 
	15.*(-169. + 12.*Lshb)*sh_five + 
	(2235. - 180.*Lshb)*sh_six + 
	(-1037. + 60.*Lshb)*sh_seven + 213.*sh_eight) - 
	0.025856729659175254*sh_m1_pow_m5*sh_m6*
	(21. - 10.*sh - 307.*sh_sq + 1247.*sh_cube + 
	(-3143. + 120.*pi_sq)*
	sh_four + (3526. - 
	240.*pi_sq)*sh_five + 
	(397. + 120.*pi_sq)*sh_six - 
	2577.*sh_seven + 1140.*sh_eight + 
	138.*sh_nine - 72.*sh_ten))*pow(z,4.);
	else
	res += (-17.996283842785974 - 135.35554718355516*(1. - sh) - 
	472.35812505085903*msh_p1_pow_2 - 
	1200.9383413133833*msh_p1_pow_3 + 
	log(z)*(-14.893476283684944 - 
	84.70664636345812*(1. - sh) - 
	265.29004630313807*msh_p1_pow_2 - 
	624.2848808911273*msh_p1_pow_3 - 
	1241.9208884413238*msh_p1_pow_4) - 
	2545.7015752257666*msh_p1_pow_4)*pow(z,4.);}

	/*if(9<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(10<=maxpow){
	if(sh<.900001)
	res += (39.09537524467298*Li2sh*sh_m1_pow_m3*sh_m3 + 
	39.09537524467298*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 - lz*
	(39.09537524467298*Lsh*sh_m1_pow_m3*sh_m3 + 
	0.09308422677303091*sh_m1_pow_m2*sh_m8*
	(6. + 51.*sh + 181.*sh_sq + 301.*sh_cube - 
	2324.*sh_four + 2220.*sh_five - 
	880.*sh_six - 70.*sh_seven + 95.*sh_eight))\
	+ 0.09308422677303091*Lsh*sh_m1_pow_m7*sh_m8*
	(-12. - 42.*sh + 28.*sh_sq + 308.*sh_cube + 
	2793.*sh_four + 
	2.*(-8339. + 210.*Lshb)*sh_five + 
	(37042. - 1680.*Lshb)*sh_six + 
	168.*(-266. + 15.*Lshb)*sh_seven - 
	8.*(-3959. + 210.*Lshb)*sh_eight + 
	3.*(-4351. + 140.*Lshb)*sh_nine + 2720.*sh_ten
	) - 0.001551403779550515*sh_m1_pow_m6*sh_m8*
	(-162. + 36.*sh + 2914.*sh_sq + 5862.*sh_cube - 
	97530.*sh_four + 329152.*sh_five + 
	4200.*pi_sq*pow(-1. + sh,3)*
	sh_five - 497696.*sh_six + 
	286926.*sh_seven + 114303.*sh_eight - 
	259230.*sh_nine + 141000.*sh_ten - 
	22550.*sh_eleven + 2975.*sh_twelve - 
	600.*sh_thirteen))*pow(z,5.);
	else
	res += (-45.438400411920945 - 460.45774991614974*(1. - sh) - 
	2000.8969832546916*msh_p1_pow_2 - 
	6063.926888543121*msh_p1_pow_3 + 
	log(z)*(-57.71222059927916 - 
	352.7892194697871*(1. - sh) - 
	1274.7884856566582*msh_p1_pow_2 - 
	3472.041658634053*msh_p1_pow_3 - 
	7917.744329314009*msh_p1_pow_4) - 
	14922.183969640308*msh_p1_pow_4)*pow(z,5.);}

	}
	return res;
}

#ifdef CACHE
double F_27im_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_27im(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_29re(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double Li3shb = creal(CLi3(1.-sh));
	double Li4sh = creal(CLi4(sh));
	double sqrtsh = sqrt(sh);
	double sqrt4sh = sqrt(4-sh);
	double ash = asin(sqrtsh/2);
	double cl2 = Cl2(2*ash);
	double cl3 = Cl3(2*ash);
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double sh_fourteen = sh*sh_thirteen;
	double sh_fifteen = sh*sh_fourteen;
	double sh_sixteen = sh*sh_fifteen;
	double sh_seventeen = sh*sh_sixteen;
	double sh_eighteen = sh*sh_seventeen;
	double sh_nineteen = sh*sh_eighteen;
	double sh_twenty = sh*sh_nineteen;
	double sh_twenty_one = sh*sh_twenty;
	double lz_sq = pow(lz,2);
	double lz_cube = lz*lz_sq;       
	double lz_four = lz*lz_cube; 
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;      
	double z_m3 = z_m2*z_m1;   
	double z_sqrt = pow(z,0.5);
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = -24.2913 + 0.296296*Lsh + 0.3950617283950617*L*Lsh + 
	L*(1.0534979423868314 - 7.111111111111111*lz) - 
	11.55555*lz - 26.8464*sh - 4.444445*lz*sh - 86.7723*z + 
	3.55556*Lsh*z - 28.77965*lz*z - 428.313*sh*z + 
	16.*Lsh*sh*z - 25.4303*lz*sh*z + 3.555555*Lsh*lz*sh*z + 
	6.320987654320987*pow(L,2) + 1.7777775*z*lz_sq - 
	18.808425*sh*z*lz_sq + 0.5925925*z*lz_cube + 
	0.5925925*sh*z*lz_cube - 0.29629625*sh*z*lz_four - 
	222.769*sh_sq - 0.0987654*Lsh*sh_sq - 
	66.186*lz*sh_sq - 1276.44*z*sh_sq + 
	32.*Lsh*z*sh_sq + 65.7645*lz*z*sh_sq + 
	7.1111*Lsh*lz*z*sh_sq - 8.*lz_sq*sh_sq - 
	85.7585*z*lz_sq*sh_sq - 
	0.5925925*lz_cube*sh_sq + 
	1.777775*z*lz_cube*sh_sq - 
	0.8888875*z*lz_four*sh_sq - 478.485*sh_cube - 
	0.131687*Lsh*sh_cube - 171.951*lz*sh_cube - 
	2553.47*z*sh_cube + 49.7778*Lsh*z*sh_cube + 
	206.4045*lz*z*sh_cube + 10.66665*Lsh*lz*z*sh_cube - 
	16.8889*lz_sq*sh_cube - 
	194.03575*z*lz_sq*sh_cube - 
	1.777775*lz_cube*sh_cube + 
	4.14815*z*lz_cube*sh_cube - 
	1.777775*z*lz_four*sh_cube + 
	L*sh_cube*(-0.00125416421712718 + 
	0.022574955908289243*z_m3) - 
	0.0142243*sh_cube*z_m3 + 
	L*sh_sq*(-0.008465608465608466 + 
	0.1523809523809524*z_m2) - 
	0.0132191*sh_sq*z_m2 + 
	0.0288536*sh_cube*z_m2 + 0.8462*sh*z_m1 + 
	0.367901*sh_sq*z_m1 + 
	10.8601*sh_cube*z_m1 + 
	2.758375*lz*sh_cube*z_m1 + 
	0.2962975*lz_sq*sh_cube*z_m1 + 
	L*sh*(-0.07901234567901234 + 
	1.4222222222222223*z_m1) + 
	140.368*sh*z_sqrt + 842.206*sh_sq*z_sqrt + 
	sh_cube*(-46.7892*pow(z,-0.5) + 1918.36*z_sqrt) - 
	23.3946*pow(z,1.5) + 96.5187*z_sq - 
	29.7586*Lsh*z_sq - 162.7315*lz*z_sq + 
	3.555555*Lsh*lz*z_sq + 108.781*sh*z_sq - 
	111.923*Lsh*sh*z_sq - 448.7875*lz*sh*z_sq + 
	14.2222*Lsh*lz*sh*z_sq + 
	23.1111*lz_sq*z_sq + 
	3.55555*Lsh*lz_sq*z_sq + 
	73.77775*sh*lz_sq*z_sq + 
	14.222225*Lsh*sh*lz_sq*z_sq + 
	2.9629625*lz_cube*z_sq + 
	8.2963*sh*lz_cube*z_sq + 
	268.098*sh_sq*z_sq - 
	244.716*Lsh*sh_sq*z_sq - 
	816.045*lz*sh_sq*z_sq + 
	42.66665*Lsh*lz*sh_sq*z_sq + 
	178.66675*lz_sq*sh_sq*z_sq + 
	32.*Lsh*lz_sq*sh_sq*z_sq + 
	16.*lz_cube*sh_sq*z_sq + 
	527.368*sh_cube*z_sq - 
	421.619*Lsh*sh_cube*z_sq - 
	1252.835*lz*sh_cube*z_sq + 
	94.815*Lsh*lz*sh_cube*z_sq + 
	340.74*lz_sq*sh_cube*z_sq + 
	56.889*Lsh*lz_sq*sh_cube*z_sq + 
	26.074125*lz_cube*sh_cube*z_sq + 
	88.3801*z_cube + 55.2172*Lsh*z_cube + 
	85.7285*lz*z_cube + 3.160495*Lsh*lz*z_cube + 
	437.34*sh*z_cube + 249.663*Lsh*sh*z_cube + 
	204.405*lz*sh*z_cube + 7.1111*Lsh*lz*sh*z_cube - 
	30.22225*lz_sq*z_cube - 
	4.74075*Lsh*lz_sq*z_cube - 
	95.4075*sh*lz_sq*z_cube - 
	21.333325*Lsh*sh*lz_sq*z_cube + 
	823.218*sh_sq*z_cube + 
	668.137*Lsh*sh_sq*z_cube + 
	298.311*lz*sh_sq*z_cube - 
	253.3325*lz_sq*sh_sq*z_cube - 
	56.889*Lsh*lz_sq*sh_sq*z_cube + 
	1675.61*sh_cube*z_cube + 
	1391.36*Lsh*sh_cube*z_cube + 
	440.5585*lz*sh_cube*z_cube - 
	35.55555*Lsh*lz*sh_cube*z_cube - 
	496.89*lz_sq*sh_cube*z_cube - 
	118.5185*Lsh*lz_sq*sh_cube*z_cube;
	else{
		
	if(0<=maxpow){
	if(sh<.900001)
	res += (8.547960200246003 + 6.320987654320987*pow(L,2) + 
	0.19753086419753085*Li2sh*(-8. + 17.*sh)*
	sh_m1_pow_m1 + 0.19753086419753085*cl2*(2. + sh)*
	sqrt4sh*sqrtsh*sh_m2 - 
	1.*L*(6.716049382716049*Lsh - 
	0.7901234567901234*ash*(2. + sh)*sqrt4sh*sqrtsh*
	sh_m2 + 0.3950617283950617*(4. - 31.*sh)*
	sh_m1) - 0.8888888888888888*pow(Lsh,2)*
	sh_m1_pow_m3*sh_sq - 
	1.*Lsh*(-0.19753086419753085*ash*(2. + sh)*sqrt4sh*
	sqrtsh*sh_m2 - 
	0.09876543209876543*sh_m1_pow_m2*
	(-114. + 199.*sh - 67.*sh_sq + 
	2.*Lshb*(8. - 25.*sh + 17.*sh_sq))) + 
	0.19753086419753085*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m2*(-4. + 3.*sh + 18.*sh_sq - 
	16.*sh_cube + 5.*sh_four) + 
	0.5925925925925926*pow(ash,2)*sh_m1_pow_m4*
	sh_m1*(12. - 33.*sh + 18.*sh_sq - 
	4.*sh_four + sh_five) - 
	0.0054869684499314125*sh_m1_pow_m4*sh_m1*
	(288. + 5.*sh*(79. + 60.*pi_sq) - 
	2.*(2257. + 618.*pi_sq)*
	sh_sq + 6.*(1391. + 
	306.*pi_sq)*sh_cube - 
	10.*(617. + 123.*pi_sq)*
	sh_four + (1655. + 
	312.*pi_sq)*sh_five));
	else
	res += 1.*(-4.256813794434612 + 
	6.247387287113509*(1. - sh) + 
	4.358260258488125*msh_p1_pow_2 + 
	3.2040567770853476*msh_p1_pow_3 + 
	2.512816857053804*msh_p1_pow_4 + 
	L*(12.81635480205537 + 
	6.816868131135554*(1. - sh) + 
	3.34395932171576*msh_p1_pow_2 + 
	2.241420882666335*msh_p1_pow_3 + 
	1.6783963168171898*msh_p1_pow_4) + 
	6.320987654320987*pow(L,2));}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (0.5925925925925926*sh*pow(Lsh,4)*sh_m1_pow_m3 - 
	3.5555555555555554*Li2sh*(1. + 2.*sh)*
	sh_m1_pow_m2 + 1.1851851851851851*(-3. + sh)*
	pow(Lsh,3)*sh_m1_pow_m2 + 
	42.666666666666664*L*sh_m1 - 
	1.*lz*(64.*Li2sh*pow(9 - 9*sh,-1) - 
	3.5555555555555554*Lsh*
	(2. + 2.*Lshb*(-1. + sh) + sh)*sh_m1_pow_m2 + 
	3.5555555555555554*pow(Lsh,2)*sh_m1_pow_m1 + 
	1.1851851851851851*
	(-18. + sh*(27. + pi_sq))*
	sh_m1_pow_m1*sh_m1) - 
	4.273980100123001*sh_m1_pow_m3*
	(-3. + 3.*sh + 2.*sh_sq) + 
	3.5555555555555554*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m1*(2. - 11.*sh + 3.*sh_sq) - 
	3.5555555555555554*Li3sh*sh_m1_pow_m3*
	(-3. - 3.*sh + 4.*sh_sq) - 
	1.*pow(Lsh,2)*(-7.111111111111111*Li2sh*sh*
	sh_m1_pow_m3 - 
	0.5925925925925926*sh_m1_pow_m3*
	(12. + 3.*Lshb*(-7. + 5.*sh) + 
	2.*sh*(-6. + pi_sq) + 
	3.*sh_sq)) + 
	21.333333333333332*pow(ash,2)*sh_m1_pow_m4*
	(-1. + 8.*sh - 5.*sh_sq + sh_cube) - 
	1.*Lsh*(-34.19184080098401*sh*sh_m1_pow_m3 + 
	28.444444444444443*Li3sh*sh*sh_m1_pow_m3 - 
	3.5555555555555554*Li2sh*sh_m1_pow_m3*
	(-5. + sh + 2.*sh_sq) + 
	0.5925925925925926*sh_m1_pow_m3*sh_m1*
	(36. - sh*(66. + 6.*Lshb + 
	7.*pi_sq) + 
	(69. - 6.*Lshb + 5.*pi_sq)*
	sh_sq + 3.*(-13. + 4.*Lshb)*sh_cube)) + 
	0.11851851851851852*sh_m1_pow_m4*sh_m1*
	(-4.*(-1. + sh)*pow(pi,4.)*sh_sq + 
	15.*(-1. + sh)*(-56. + 156.*sh + 
	3.*(-49. + 8.*Li4sh)*sh_sq + 47.*sh_cube)\
	+ 5.*pi_sq*
	(-3. + 8.*sh - 10.*sh_sq + 2.*sh_four)))*
	z;
	else
	res += (67.82325420628487 + 57.04716984919418*(1. - sh) + 
	51.233916954097545*msh_p1_pow_2 + 
	45.87999913302814*msh_p1_pow_3 + 
	log(z)*(-16. - 17.77777777777778*(1. - sh) - 
	18.469135802469136*msh_p1_pow_2 - 
	18.874074074074073*msh_p1_pow_3 - 
	19.152592592592594*msh_p1_pow_4) + 
	41.468913753494085*msh_p1_pow_4 + 
	L*(42.666666666666664 + 
	42.666666666666664*(1. - sh) + 
	42.666666666666664*msh_p1_pow_2 + 
	42.666666666666664*msh_p1_pow_3 + 
	42.666666666666664*msh_p1_pow_4))*z;}

	if(3<=maxpow){
	if(sh<.900001)
	res += -2.3703703703703702*(-2. + sh)*pi_sq*
	sh_m1*pow(z,1.5);
	else
	res += (23.39461783961922 + 46.78923567923844*(1. - sh) + 
	46.78923567923844*msh_p1_pow_2 + 
	46.78923567923844*msh_p1_pow_3 + 
	46.78923567923844*msh_p1_pow_4)*pow(z,1.5);}

	if(4<=maxpow){
	if(sh<.900001)
	res += (3.5555555555555554*Li3shb*(6. + sh)*sh_m1_pow_m3 + 
	3.5555555555555554*Li3sh*(-11. + 3.*sh)*
	sh_m1_pow_m3 + 0.5925925925925926*(2. + sh)*
	pow(Lsh,3)*sh_m1_pow_m3 - 
	1.*L*(-21.333333333333332*sh_m2 - 
	42.666666666666664*Lsh*sh_m2) - 
	3.5555555555555554*Li2sh*sh_m1_pow_m3*sh_m2*
	(-3. + 10.*sh - 9.*sh_sq + sh_cube) - 
	3.5555555555555554*cl2*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	sh_m1*(2. + 8.*sh_sq - 5.*sh_cube + 
	sh_four) + 3.5555555555555554*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(-2. - 4.*sh + 24.*sh_sq - 15.*sh_cube + 
	3.*sh_four) - 4.273980100123001*sh_m1_pow_m4*
	sh_m2*(8. - 32.*sh + 57.*sh_sq - 
	52.*sh_cube + 13.*sh_four) - 
	0.8888888888888888*pow(Lsh,2)*sh_m1_pow_m4*
	sh_m2*(-31. + 113.*sh + 
	3.*(-53. + 2.*Lshb)*sh_sq + 
	(87. - 8.*Lshb)*sh_cube + 
	2.*(-8. + Lshb)*sh_four) - 
	1.*Lsh*(3.5555555555555554*Li2sh*(-7. + 2.*sh)*
	sh_m1_pow_m3 + 
	21.333333333333332*pow(ash,2)*sh_m1_pow_m4*
	(-1. - 3.*sh + sh_sq) + 
	0.2962962962962963*sh_m1_pow_m3*sh_m2*
	(333. - 1032.*sh + 
	(1089. - 46.*pi_sq)*
	sh_sq + 2.*
	(-153. + pi_sq)*sh_cube + 
	12.*Lshb*(-3. + 10.*sh - 9.*sh_sq + 
	sh_cube)) + 
	3.5555555555555554*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(2. + 8.*sh_sq - 5.*sh_cube + sh_four)) + 
	10.666666666666666*pow(ash,2)*sh_m1_pow_m4*
	(4. - 2.*sh - 24.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five) - 
	1.*lz_sq*(1.7777777777777777*Lsh*(6. + sh)*
	sh_m1_pow_m3 + 
	0.8888888888888888*sh_m1_pow_m2*sh_m2*
	(-28. + 54.*sh - 47.*sh_sq + 12.*sh_cube - 
	6.*sh_four + sh_five)) - 
	1.*lz*(-3.5555555555555554*Li2sh*(6. + sh)*
	sh_m1_pow_m3 + 
	1.7777777777777777*pow(Lsh,2)*sh_m1_pow_m2 + 
	42.666666666666664*L*sh_m2 + 
	42.666666666666664*pow(ash,2)*sh_m1_pow_m4*
	(-1. - 3.*sh + sh_sq) - 
	1.7777777777777777*Lsh*sh_m1_pow_m3*sh_m2*
	(29. - 82.*sh + 3.*(25. + 4.*Lshb)*sh_sq + 
	2.*(-12. + Lshb)*sh_cube) + 
	7.111111111111111*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*sh_m1*
	(2. + 8.*sh_sq - 5.*sh_cube + sh_four) - 
	0.2962962962962963*sh_m1_pow_m4*sh_m2*
	(2.*pi_sq*sh_sq*
	(4. - 11.*sh + sh_sq) + 
	3.*pow(-1. + sh,2)*
	(-108. + 216.*sh - 135.*sh_sq + 
	32.*sh_cube - 20.*sh_four + 3.*sh_five))
	) - 0.14814814814814814*sh_m1_pow_m4*
	(-144.*cl3*(-1. - 3.*sh + sh_sq) + 
	sh_m2*(pow(-1. + sh,2)*
	(-232. + 560.*sh - 757.*sh_sq + 
	336.*sh_cube - 168.*sh_four + 
	21.*sh_five) + 
	2.*pi_sq*
	(45. - 151.*sh + 211.*sh_sq - 
	135.*sh_cube - 14.*sh_four + 
	52.*sh_five - 16.*sh_six + 2.*sh_seven))
	))*z_sq;
	else
	res += z_sq*(-165.12748700917908 - 
	411.3508120821892*(1. - sh) - 
	670.1268907046218*msh_p1_pow_2 - 
	925.5677817546755*msh_p1_pow_3 + 
	L*(21.333333333333332 - 
	42.666666666666664*msh_p1_pow_2 - 
	99.55555555555556*msh_p1_pow_3 - 
	167.11111111111111*msh_p1_pow_4) + 
	log(z)*(-92.67670236481818 - 
	123.95171074367174*(1. - sh) - 
	143.75841383606902*msh_p1_pow_2 - 
	143.67298964675717*msh_p1_pow_3 + 
	L*(-42.666666666666664 - 
	85.33333333333333*(1. - sh) - 
	128.*msh_p1_pow_2 - 
	170.66666666666666*msh_p1_pow_3 - 
	213.33333333333334*msh_p1_pow_4) - 
	130.3275572508033*msh_p1_pow_4) - 
	1177.4249452560534*msh_p1_pow_4 + 
	(26.074074074074073 + 
	49.925925925925924*(1. - sh) + 
	74.4*msh_p1_pow_2 + 
	99.61481481481482*msh_p1_pow_3 + 
	124.74074074074075*msh_p1_pow_4)*pow(log(z),2));}

	if(5<=maxpow){
	if(sh<.900001)
	res += 0.3160493827160494*pi_sq*sh_m3*
	(3. + 14.*sh + 3.*sh_sq)*pow(z,2.5);
	else
	res += (62.385647572317914 + 124.77129514463583*(1. - sh) + 
	196.51478985280144*msh_p1_pow_2 + 
	277.6161316968147*msh_p1_pow_3 + 
	368.0753206766757*msh_p1_pow_4)*pow(z,2.5);}

	if(6<=maxpow){
	if(sh<.900001)
	res += (4.7407407407407405*Li3shb*(2. + sh)*sh_m1_pow_m4 - 
	1.*L*(37.925925925925924*sh_m3 - 
	56.888888888888886*Lsh*sh_m3) - 
	7.111111111111111*Li3sh*sh_m1_pow_m3*sh_m1 + 
	2.3703703703703702*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m1 - 2.8493200667486676*sh_m1_pow_m4*
	sh_m1*(3. + sh + 2.*sh_sq) + 
	2.3703703703703702*cl2*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	(4. + 15.*sh - 29.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five) + 
	0.7901234567901234*ash*sqrt4sh*sqrtsh*sh_m1_pow_m3*
	(-8. - 21.*sh - 17.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five) + 
	0.3950617283950617*Li2sh*sh_m1_pow_m4*sh_m3*
	(-36. + 144.*sh - 234.*sh_sq + 208.*sh_cube - 
	73.*sh_four + 9.*sh_five) - 
	1.*Lsh*(14.222222222222221*(2. + sh)*pow(ash,2)*
	sh_m1_pow_m4 - 
	7.111111111111111*Li2sh*sh_m1_pow_m3*
	sh_m1 - 2.3703703703703702*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*
	(4. + 15.*sh - 29.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five) - 
	0.03292181069958848*sh_m1_pow_m4*sh_m4*
	(30. + (1508. - 432.*Lshb)*sh + 
	2.*(-3427. + 864.*Lshb)*sh_sq + 
	36.*(296. - 78.*Lshb + 
	5.*pi_sq)*sh_cube + 
	(-6541. + 2496.*Lshb - 
	36.*pi_sq)*sh_four + 
	(1552. - 876.*Lshb + 
	72.*pi_sq)*sh_five + 
	27.*(3. + 4.*Lshb)*sh_six)) - 
	2.3703703703703702*pow(ash,2)*sh_m1_pow_m4*
	(16. - 22.*sh - 171.*sh_sq + 336.*sh_cube - 
	300.*sh_four + 141.*sh_five - 33.*sh_six + 
	3.*sh_seven) - 0.09876543209876543*pow(Lsh,2)*
	sh_m1_pow_m5*sh_m4*
	(12. + 494.*sh - 2533.*sh_sq + 
	(4871. - 36.*Lshb)*sh_cube + 
	(-4634. + 72.*Lshb)*sh_four - 
	4.*(-541. + 9.*Lshb)*sh_five - 428.*sh_six + 
	18.*sh_seven) - 
	1.*lz_sq*(2.3703703703703702*Lsh*(2. + sh)*
	sh_m1_pow_m4 - 
	0.09876543209876543*sh_m1_pow_m3*sh_m4*
	(-3. - 491.*sh + 1482.*sh_sq - 
	1467.*sh_cube + 572.*sh_four - 
	132.*sh_five + 321.*sh_six - 
	378.*sh_seven + 222.*sh_eight - 
	60.*sh_nine + 6.*sh_ten)) - 
	1.*lz*(-4.7407407407407405*Li2sh*(2. + sh)*
	sh_m1_pow_m4 + 
	28.444444444444443*(2. + sh)*pow(ash,2)*
	sh_m1_pow_m4 + 
	56.888888888888886*L*sh_m3 + 
	3.5555555555555554*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m1 - 4.7407407407407405*ash*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*
	(4. + 15.*sh - 29.*sh_sq + 23.*sh_cube - 
	8.*sh_four + sh_five) + 
	0.19753086419753085*Lsh*sh_m1_pow_m4*sh_m4*
	(6. + 494.*sh - 1991.*sh_sq + 
	2886.*sh_cube - 
	6.*(293. + 8.*Lshb)*sh_four - 
	6.*(-65. + 4.*Lshb)*sh_five + 9.*sh_six) - 
	0.03292181069958848*sh_m1_pow_m3*sh_m4*
	(15. + 1538.*sh - 4641.*sh_sq + 
	4446.*sh_cube - 1361.*sh_four + 
	498.*sh_five - 1215.*sh_six + 
	606.*sh_seven - 84.*sh_eight - 24.*sh_nine + 
	6.*sh_ten)) + 
	0.01646090534979424*sh_m1_pow_m4*sh_m4*
	(19. - 7836.*sh + 30701.*sh_sq - 
	45359.*sh_cube + (31147. + 1728.*cl3)*sh_four + 
	12.*(-701. + 72.*cl3)*sh_five - 343.*sh_six - 
	2371.*sh_seven + 4810.*sh_eight - 
	3164.*sh_nine + 902.*sh_ten - 
	94.*sh_eleven + 
	6.*pi_sq*
	(-6. - 338.*sh + 1349.*sh_sq - 
	1798.*sh_cube + 914.*sh_four - 
	42.*sh_five - 299.*sh_six + 
	440.*sh_seven - 400.*sh_eight + 
	188.*sh_nine - 44.*sh_ten + 4.*sh_eleven))
	)*z_cube;
	else
	res += pow(z,3.)*(-563.8946426367565 - 
	1680.2345038491285*(1. - sh) - 
	3329.6412593857735*msh_p1_pow_2 - 
	5423.45927216734*msh_p1_pow_3 + 
	L*(-37.925925925925924 - 
	170.66666666666666*(1. - sh) - 
	426.6666666666667*msh_p1_pow_2 - 
	824.8888888888889*msh_p1_pow_3 - 
	1379.5555555555557*msh_p1_pow_4) - 
	7956.241038705334*msh_p1_pow_4 + 
	log(z)*(-34.97995438478488 - 
	30.044118068027284*(1. - sh) + 
	62.84401734242212*msh_p1_pow_2 + 
	320.43564180164736*msh_p1_pow_3 + 
	L*(-56.888888888888886 - 
	170.66666666666666*(1. - sh) - 
	341.3333333333333*msh_p1_pow_2 - 
	568.8888888888889*msh_p1_pow_3 - 
	853.3333333333334*msh_p1_pow_4) + 
	749.5923384459219*msh_p1_pow_4) + 
	(54.22222222222222 + 150.45925925925926*(1. - sh) + 
	299.97037037037035*msh_p1_pow_2 + 
	503.0405643738977*msh_p1_pow_3 + 
	755.8095238095239*msh_p1_pow_4)*pow(log(z),2));}

	if(7<=maxpow){
	if(sh<.900001)
	res += 0.02708994708994709*pi_sq*sh_m5*
	(15. + 108.*sh + 314.*sh_sq + 108.*sh_cube + 
	15.*sh_four)*pow(z,3.5);
	else
	res += (149.725554173563 + 449.17666252068904*(1. - sh) + 
	943.270991293447*msh_p1_pow_2 + 
	1676.9262067439056*msh_p1_pow_3 + 
	2699.0704810394973*msh_p1_pow_4)*pow(z,3.5);}

	if(8<=maxpow){
	if(sh<.900001)
	res += (-1.*L*(117.33333333333333*sh_m4 - 
	128.*Lsh*sh_m4) + 
	14.246600333743338*sh_m1_pow_m3*sh_m2 - 
	11.851851851851851*Li3sh*sh_m1_pow_m3*sh_m2 + 
	3.950617283950617*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m2 - 3.5555555555555554*cl2*sqrt4sh*sqrtsh*
	(-6. - 13.*sh + 20.*sh_sq - 8.*sh_cube + 
	sh_four) + 0.5925925925925926*Li2sh*
	sh_m1_pow_m5*sh_m4*
	(54. - 268.*sh + 544.*sh_sq - 570.*sh_cube + 
	292.*sh_four - 58.*sh_five - 6.*sh_six + 
	3.*sh_seven) - 0.5925925925925926*ash*sqrt4sh*
	sqrtsh*sh_m1_pow_m3*
	(-6. - 41.*sh - 239.*sh_sq + 671.*sh_cube - 
	680.*sh_four + 329.*sh_five - 77.*sh_six + 
	7.*sh_seven) + 10.666666666666666*sh*pow(ash,2)*
	sh_m1_pow_m4*(-24. - 3.*sh + 174.*sh_sq - 
	390.*sh_cube + 411.*sh_four - 241.*sh_five + 
	80.*sh_six - 14.*sh_seven + sh_eight) - 
	1.*Lsh*(-11.851851851851851*Li2sh*sh_m1_pow_m3*
	sh_m2 + 3.5555555555555554*ash*sqrt4sh*sqrtsh*
	(-6. - 13.*sh + 20.*sh_sq - 8.*sh_cube + 
	sh_four) - 0.01646090534979424*sh_m1_pow_m5*
	sh_m6*(-21. + 10.*sh + 
	2.*(-847. + 972.*Lshb)*sh_sq + 
	(11470. - 9648.*Lshb)*sh_cube + 
	(-27667. + 19584.*Lshb - 
	600.*pi_sq)*sh_four + 
	8.*(3982. - 2565.*Lshb + 
	150.*pi_sq)*sh_five + 
	2.*(-9179. + 5256.*Lshb - 
	300.*pi_sq)*sh_six - 
	6.*(-607. + 348.*Lshb)*sh_seven - 
	6.*(7. + 36.*Lshb)*sh_eight + 
	6.*(5. + 18.*Lshb)*sh_nine)) - 
	0.09876543209876543*pow(Lsh,2)*sh_m1_pow_m6*
	sh_m6*(-6. - 14.*sh - 1215.*sh_sq + 
	7695.*sh_cube - 19257.*sh_four - 
	60.*Lshb*pow(-1. + sh,3)*sh_four + 
	25233.*sh_five - 18375.*sh_six + 
	6983.*sh_seven - 1056.*sh_eight - 27.*sh_nine + 
	9.*sh_ten) - 
	0.04938271604938271*lz_sq*sh_m6*
	(-3. - 25.*sh - 2580.*sh_sq + 73.*sh_four + 
	3.*sh_five + 180.*sh_six + 36.*sh_seven - 
	738.*sh_eight + 612.*sh_nine - 180.*sh_ten + 
	18.*sh_eleven) - 
	1.*lz*(128.*L*sh_m4 + 
	5.925925925925926*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 + 7.111111111111111*ash*sqrt4sh*sqrtsh*
	(-6. - 13.*sh + 20.*sh_sq - 8.*sh_cube + 
	sh_four) + 0.19753086419753085*Lsh*
	sh_m1_pow_m5*sh_m6*
	(-3. - 10.*sh - 1240.*sh_sq + 
	6380.*sh_cube - 12722.*sh_four + 
	12478.*sh_five - 6074.*sh_six + 
	1158.*sh_seven + 3.*sh_eight + 3.*sh_nine) + 
	0.0016460905349794238*sh_m1_pow_m4*sh_m6*
	(105. + 55.*sh + 18530.*sh_sq - 
	72330.*sh_cube + 91890.*sh_four - 
	26202.*sh_five - 21156.*sh_six + 
	1154.*sh_seven + 54349.*sh_eight - 
	76131.*sh_nine + 4716.*sh_ten + 
	82440.*sh_eleven - 83970.*sh_twelve + 
	36900.*sh_thirteen - 7740.*sh_fourteen + 
	630.*sh_fifteen)) - 
	0.000027434842249657064*sh_m1_pow_m5*sh_m6*
	(2775. - 5550.*sh - 9.841575e6*sh_sq + 
	4.85214e7*sh_cube - 9.54242e7*sh_four + 
	9.3237052e7*sh_five - 4.3743086e7*sh_six + 
	4.46555e6*sh_seven + 7.670645e6*sh_eight - 
	4.88458e6*sh_nine - 1.3483707e7*sh_ten + 
	3.4549776e7*sh_eleven - 3.627225e7*sh_twelve + 
	2.100825e7*sh_thirteen - 6.9246e6*sh_fourteen + 
	1.21095e6*sh_fifteen - 86850.*sh_sixteen + 
	1200.*pi_sq*
	(-9. - 30.*sh - 2694.*sh_sq + 
	13896.*sh_cube - 27086.*sh_four + 
	25318.*sh_five - 11444.*sh_six + 
	2810.*sh_seven - 1033.*sh_eight - 
	2857.*sh_nine + 10122.*sh_ten - 
	14418.*sh_eleven + 11736.*sh_twelve - 
	5778.*sh_thirteen + 1692.*sh_fourteen - 
	270.*sh_fifteen + 18.*sh_sixteen)))*pow(z,4.);
	else
	res += pow(z,4.)*(-1424.5308091578474 - 
	5728.7251525464435*(1. - sh) - 
	14040.07057463595*msh_p1_pow_2 - 
	27178.850860660903*msh_p1_pow_3 + 
	L*(-117.33333333333333 - 
	597.3333333333334*(1. - sh) - 
	1749.3333333333333*msh_p1_pow_2 - 
	3925.3333333333335*msh_p1_pow_3 - 
	7509.333333333333*msh_p1_pow_4) - 
	46368.872894680535*msh_p1_pow_4 + 
	log(z)*(34.25652635469211 + 
	169.5621944588141*(1. - sh) + 
	865.728317571408*msh_p1_pow_2 + 
	2770.286236159471*msh_p1_pow_3 + 
	L*(-128. - 512.*(1. - sh) - 
	1280.*msh_p1_pow_2 - 
	2560.*msh_p1_pow_3 - 4480.*msh_p1_pow_4
	) + 6293.656266280077*msh_p1_pow_4) + 
	(128.59259259259258 + 
	497.77777777777777*(1. - sh) + 
	1274.962962962963*msh_p1_pow_2 + 
	2588.641975308642*msh_p1_pow_3 + 
	4550.617283950617*msh_p1_pow_4)*pow(log(z),2));}

	if(9<=maxpow){
	if(sh<.900001)
	res += 0.006449987402368355*pi_sq*sh_m7*
	(35. + 330.*sh + 1389.*sh_sq + 3212.*sh_cube + 
	1389.*sh_four + 330.*sh_five + 35.*sh_six)*
	pow(z,4.5);
	else
	res += (427.7872976387514 + 1711.1491905550056*(1. - sh) + 
	4470.377260324953*msh_p1_pow_2 + 
	9518.267372462218*msh_p1_pow_3 + 
	17894.4954414067*msh_p1_pow_4)*pow(z,4.5);}

	if(10<=maxpow){
	if(sh<.900001)
	res += (-1.*L*(352.7111111111111*sh_m5 - 
	341.3333333333333*Lsh*sh_m5) + 
	29.917860700861013*sh_m1_pow_m3*sh_m3 - 
	24.88888888888889*Li3sh*sh_m1_pow_m3*sh_m3 + 
	8.296296296296296*pow(Lsh,3)*sh_m1_pow_m3*
	sh_m3 + 2.3703703703703702*cl2*sh*sqrt4sh*sqrtsh*
	(44. + 108.*sh - 221.*sh_sq + 132.*sh_cube - 
	33.*sh_four + 3.*sh_five) + 
	0.23703703703703705*ash*sh*sqrt4sh*sqrtsh*
	sh_m1_pow_m3*(-496. + 326.*sh + 6147.*sh_sq - 
	16495.*sh_cube + 18680.*sh_four - 
	11251.*sh_five + 3760.*sh_six - 
	658.*sh_seven + 47.*sh_eight) + 
	1.1851851851851851*Li2sh*sh_m1_pow_m6*sh_m5*
	(-72. + 429.*sh - 1073.*sh_sq + 1442.*sh_cube - 
	1092.*sh_four + 435.*sh_five - 47.*sh_six - 
	6.*sh_seven - sh_eight + sh_nine) - 
	7.111111111111111*sh*pow(ash,2)*sh_m1_pow_m4*
	(-24. + 162.*sh + 114.*sh_sq - 1809.*sh_cube + 
	4224.*sh_four - 4957.*sh_five + 
	3423.*sh_six - 1445.*sh_seven + 
	366.*sh_eight - 51.*sh_nine + 3.*sh_ten) - 
	1.*Lsh*(-24.88888888888889*Li2sh*sh_m1_pow_m3*
	sh_m3 - 2.3703703703703702*ash*sh*sqrt4sh*
	sqrtsh*(44. + 108.*sh - 221.*sh_sq + 
	132.*sh_cube - 33.*sh_four + 3.*sh_five) + 
	0.0003292181069958848*sh_m1_pow_m7*sh_m8*
	(486. - 594.*sh - 8634.*sh_sq - 
	132632.*sh_cube + 806792.*sh_four - 
	1.649394e6*sh_five + 
	63000.*pi_sq*pow(-1. + sh,4)*
	sh_five + 1.299724e6*sh_six + 
	289454.*sh_seven - 1.252893e6*sh_eight + 
	839053.*sh_nine - 
	3600.*Lshb*(-1. + sh)*sh_cube*
	(-72. + 429.*sh - 1073.*sh_sq + 
	1442.*sh_cube - 1092.*sh_four + 
	435.*sh_five - 47.*sh_six - 
	6.*sh_seven - sh_eight + sh_nine) - 
	206662.*sh_ten + 15750.*sh_eleven + 
	75.*sh_twelve - 525.*sh_thirteen)) - 
	0.009876543209876543*pow(Lsh,2)*sh_m1_pow_m7*
	sh_m8*(36. + 126.*sh - 84.*sh_sq + 
	37024.*sh_cube - 269605.*sh_four + 
	818292.*sh_five - 
	1260.*Lshb*pow(-1. + sh,4)*sh_five - 
	1.360586e6*sh_six + 1.344224e6*sh_seven - 
	788604.*sh_eight + 252535.*sh_nine - 
	33268.*sh_ten - 300.*sh_eleven - 
	120.*sh_twelve + 60.*sh_thirteen) + 
	0.0049382716049382715*lz_sq*sh_m8*
	(18. + 189.*sh + 903.*sh_sq + 78416.*sh_cube - 
	630.*sh_four - 1383.*sh_five - 
	1839.*sh_six - 18.*sh_seven - 3120.*sh_eight - 
	8160.*sh_nine - 1440.*sh_ten + 
	46560.*sh_eleven - 51720.*sh_twelve + 
	23040.*sh_thirteen - 4680.*sh_fourteen + 
	360.*sh_fifteen) - 
	1.*lz*(341.3333333333333*L*sh_m5 + 
	12.444444444444445*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 - 4.7407407407407405*ash*sh*sqrt4sh*
	sqrtsh*(44. + 108.*sh - 221.*sh_sq + 
	132.*sh_cube - 33.*sh_four + 3.*sh_five) + 
	0.019753086419753086*Lsh*sh_m1_pow_m6*sh_m8*
	(18. + 81.*sh + 39.*sh_sq + 37525.*sh_cube - 
	231198.*sh_four + 582087.*sh_five - 
	774221.*sh_six + 575781.*sh_seven - 
	227535.*sh_eight + 38728.*sh_nine - 
	450.*sh_ten + 90.*sh_eleven + 15.*sh_twelve
	) - 0.000047031158142269256*sh_m1_pow_m5*
	sh_m8*(1701. + 1323.*sh - 29274.*sh_sq - 
	957341.*sh_cube + 5.622085e6*sh_four - 
	1.3372289e7*sh_five + 1.7244458e7*sh_six - 
	1.2810919e7*sh_seven + 6.090011e6*sh_eight - 
	4.03983e6*sh_nine + 4.371768e6*sh_ten + 
	4.551324e6*sh_eleven - 1.1897697e7*sh_twelve - 
	1.323264e7*sh_thirteen + 6.277908e7*sh_fourteen - 
	8.508528e7*sh_fifteen + 6.212682e7*sh_sixteen - 
	2.701314e7*sh_seventeen + 6.98418e6*sh_eighteen - 
	990360.*sh_nineteen + 59220.*sh_twenty)) - 
	1.679684219366759e-7*sh_m1_pow_m6*sh_m8*
	(-161406. + 318843.*sh + 1.910657e6*sh_sq + 
	3.755399935e9*sh_cube - 
	2.2256359629e10*sh_four + 
	5.4715328248e10*sh_five - 
	7.1399970054e10*sh_six + 
	5.2256870574e10*sh_seven - 
	2.1211105332e10*sh_eight + 
	5.981785659e9*sh_nine - 
	2.944201467e9*sh_ten - 
	1.246759469e9*sh_eleven + 
	2.466843039e9*sh_twelve + 
	1.8171607842e10*sh_thirteen - 
	6.053910752e10*sh_fourteen + 
	8.996385496e10*sh_fifteen - 
	7.96292826e10*sh_sixteen + 
	4.509057616e10*sh_seventeen - 
	1.647058952e10*sh_eighteen + 
	3.74875284e9*sh_nineteen - 4.8245988e8*sh_twenty + 
	2.674812e7*sh_twenty_one - 
	58800.*pi_sq*
	(-18. - 81.*sh - 39.*sh_sq - 
	28525.*sh_cube + 176088.*sh_four - 
	439209.*sh_five + 573427.*sh_six - 
	414687.*sh_seven + 157597.*sh_eight - 
	24058.*sh_nine + 3132.*sh_ten - 
	3736.*sh_eleven - 77479.*sh_twelve + 
	318588.*sh_thirteen - 608560.*sh_fourteen + 
	702440.*sh_fifteen - 529920.*sh_sixteen + 
	267160.*sh_seventeen - 89120.*sh_eighteen + 
	18840.*sh_nineteen - 2280.*sh_twenty + 
	120.*sh_twenty_one)))*pow(z,5.);
	else
	res += pow(z,5.)*(-4202.986795121401 - 
	21563.981801901395*(1. - sh) - 
	63115.765734768596*msh_p1_pow_2 - 
	142119.7662710934*msh_p1_pow_3 + 
	L*(-352.7111111111111 - 
	2104.8888888888887*(1. - sh) - 
	7168.*msh_p1_pow_2 - 
	18432.*msh_p1_pow_3 - 
	39850.666666666664*msh_p1_pow_4) - 
	277830.4069154384*msh_p1_pow_4 + 
	log(z)*(251.91238577751216 + 
	956.915478778669*(1. - sh) + 
	4881.515489147685*msh_p1_pow_2 + 
	16699.40413911222*msh_p1_pow_3 + 
	L*(-341.3333333333333 - 
	1706.6666666666667*(1. - sh) - 
	5120.*msh_p1_pow_2 - 
	11946.666666666666*msh_p1_pow_3 - 
	23893.333333333332*msh_p1_pow_4) + 
	41147.84749220197*msh_p1_pow_4) + 
	(377.758024691358 + 1862.716049382716*(1. - sh) + 
	5810.814814814815*msh_p1_pow_2 + 
	13779.160493827161*msh_p1_pow_3 + 
	27666.172839506173*msh_p1_pow_4)*pow(log(z),2));}

	}
	return res;
}

#ifdef CACHE
double F_29re_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_29re(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double F_29im(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
	double Lshb = log(1.-sh);
	double Li2sh = creal(CLi2(sh));
	double Li3sh = creal(CLi3(sh));
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;
	double sh_m6 = sh_m1*sh_m5;
	double sh_m7 = sh_m1*sh_m6;
	double sh_m8 = sh_m1*sh_m7;
	double sh_sq = pow(sh,2.);
	double sh_cube = sh*sh_sq;
	double sh_four = sh*sh_cube;
	double sh_five = sh*sh_four;
	double sh_six = sh*sh_five;
	double sh_seven = sh*sh_six;
	double sh_eight = sh*sh_seven;
	double sh_nine = sh*sh_eight;
	double sh_ten = sh*sh_nine;       
	double sh_eleven = sh*sh_ten;
	double sh_twelve = sh*sh_eleven;
	double sh_thirteen = sh*sh_twelve;
	double lz_sq = pow(lz,2);
	double z_m1 = pow(z,-1.);      
	double z_m2 = z_m1*z_m1;      
	double z_m3 = z_m2*z_m1;    
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double msh_p1_pow_2 = pow(1. - sh,2.);
	double msh_p1_pow_3 = msh_p1_pow_2*(1.-sh);
	double msh_p1_pow_4 = msh_p1_pow_2*msh_p1_pow_2*(1.-sh);
	double sh_m1_pow_m1 = pow(-1. + sh,-1.);
	double sh_m1_pow_m2 = sh_m1_pow_m1*sh_m1_pow_m1;
	double sh_m1_pow_m3 = sh_m1_pow_m1*sh_m1_pow_m2;
	double sh_m1_pow_m4 = sh_m1_pow_m1*sh_m1_pow_m3;
	double sh_m1_pow_m5 = sh_m1_pow_m1*sh_m1_pow_m4;
	double sh_m1_pow_m6 = sh_m1_pow_m1*sh_m1_pow_m5;
	double sh_m1_pow_m7 = sh_m1_pow_m1*sh_m1_pow_m6;
	const double pi_sq = pow(pi,2.);

	double res=0;

	if(sh<0.4)
	res = -22.0299 - 1.241123023640412*L + 0.620562*Lsh - 
	5.58505*lz + 1.86168*sh + 5.58505*lz*sh + 97.2931*z + 
	33.5103*lz*z + 184.792*sh*z + 11.1701*Lsh*sh*z + 
	100.531*lz*sh*z + 11.1701*z*lz_sq + 
	11.1701*sh*z*lz_sq + 8.19141*sh_sq + 
	5.58505*lz*sh_sq + 313.249*z*sh_sq + 
	22.3402*Lsh*z*sh_sq + 156.3815*lz*z*sh_sq + 
	11.1701*z*lz_sq*sh_sq + 10.3323*sh_cube + 
	5.58505*lz*sh_cube + 459.887*z*sh_cube + 
	33.5103*Lsh*z*sh_cube + 208.5085*lz*z*sh_cube + 
	11.1701*z*lz_sq*sh_cube + 
	0.0177303*sh_cube*z_m3 + 
	0.11968*sh_sq*z_m2 - 
	0.0221629*sh_cube*z_m2 + 1.11701*sh*z_m1 - 
	0.372337*sh_sq*z_m1 - 
	0.523045*sh_cube*z_m1 - 160.51*z_sq + 
	11.1701*Lsh*z_sq + 44.68045*lz*z_sq + 
	22.3402*Lsh*lz*z_sq - 396.864*sh*z_sq + 
	44.6804*Lsh*sh*z_sq + 201.062*lz*sh*z_sq + 
	89.361*Lsh*lz*sh*z_sq + 44.6805*lz_sq*z_sq + 
	122.87125*sh*lz_sq*z_sq - 
	652.279*sh_sq*z_sq + 
	134.041*Lsh*sh_sq*z_sq + 
	491.4845*lz*sh_sq*z_sq + 
	201.062*Lsh*lz*sh_sq*z_sq + 
	234.57225*lz_sq*sh_sq*z_sq - 
	890.889*sh_cube*z_sq + 
	297.87*Lsh*sh_cube*z_sq + 
	934.565*lz*sh_cube*z_sq + 
	357.4435*Lsh*lz*sh_cube*z_sq + 
	379.7825*lz_sq*sh_cube*z_sq + 
	142.135*z_cube + 9.92898*Lsh*z_cube - 
	104.2545*lz*z_cube - 29.78695*Lsh*lz*z_cube + 
	382.697*sh*z_cube + 22.3402*Lsh*sh*z_cube - 
	402.124*lz*sh*z_cube - 134.0415*Lsh*lz*sh*z_cube + 
	640.989*sh_sq*z_cube - 
	990.415*lz*sh_sq*z_cube - 
	357.4435*Lsh*lz*sh_sq*z_cube + 
	810.709*sh_cube*z_cube - 
	111.701*Lsh*sh_cube*z_cube - 
	1958.49*lz*sh_cube*z_cube - 
	744.675*Lsh*lz*sh_cube*z_cube;
	else{
	
	if(0<=maxpow){
	if(sh<.900001)
	res += (21.099091401887005*L + 
	0.310280755910103*sh_m1_pow_m2*
	(76. - 161.*sh + 67.*sh_sq) - 
	0.620561511820206*Lsh*sh_m1_pow_m3*
	(-8. + 33.*sh - 51.*sh_sq + 17.*sh_cube));
	else
	res += (17.065441575055665 + 8.22244003161773*(1. - sh) + 
	21.099091401887005*L + 
	3.5992567685571952*msh_p1_pow_2 + 
	2.2133360588254014*msh_p1_pow_3 + 
	1.573566690686951*msh_p1_pow_4);}

	/*if(1<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(2<=maxpow){
	if(sh<.900001)
	res += (-107.41683587308312*sh*sh_m1_pow_m3 + 
	89.36085770210967*Li3sh*sh*sh_m1_pow_m3 - 
	11.170107212763709*Li2sh*(-7. + 5.*sh)*
	sh_m1_pow_m3 - 7.446738141842472*sh*pow(Lsh,3)*
	sh_m1_pow_m3 + 22.340214425527417*pow(Lsh,2)*
	sh_m1_pow_m2 - Lsh*
	(44.680428851054835*Li2sh*sh*sh_m1_pow_m3 + 
	11.170107212763709*sh_m1_pow_m3*
	(5. + Lshb*(-7. + 5.*sh) + 
	sh*(-3. + 2.*pi_sq) - 
	1.*sh_sq)) + 
	1.861684535460618*sh_m1_pow_m3*sh_m1*
	(18. - sh*(60. + 11.*pi_sq) + 
	(72. + 5.*pi_sq)*sh_sq + 
	(-30. + 4.*pi_sq)*sh_cube))*
	z;
	else
	res += (5.5716604023132 - 7.453434743876799*(1. - sh) - 
	13.15061358913166*msh_p1_pow_2 - 
	16.48686093313114*msh_p1_pow_3 - 
	18.724591847033217*msh_p1_pow_4)*z;}

	/*if(3<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(4<=maxpow){
	if(sh<.900001)
	res += (11.170107212763709*Li2sh*(3. + 2.*sh)*sh_m1_pow_m3 - 
	11.170107212763709*pow(Lsh,2)*sh_m1_pow_m2 - 
	134.0412865531645*L*sh_m2 - 
	1.*lz*(78.19075048934596*Lsh*sh_m1_pow_m3 - 
	5.585053606381854*sh_m1_pow_m2*sh_m2*
	(23. - 39.*sh + 30.*sh_sq)) - 
	0.930842267730309*sh_m1_pow_m3*sh_m2*
	(-315. + 900.*sh + 
	(-957. + 6.*pi_sq)*sh_sq + 
	4.*(78. + pi_sq)*sh_cube + 
	24.*sh_four) + 
	5.585053606381854*Lsh*sh_m1_pow_m4*sh_m2*
	(-25. + 87.*sh - (121. + 6.*Lshb)*sh_sq + 
	(67. + 2.*Lshb)*sh_cube + 
	2.*(-7. + 2.*Lshb)*sh_four))*z_sq;
	else
	res += (-255.36106211401477 - 376.7273797882516*(1. - sh) - 
	438.1474554206565*msh_p1_pow_2 - 
	458.2075495484674*msh_p1_pow_3 + 
	L*(-134.0412865531645 - 
	268.082573106329*(1. - sh) - 
	402.1238596594935*msh_p1_pow_2 - 
	536.165146212658*msh_p1_pow_3 - 
	670.2064327658226*msh_p1_pow_4) - 
	446.58202617315146*msh_p1_pow_4 + 
	log(z)*(141.48802469500697 + 
	276.46015351590177*(1. - sh) + 
	408.8259239871517*msh_p1_pow_2 + 
	539.8885152835793*msh_p1_pow_3 + 
	670.2064327658226*msh_p1_pow_4))*z_sq;}

	/*if(5<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(6<=maxpow){
	if(sh<.900001)
	res += (-178.72171540421934*L*sh_m3 - 
	22.340214425527417*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m1 + 7.446738141842472*Li2sh*sh_m1_pow_m4*
	sh_m1*(3. + sh + 2.*sh_sq) + 
	0.620561511820206*Lsh*sh_m1_pow_m5*sh_m4*
	(12. + 422.*sh - 2173.*sh_sq + 
	(4115. - 36.*Lshb)*sh_cube + 
	6.*(-625. + 4.*Lshb)*sh_four - 
	6.*(-267. + 2.*Lshb)*sh_five + 
	24.*(-11. + Lshb)*sh_six) - 
	1.*lz*(7.446738141842472*Lsh*sh_m1_pow_m4*
	sh_m1*(3. + sh + 2.*sh_sq) - 
	0.620561511820206*sh_m1_pow_m3*sh_m4*
	(-6. - 428.*sh + 1275.*sh_sq - 
	1143.*sh_cube + 335.*sh_four + 
	27.*sh_five + 12.*sh_six)) - 
	0.10342691863670102*sh_m1_pow_m4*sh_m4*
	(30. + 1796.*sh - 7538.*sh_sq + 
	36.*(299. + pi_sq)*sh_cube + 
	(-6521. + 12.*pi_sq)*sh_four + 
	8.*(106. + 3.*pi_sq)*sh_five + 
	381.*sh_six + 24.*sh_seven))*z_cube;
	else
	res += (-139.59531208395535 - 127.73017597795302*(1. - sh) + 
	172.3157475693439*msh_p1_pow_2 + 
	862.5554678942696*msh_p1_pow_3 + 
	L*(-178.72171540421934 - 
	536.165146212658*(1. - sh) - 
	1072.330292425316*msh_p1_pow_2 - 
	1787.2171540421934*msh_p1_pow_3 - 
	2680.8257310632903*msh_p1_pow_4) + 
	2025.5930676420699*msh_p1_pow_4 + 
	log(z)*(322.07142463468693 + 
	919.2998236104531*(1. - sh) + 
	1803.5999779542467*msh_p1_pow_2 + 
	2978.9966723284447*msh_p1_pow_3 + 
	4449.37284876415*msh_p1_pow_4))*z_cube;}

	/*if(7<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(8<=maxpow){
	if(sh<.900001)
	res += (-402.1238596594935*L*sh_m4 - 
	37.23369070921236*Li2sh*sh_m1_pow_m3*sh_m2 - 
	37.23369070921236*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m2 - lz*
	(-37.23369070921236*Lsh*sh_m1_pow_m3*sh_m2 - 
	0.620561511820206*sh_m1_pow_m2*sh_m6*
	(3. + 19.*sh + 1126.*sh_sq - 2252.*sh_cube + 
	975.*sh_four + 39.*sh_five + 30.*sh_six))\
	- 0.620561511820206*Lsh*sh_m1_pow_m6*sh_m6*
	(6. + 14.*sh + 1053.*sh_sq - 6729.*sh_cube + 
	(16821. - 60.*Lshb)*sh_four + 
	3.*(-7297. + 60.*Lshb)*sh_five - 
	3.*(-5263. + 60.*Lshb)*sh_six + 
	(-5933. + 60.*Lshb)*sh_seven + 900.*sh_eight) + 
	0.05171345931835051*sh_m1_pow_m5*sh_m6*
	(21. - 10.*sh + 3476.*sh_sq - 18346.*sh_cube + 
	(36523. + 120.*pi_sq)*
	sh_four - 16.*
	(2108. + 15.*pi_sq)*sh_five + 
	2.*(6731. + 60.*pi_sq)*
	sh_six - 96.*sh_seven - 678.*sh_eight - 
	408.*sh_nine + 144.*sh_ten))*pow(z,4.);
	else
	res += (-24.201898960988036 + 587.5698023025075*(1. - sh) + 
	2943.504986434556*msh_p1_pow_2 + 
	8455.573910153013*msh_p1_pow_3 + 
	L*(-402.1238596594935 - 
	1608.495438637974*(1. - sh) - 
	4021.238596594935*msh_p1_pow_2 - 
	8042.47719318987*msh_p1_pow_3 - 
	14074.335088082273*msh_p1_pow_4) + 
	18802.085498788645*msh_p1_pow_4 + 
	log(z)*(819.141195602672 + 
	3200.2357164568025*(1. - sh) + 
	7927.052751991311*msh_p1_pow_2 + 
	15806.942829084288*msh_p1_pow_3 + 
	27672.433541664814*msh_p1_pow_4))*pow(z,4.);}

	/*if(9<=maxpow){
	if(sh<.900001)
	res += 0.;
	else
	res += 0.;}*/

	if(10<=maxpow){
	if(sh<.900001)
	res += (-1072.330292425316*L*sh_m5 - 
	78.19075048934596*Li2sh*sh_m1_pow_m3*sh_m3 - 
	78.19075048934596*pow(Lsh,2)*sh_m1_pow_m3*
	sh_m3 - lz*
	(-78.19075048934596*Lsh*sh_m1_pow_m3*sh_m3 - 
	0.062056151182020604*sh_m1_pow_m2*sh_m8*
	(18. + 153.*sh + 543.*sh_sq + 
	34531.*sh_cube - 69998.*sh_four + 
	32548.*sh_five + 60.*sh_six + 
	690.*sh_seven + 195.*sh_eight)) - 
	0.062056151182020604*Lsh*sh_m1_pow_m7*sh_m8*
	(-36. - 126.*sh + 84.*sh_sq - 32704.*sh_cube + 
	239545.*sh_four + 
	36.*(-20227. + 35.*Lshb)*sh_five + 
	(1.209686e6 - 5040.*Lshb)*sh_six + 
	56.*(-21289. + 135.*Lshb)*sh_seven + 
	(696984. - 5040.*Lshb)*sh_eight + 
	35.*(-6389. + 36.*Lshb)*sh_nine + 
	30808.*sh_ten) - 
	0.0010342691863670101*sh_m1_pow_m6*sh_m8*
	(486. - 108.*sh - 8742.*sh_sq + 
	126466.*sh_cube - 671262.*sh_four + 
	1.573824e6*sh_five - 
	12600.*pi_sq*pow(-1. + sh,3)*
	sh_five - 1.781252e6*sh_six + 
	797802.*sh_seven + 214029.*sh_eight - 
	253718.*sh_nine - 92700.*sh_ten + 
	88650.*sh_eleven - 4275.*sh_twelve - 
	5400.*sh_thirteen))*pow(z,5.);
	else
	res += (351.27268533709133 + 3865.7340081335415*(1. - sh) + 
	17368.68228224328*msh_p1_pow_2 + 
	52713.83107089804*msh_p1_pow_3 + 
	L*(-1072.330292425316 - 
	5361.651462126581*(1. - sh) - 
	16084.95438637974*msh_p1_pow_2 - 
	37531.56023488606*msh_p1_pow_3 - 
	75063.12046977213*msh_p1_pow_4) + 
	127763.72928220684*msh_p1_pow_4 + 
	log(z)*(2503.8415878921674 + 
	12312.560956024709*(1. - sh) + 
	36728.24335783481*msh_p1_pow_2 + 
	85657.34659956669*msh_p1_pow_3 + 
	171742.8806422893*msh_p1_pow_4))*pow(z,5.);}

	}
	return res;
}

#ifdef CACHE
double F_29im_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=F_29im(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double DeltaF_29re(double L, double z, double sh, int maxpow)
{
	double lz = log(z);
	double Lsh = log(sh);
        double sh_sq = pow(sh,2.);
        double sh_cube = sh*sh_sq;
        double lz_sq = pow(lz,2);
        double z_m1 = pow(z,-1.);      
        double z_m2 = z_m1*z_m1;      
        double z_m3 = z_m2*z_m1;   
        double z_sq = pow(z,2.);
        double z_cube = z*z_sq;
        double sh_m1 = pow(sh,-1.);
        double sh_m2 = sh_m1*sh_m1;
        double sh_m3 = sh_m1*sh_m2;
        double sh_m4 = sh_m1*sh_m3;
        double sh_m5 = sh_m1*sh_m4;
        
	double res=0;

	if(sh<0.4)
	res = 0.06772486772486773*(0.6666666666666666 + L - 0.5*lz)*
	(105. + sh_cube*z_m3 + 4.5*sh_sq*z_m2 + 
	21.*sh*z_m1);
	else{
	
	/*if(0<=maxpow)
	res += 0.;

	if(1<=maxpow)
	res += 0.;*/

	if(2<=maxpow)
	res += (-28.444444444444443*sh_m1 - 
	42.666666666666664*L*sh_m1 + 
	21.333333333333332*lz*sh_m1)*z;

	/*if(3<=maxpow)
	res += 0.;*/

	if(4<=maxpow)
	res += (-56.888888888888886*Lsh*sh_m2 - 
	85.33333333333333*L*Lsh*sh_m2 - 
	42.666666666666664*lz_sq*sh_m2 - 
	1.*lz*(-56.888888888888886*sh_m2 - 
	85.33333333333333*L*sh_m2 - 
	42.666666666666664*Lsh*sh_m2))*z_sq;

	/*if(5<=maxpow)
	res += 0.;*/

	if(6<=maxpow)
	res += (113.77777777777777*sh_m3 - 
	113.77777777777777*Lsh*sh_m3 - 
	85.33333333333333*lz_sq*sh_m3 - 
	1.*lz*(-28.444444444444443*sh_m3 - 
	170.66666666666666*L*sh_m3 - 
	85.33333333333333*Lsh*sh_m3) - 
	1.*L*(-170.66666666666666*sh_m3 + 
	170.66666666666666*Lsh*sh_m3))*z_cube;

	/*if(7<=maxpow)
	res += 0.;*/

	if(8<=maxpow)
	res += (398.22222222222223*sh_m4 - 
	341.3333333333333*Lsh*sh_m4 - 
	256.*lz_sq*sh_m4 - 
	1.*lz*(-42.666666666666664*sh_m4 - 
	512.*L*sh_m4 - 256.*Lsh*sh_m4) - 
	1.*L*(-597.3333333333334*sh_m4 + 
	512.*Lsh*sh_m4))*pow(z,4.);

	/*if(9<=maxpow)
	res += 0.;*/

	if(10<=maxpow)
	res += (1403.2592592592594*sh_m5 - 
	1137.7777777777778*Lsh*sh_m5 - 
	853.3333333333334*lz_sq*sh_m5 - 
	1.*lz*(-85.33333333333333*sh_m5 - 
	1706.6666666666667*L*sh_m5 - 
	853.3333333333334*Lsh*sh_m5) - 
	1.*L*(-2104.8888888888887*sh_m5 + 
	1706.6666666666667*Lsh*sh_m5))*pow(z,5.);

	}
	return res;
}

#ifdef CACHE
double DeltaF_29re_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=DeltaF_29re(L,z,sh,maxpow);
		}
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif

/*--------------------------------------------------------------------*/

double DeltaF_29im(double L, double z, double sh, int maxpow)
{
	double lz = log(z); 
	double z_sq = pow(z,2.);
	double z_cube = z*z_sq;
	double sh_m1 = pow(sh,-1.);
	double sh_m2 = sh_m1*sh_m1;
	double sh_m3 = sh_m1*sh_m2;
	double sh_m4 = sh_m1*sh_m3;
	double sh_m5 = sh_m1*sh_m4;

	double res=0;

	if(sh<0.4)
	res = 0.;
	else{
		
	/*if(0<=maxpow)
	res += 0.;

	if(1<=maxpow)
	res += 0.;

	if(2<=maxpow)
	res += 0.;

	if(3<=maxpow)
	res += 0.;*/

	if(4<=maxpow)
	res += (178.72171540421934*sh_m2 + 
	268.082573106329*L*sh_m2 - 
	134.0412865531645*lz*sh_m2)*z_sq;

	/*if(5<=maxpow)
	res += 0.;*/

	if(6<=maxpow)
	res += (357.4434308084387*sh_m3 + 
	536.165146212658*L*sh_m3 - 
	268.082573106329*lz*sh_m3)*z_cube;

	/*if(7<=maxpow)
	res += 0.;*/

	if(8<=maxpow)
	res += (1072.330292425316*sh_m4 + 
	1608.495438637974*L*sh_m4 - 
	804.247719318987*lz*sh_m4)*pow(z,4.);

	/*if(9<=maxpow)
	res += 0.;*/

	if(10<=maxpow)
	res += (3574.434308084387*sh_m5 + 
	5361.651462126581*L*sh_m5 - 
	2680.8257310632903*lz*sh_m5)*pow(z,5.);

	}
	return res;
}

#ifdef CACHE
double DeltaF_29im_cache(double L, double z, double shat, int maxpow)
{
	static double complex Lzstat[NVAL];
	static double table[NVAL][NTAB+1];
	complex double Lz=L+I*z;
	int test=0;
	static int currentmax=-1;
	static int loop=0;
	
	int je=0;
	int max;
	if(loop>0) max=NVAL-1; else max=min(currentmax,NVAL-1);
	
	while(je<=min(currentmax,NVAL-1)&&test==0)
	{
		if((fabs(1.-creal(Lzstat[je])/L)<1.e-4)&&(fabs(1.-cimag(Lzstat[je])/z)<1.e-4)) test=1;
		je++;
	}
	je--;
		
	if(test==0)
	{
		if(currentmax<NVAL-1) currentmax++; else currentmax=0;
		Lzstat[currentmax]=Lz;
		
		int ie;
		for(ie=0;ie<=NTAB;ie++)
		{
			double sh;
			
			if(ie==0)sh=0.00001;
			else if(ie==NTAB)sh=0.99999;
			else sh=(double)ie/NTAB;
			
			table[currentmax][ie]=DeltaF_29im(L,z,sh,maxpow);
		}			
		je=currentmax;			
	}
	
	return interpol_fromtable(shat,table[je],NTAB,3);
}
#endif
