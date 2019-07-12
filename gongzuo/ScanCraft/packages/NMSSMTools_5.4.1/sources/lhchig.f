      SUBROUTINE LHCHIG(PAR,PROB)

*   Subroutine to check LHC Higgs constraints
*      PROB(45) =/= 0  excluded by t -> bH+ (LHC)

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PAR(*),PROB(*),SIG(5,8)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BREESM,BRMMSM,BRLLSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),LHC_TBH
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH

      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .       HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .       HCBRSUSY,HCWIDTH
      COMMON/LHCSIG/SIG

* Loop over H1, H2, H3

      DO I=1,3

       DO J=1,8
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(SMASS(I),BRJJSM,BREESM,BRMMSM,BRLLSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VBF/VH
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF
       IF(BRLLSM.NE.0d0)SIG(I,2)=CJ(I)**2*BRLL(I)/BRLLSM
       
*   H -> bb
* VBF/VH
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH
       IF(BRBBSM.NE.0d0)SIG(I,4)=CU(I)**2*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF/VH
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF
       IF(BRZZSM.NE.0d0)SIG(I,6)=CJ(I)**2*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF/VH
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF
       IF(BRGGSM.NE.0d0)SIG(I,8)=CJ(I)**2*BRGG(I)/BRGGSM

      ENDDO

* Loop over A1, A2

      DO I=1,2

       DO J=1,8
        SIG(I+3,J)=0d0
       ENDDO

       CALL HDECAY(PMASS(I),BRJJSM,BREESM,BRMMSM,BRLLSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   A -> tautau
* ggF
       IF(BRLLSM.NE.0d0)SIG(I+3,2)=CJ(I+3)**2*BRLL(I+3)/BRLLSM
       
*   A -> bb
* ttH
       IF(BRBBSM.NE.0d0)SIG(I+3,4)=CU(I+3)**2*BRBB(I+3)/BRBBSM

*   A -> ZZ/WW
* ggF
       IF(BRZZSM.NE.0d0)SIG(I+3,6)=CJ(I+3)**2*BRZZ(I+3)/BRZZSM
       
*   A -> gammagamma
* ggF
       IF(BRGGSM.NE.0d0)SIG(I+3,8)=CJ(I+3)**2*BRGG(I+3)/BRGGSM

      ENDDO

* Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/LHC_TBH(CMASS),1d0)

* Other LHC bounds

      CALL Higgs_CHI2(PAR,PROB)
      CALL Higgs_KAPPAS(PAR,PROB)
      CALL HEAVYHIGGS_TAUTAU(PAR,PROB)
      CALL ATLAS_H_GAMGAM(PROB)
      CALL LHC_HSMAA_LEPTONS(PROB)
      CALL LHC_HSMAA_PHOTONS(PROB)

      END


      DOUBLE PRECISION FUNCTION LHC_TBH(M)

* ATLAS constraints on BR(t->bH+)*BR(H+->taunu), ATLAS-CONF-2011-151 tab.5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=8)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0/ 
      DATA Y/.104d0,.098d0,.095d0,.077d0,.066d0,.071d0,.052d0,.141d0/ 

      LHC_TBH=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_TBH=(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END


      SUBROUTINE Higgs_CHI2(PAR,PROB)

*      If HFLAG=/=0:
*      If HFLAG = 1: H1 is SM-like
*      If HFLAG = 2: H2 is SM-like
*      If SMASS(X) is in the window MHMIN...MHMAX: DX=0
*      If SMASS(X) < MHMIN: DX = SMASS(X)/MHMIN - 1 < 0
*      If SMASS(X) > MHMAX: DX = SMASS(X)/MHMAX - 1 > 0
*
*      PROB(46) =/= 0  no Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0  chi2gam > chi2max
*      PROB(48) =/= 0  chi2bb > chi2max
*      PROB(49) =/= 0  chi2zz > chi2max


      IMPLICIT NONE
      INTEGER GMUFLAG,HFLAG
      DOUBLE PRECISION PAR(*),PROB(*),SIG(5,8),D1,D2
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION agg,bgg,cgg,mugcengg,muvcengg
      DOUBLE PRECISION abb,bbb,cbb,mugcenbb,muvcenbb
      DOUBLE PRECISION azz,bzz,czz,mugcenzz,muvcenzz

      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/LHCSIG/SIG
      COMMON/GMUFLAG/GMUFLAG,HFLAG

* From J. Bernon with run-II data (sep. 2016):

c Chi^2 from gammagamma:
      agg=17.47d0
      bgg=3.17d0
      cgg=6.26d0
      mugcengg=1.18d0
      muvcengg=1.07d0

c Chi^2 from bb/tautau:
      abb=4.81d0
      bbb=2.69d0
      cbb=21.57d0
      mugcenbb=1.27d0
      muvcenbb=0.86d0

c Chi^2 from ZZ/WW:
      azz=36.40d0
      bzz=4.56d0
      czz=8.16d0
      mugcenzz=1.11d0
      muvcenzz=1.37d0

      IF(HFLAG.EQ.2)THEN
       D1=1d99
      ELSE
       D1=DDIM(SMASS(1)/MHMAX,1d0)-DDIM(1d0,SMASS(1)/MHMIN)
      ENDIF
      IF(HFLAG.EQ.1)THEN
       D2=1d99
      ELSE
       D2=DDIM(SMASS(2)/MHMAX,1d0)-DDIM(1d0,SMASS(2)/MHMIN)
      ENDIF

      IF(D1.EQ.0d0 .AND. D2.EQ.0d0)THEN
       IF(DABS(SMASS(1)-SMASS(2)).LE.3d0) THEN
        chi2gam=agg*(SIG(1,8)+SIG(2,8)-mugcengg)**2 
     .     +cgg*(SIG(1,7)+SIG(2,7)-muvcengg)**2
     .     +2d0*bgg*(SIG(1,8)+SIG(2,8)-mugcengg)
     .       *(SIG(1,7)+SIG(2,7)-muvcengg)
        chi2zz=azz*(SIG(1,6)+SIG(2,6)-mugcenzz)**2 
     .    +czz*(SIG(1,5)+SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)+SIG(2,6)-mugcenzz)
     .       *(SIG(1,5)+SIG(2,5)-muvcenzz)
       ELSE
        chi2gam=MIN(
     .    agg*(SIG(1,8)-mugcengg)**2 +cgg*(SIG(1,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(1,8)-mugcengg)*(SIG(1,7)-muvcengg),
     .    agg*(SIG(2,8)-mugcengg)**2+cgg*(SIG(2,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(2,8)-mugcengg)*(SIG(2,7)-muvcengg))
        chi2zz=MIN(
     .    azz*(SIG(1,6)-mugcenzz)**2+czz*(SIG(1,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)-mugcenzz)*(SIG(1,5)-muvcenzz),
     .    azz*(SIG(2,6)-mugcenzz)**2+czz*(SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(2,6)-mugcenzz)*(SIG(2,5)-muvcenzz))
       ENDIF
       chi2bb=abb*(SIG(1,2)+SIG(2,2)-mugcenbb)**2 
     .    +cbb*(SIG(1,3)+SIG(2,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(1,2)+SIG(2,2)-mugcenbb)
     .       *(SIG(1,3)+SIG(2,3)-muvcenbb)
      ELSEIF(D1.EQ.0d0)THEN
       chi2gam=agg*(SIG(1,8)-mugcengg)**2 +cgg*(SIG(1,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(1,8)-mugcengg)*(SIG(1,7)-muvcengg)
       chi2bb=abb*(SIG(1,2)-mugcenbb)**2+cbb*(SIG(1,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(1,2)-mugcenbb)*(SIG(1,3)-muvcenbb)
       chi2zz=azz*(SIG(1,6)-mugcenzz)**2+czz*(SIG(1,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)-mugcenzz)*(SIG(1,5)-muvcenzz)
      ELSEIF(D2.EQ.0d0)THEN
       chi2gam=agg*(SIG(2,8)-mugcengg)**2+cgg*(SIG(2,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(2,8)-mugcengg)*(SIG(2,7)-muvcengg)
       chi2bb=abb*(SIG(2,2)-mugcenbb)**2+cbb*(SIG(2,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(2,2)-mugcenbb)*(SIG(2,3)-muvcenbb)
       chi2zz=azz*(SIG(2,6)-mugcenzz)**2+czz*(SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(2,6)-mugcenzz)*(SIG(2,5)-muvcenzz)
      ELSE
       chi2gam=0d0
       chi2bb=0d0
       chi2zz=0d0
       IF(DABS(D1).LT.DABS(D2))THEN
        PROB(46)=D1
        chi2gam=agg*(SIG(1,8)-mugcengg)**2 +cgg*(SIG(1,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(1,8)-mugcengg)*(SIG(1,7)-muvcengg)
        chi2bb=abb*(SIG(1,2)-mugcenbb)**2+cbb*(SIG(1,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(1,2)-mugcenbb)*(SIG(1,3)-muvcenbb)
        chi2zz=azz*(SIG(1,6)-mugcenzz)**2+czz*(SIG(1,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)-mugcenzz)*(SIG(1,5)-muvcenzz)
       ELSE
        PROB(46)=D2
        chi2gam=agg*(SIG(2,8)-mugcengg)**2+cgg*(SIG(2,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(2,8)-mugcengg)*(SIG(2,7)-muvcengg)
        chi2bb=abb*(SIG(2,2)-mugcenbb)**2+cbb*(SIG(2,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(2,2)-mugcenbb)*(SIG(2,3)-muvcenbb)
        chi2zz=azz*(SIG(2,6)-mugcenzz)**2+czz*(SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(2,6)-mugcenzz)*(SIG(2,5)-muvcenzz)
       ENDIF
      ENDIF

c      PROB(47)=DDIM(chi2gam/chi2MAX,1d0)
c      PROB(48)=DDIM(chi2bb/chi2MAX,1d0)
c      PROB(49)=DDIM(chi2zz/chi2MAX,1d0)

      END


      SUBROUTINE Higgs_KAPPAS(PAR,PROB)

*      If HFLAG=/=0:
*      If HFLAG = 1: H1 is SM-like
*      If HFLAG = 2: H2 is SM-like
*      If SMASS(X) is in the window MHMIN...MHMAX: DX=0
*      If SMASS(X) < MHMIN: DX = SMASS(X)/MHMIN - 1 < 0
*      If SMASS(X) > MHMAX: DX = SMASS(X)/MHMAX - 1 > 0
*
*      CV=C(1), CU=C(2), CB=C(3), CJ=C(4) (2gluons), CG=C(5) (2photons),
*      CL=C(6) (2taus), BR_bsm=C(7) are constrained by the combination of kappas from
*      a) 8 TeV ATLAS+CMS combination, 1606.02266 Table 17 l.h.s
*      b) CMS 1809.10733 Table 8 r.h.s
*      c) ATLAS-CONF-2018-031 Table 9(b)
*      where the H-gamgam and H-gluongluon loops are not resolved, BR_bsm are allowed,
*      and |CV|<=1 is imposed as in the NMSSM.
*      The combinations and 2-sigma uncertainties assume Gaussian uncertainties.
*      For C(7)=BR_bsm the BR_inv and BR_undet from CMS 1809.10733 Table 8 r.h.s are combined.
*
*      PROB(67) =/= 0 k_WZ(H_SM) 2 sigma away from LHC measured value
*      PROB(68) =/= 0 k_top(H_SM) 2 sigma away from LHC measured value
*      PROB(69) =/= 0 k_bot(H_SM) 2 sigma away from LHC measured value
*      PROB(70) =/= 0 k_glu(H_SM) 2 sigma away from LHC measured value
*      PROB(71) =/= 0 k_gam(H_SM) 2 sigma away from LHC measured value
*      PROB(72) =/= 0 k_tau(H_SM) 2 sigma away from LHC measured value
*      PROB(73) =/= 0 BR_bsm(H_SM) 2 sigma away from LHC measured value

      IMPLICIT NONE
      INTEGER I,GMUFLAG,HFLAG
      DOUBLE PRECISION PAR(*),PROB(*),D1,D2,C(7),CMIN(7),CMAX(7)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5),CB(5),CL(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CI(3)

      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL

      IF(HFLAG.EQ.2)THEN
       D1=1d99
      ELSE
       D1=DDIM(SMASS(1)/MHMAX,1d0)-DDIM(1d0,SMASS(1)/MHMIN)
      ENDIF
      IF(HFLAG.EQ.1)THEN
       D2=1d99
      ELSE
       D2=DDIM(SMASS(2)/MHMAX,1d0)-DDIM(1d0,SMASS(2)/MHMIN)
      ENDIF

      CMIN(1)=0.900d0
      CMAX(1)=1.000d0
      CMIN(2)=0.941d0
      CMAX(2)=1.323d0
      CMIN(3)=0.607d0
      CMAX(3)=0.915d0
      CMIN(4)=0.833d0
      CMAX(4)=1.050d0
      CMIN(5)=0.858d0
      CMAX(5)=1.020d0
      CMIN(6)=0.772d0
      CMAX(6)=0.970d0
      CMIN(7)=0.000d0
      CMAX(7)=0.194d0

      IF(D1.EQ.0d0 .AND. D2.EQ.0d0)THEN
       IF(CV(1)**2.GT.CV(2)**2)THEN
        C(1)=CV(1)
        C(2)=CU(1)
        C(3)=CB(1)
        C(4)=CJ(1)
        C(5)=CG(1)
        C(6)=CL(1)
        C(7)=BRHIGGS(1)+BRSUSY(1)
       ELSE
        C(1)=CV(2)
        C(2)=CU(2)
        C(3)=CB(2)
        C(4)=CJ(2)
        C(5)=CG(2)
        C(6)=CL(2)
        C(7)=BRHIGGS(2)+BRSUSY(2)
       ENDIF
      ELSEIF(D1.EQ.0d0)THEN
       C(1)=CV(1)
       C(2)=CU(1)
       C(3)=CB(1)
       C(4)=CJ(1)
       C(5)=CG(1)
       C(6)=CL(1)
       C(7)=BRHIGGS(1)+BRSUSY(1)
      ELSEIF(D2.EQ.0d0)THEN
       C(1)=CV(2)
       C(2)=CU(2)
       C(3)=CB(2)
       C(4)=CJ(2)
       C(5)=CG(2)
       C(6)=CL(2)
       C(7)=BRHIGGS(2)+BRSUSY(2)
      ELSE
       IF(DABS(D1).LT.DABS(D2))THEN
        PROB(46)=D1
        C(1)=CV(1)
        C(2)=CU(1)
        C(3)=CB(1)
        C(4)=CJ(1)
        C(5)=CG(1)
        C(6)=CL(1)
        C(7)=BRHIGGS(1)+BRSUSY(1)
       ELSE
        PROB(46)=D2
        C(1)=CV(2)
        C(2)=CU(2)
        C(3)=CB(2)
        C(4)=CJ(2)
        C(5)=CG(2)
        C(6)=CL(2)
        C(7)=BRHIGGS(2)+BRSUSY(2)
       ENDIF
      ENDIF

      PROB(67)=-DDIM(1d0,C(1)/CMIN(1))
      DO I=2,6
       PROB(66+I)=DDIM(C(I)/CMAX(I),1d0)-DDIM(1d0,C(I)/CMIN(I))
      ENDDO
      PROB(73)=DDIM(C(7)/CMAX(7),1d0)

      END


      SUBROUTINE HEAVYHIGGS_TAUTAU(PAR,PROB)

* Constraints from ggF/bb->H/A->tautau
* PROB(51) =/= 0: excluded by H/A->tautau

      IMPLICIT NONE

      INTEGER I,I1,J,J1,JBAR,JBARbb,NX,NA13,NC13
      PARAMETER(NX=18,NA13=15,NC13=28)

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HMAS(NX),LCMS(NX),LCMSbb(NX)
      DOUBLE PRECISION LATLAS(NX),LATLASbb(NX)
      DOUBLE PRECISION HCMS13(NC13),LCMS13(NC13),LCMSbb13(NC13)
      DOUBLE PRECISION HATLAS13(NA13),LATLAS13(NA13),LATLASbb13(NA13)
      DOUBLE PRECISION MH(5),SIG(5),SIGbb(5),LCMSH(5),LATLASH(5)
      DOUBLE PRECISION LCMSHbb(5),LATLASHbb(5)
      DOUBLE PRECISION DEL,SIGTOT,MBAR,LCMSMB,LATLASMB
      DOUBLE PRECISION SIGTOTbb,MBARbb,LCMSMBbb,LATLASMBbb
      DOUBLE PRECISION ggF8,ggF13,bbH8,bbH13

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH

* 8 TeV limits
      DATA HMAS/90d0,100d0,120d0,140d0,160d0,180d0,200d0,250d0,300d0,
     . 350d0,400d0,450d0,500d0,600d0,700d0,800d0,900d0,1000d0/

* Upper limit on ggF->H->tautau (8 TeV) from CMS-PAS-HIG-13-021, Tab. 7
      DATA LCMS/50.2d0,31.3d0,7.38d0,2.27d0,.845d0,.549d0,.517d0,.315d0,
     . .15d0,.112d0,.103d0,.607d-1,.385d-1,.193d-1,.143d-1,.115d-1,
     . .923d-2,.865d-2/

* Upper limit on bbH->tautau (8 TeV) from CMS-PAS-HIG-13-021, Tab. 8
      DATA LCMSbb/6.03d0,4.14d0,1.76d0,1.25d0,.814d0,.659d0,.553d0,
     . .217d0,
     . .975d-1,.638d-1,.613d-1,.431d-1,.320d-1,.203d-1,.173d-1,.166d-1,
     . .146d-1,.133d-1/

* Upper limit on ggF->H->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLAS/29.1d0,24.0d0,5.25d0,2.02d0,1.39d0,1.00d0,.794d0,
     .  .281d0,.127d0,.112d0,.773d-1,.400d-1,.240d-1,.177d-1,.127d-1,
     .  .993d-2,.840d-2,.735d-2/

* Upper limit on bbH->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASbb/6.32d0,6.32d0,2.73d0,1.27d0,.966d0,.606d0,.393d0,
     .  .305d0,.116d0,.101d0,.656d-1,.363d-1,.238d-1,.159d-1,.117d-1,
     .  .943d-2,.785d-2,.716d-2/

* 13 TeV limits
      DATA HATLAS13/2d2,2.5d2,3d2,3.5d2,4d2,5d2,6d2,7d2,8d2,
     . 1d3,1.2d3,1.5d3,1.75d3,2d3,2.25d3/
      DATA HCMS13/9d1,1d2,1.1d2,1.2d2,1.3d2,1.4d2,1.6d2,1.8d2,2d2,
     . 2.5d2,3.5d2,4d2,4.5d2,5d2,6d2,7d2,8d2,9d2,1d3,1.2d3,1.4d3,
     . 1.6d3,1.8d3,2d3,2.3d3,2.6d3,2.9d3,3.2d3/

* Upper limit on ggF->H->tautau (13 TeV) from ATLAS 1709.07242, Fig. 7(a)
      DATA LATLAS13/7.85d-1,2.50d-1,1.68d-1,2.02d-1,1.52d-1,8.04d-2,
     . 5.29d-2,2.69d-2,1.49d-2,9.61d-3,7.39d-3,5.84d-3,5.79d-3,
     . 6.03d-3,6.59d-3/

* Upper limit on bbH->tautau (13 TeV) from ATLAS 1709.07242, Fig. 7(b)
      DATA LATLASbb13/6.98d-1,1.72d-1,1.24d-1,8.73d-2,6.59d-2,3.22d-2,
     . 1.95d-2,1.12d-2,7.09d-3,4.42d-3,3.89d-3,3.65d-3,3.91d-3,
     . 4.29d-3,4.77d-3/

* Upper limit on ggF->H->tautau (13 TeV) from CMS 1803.06553, Fig. 7(a)
      DATA LCMS13/1.74d1,2.59d1,1.42d1,8.61d0,4.73d0,2.82d0,1.86d0,
     . 1.30d0,8.29d-1,5.56d-1,1.09d-1,7.43d-2,6.25d-2,5.91d-2,5.01d-2,
     . 3.51d-2,2.47d-2,1.73d-2,1.37d-2,8.93d-3,6.42d-3,4.99d-3,4.43d-3,
     . 4.00d-3,3.74d-3,3.64d-3,3.53d-3,3.40d-3/

* Upper limit on bbH->tautau (13 TeV) from CMS 1803.06553, Fig. 7(b)
      DATA LCMSbb13/1.39d1,7.13d0,7.13d0,5.45d0,3.15d0,2.71d0,1.45d0,
     . 1.05d0,7.86d-1,3.61d-1,5.82d-2,4.38d-2,3.59d-2,3.00d-2,1.44d-2,
     . 8.86d-3,6.03d-3,4.74d-3,3.70d-3,2.90d-3,2.84d-3,2.46d-3,2.41d-3,
     . 2.40d-3,2.16d-3,2.35d-3,2.45d-3,2.48e-3/

      PROB(51)=0d0

* Loop over 5 Higgses
      DO I=1,5
        IF(I.LE.3) MH(I)=SMASS(I)
        IF(I.GE.4) MH(I)=PMASS(I-3)

* 8 TeV limits
        J=1
        DOWHILE(HMAS(J).LE.MH(I) .AND. J.LT.NX)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HMAS(NX)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CJ(I)**2*BRLL(I)*ggF8(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMS(J)-LCMS(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSH(I)/SIG(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLAS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLAS(J)-LATLAS(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLAS(J-1)=LATLAS(7)=.96D0 and not .794d0:
          IF(J.EQ.8) THEN
            LATLASH(I)=.96D0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLAS(J)-.96D0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG(I))

* bbH Signal cross section*BR:
          SIGbb(I)=CB(I)**2*BRLL(I)*bbH8(MH(I))

* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMSbb(J)-LCMSbb(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSHbb(I)/SIGbb(I))

* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-LATLASbb(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLASbb(J-1)=LATLASbb(7)=.858D0 and not .393d0:
          IF(J.EQ.8) THEN
            LATLASHbb(I)=.858d0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-.858d0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(1d0,LATLASHbb(I)/SIGbb(I))

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HMAS(J1).LE.MH(I1) .AND. J1.LT.NX)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HMAS(1) .AND. MH(I1).LE.HMAS(NX)) THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HMAS(JBAR).LE.MBAR.AND.JBAR.LT.NX)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HMAS(JBARbb).LE.MBARbb.AND.JBARbb.LT.NX)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 90<MBAR<1000
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0D0.AND.MBAR.GE.HMAS(1).AND.
     .         MBAR.LE.HMAS(NX).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* CMS ggF limit at MBAR:
               LCMSMB=LCMS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LCMS(JBAR)-LCMS(JBAR-1))
               PROB(51)=PROB(51)+DDIM(1d0,LCMSMB/SIGTOT)

* ATLAS ggF limit at MBAR:
               LATLASMB=LATLAS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LATLAS(JBAR)-LATLAS(JBAR-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBAR=8, 
* modif. LATLAS(JBAR-1)=LATLAS(7)=.96D0 and not .794d0:
               IF(J.EQ.8) THEN
                 LATLASMB=.96D0+(MBAR-HMAS(JBAR-1))/
     .            (HMAS(JBAR)-HMAS(JBAR-1))*(LATLAS(JBAR)-.96D0)
               ENDIF
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMB/SIGTOT)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<1000
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HMAS(1).AND.
     .         MBARbb.LE.HMAS(NX).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* CMS bbH limit at MBARbb:
              LCMSMBbb=LCMSbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .         (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .         (LCMSbb(JBARbb)-LCMSbb(JBARbb-1))
              PROB(51)=PROB(51)+DDIM(1d0,LCMSMBbb/SIGTOTbb)

* ATLAS bbH limit at MBARbb:
              LATLASMBbb=LATLASbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .         (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .         (LATLASbb(JBARbb)-LATLASbb(JBARbb-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBARbb=8, 
* modif. LATLASbb(JBARbb-1)=LATLASbb(7)=.858D0 and not .393d0:
              IF(J.EQ.8) THEN
                LATLASMBbb=.858D0+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .           (LATLASbb(JBARbb)-.858D0)
              ENDIF
              PROB(51)=PROB(51)+DDIM(1D0,LATLASMBbb/SIGTOTbb)
             ENDIF
            ENDIF
          ENDDO
        ENDIF


* ATLAS 13 TeV limits
        J=1
        DOWHILE(HATLAS13(J).LE.MH(I) .AND. J.LT.NA13)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HATLAS13(NA13)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CJ(I)**2*BRLL(I)*ggF13(MH(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLAS13(J-1)+(MH(I)-HATLAS13(J-1))/
     .      (HATLAS13(J)-HATLAS13(J-1))*(LATLAS13(J)-LATLAS13(J-1))
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG(I))

* bbH Signal cross section*BR:
          SIGbb(I)=CB(I)**2*BRLL(I)*bbH13(MH(I))

* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb13(J-1)+(MH(I)-HATLAS13(J-1))/
     .      (HATLAS13(J)-HATLAS13(J-1))*(LATLASbb13(J)-LATLASbb13(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LATLASHbb(I)/SIGbb(I))

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HATLAS13(J1).LE.MH(I1) .AND. J1.LT.NA13)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HATLAS13(1) .AND. MH(I1).LE.HATLAS13(NA13))THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HATLAS13(JBAR).LE.MBAR .AND. JBAR.LT.NA13)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HATLAS13(JBARbb).LE.MBARbb .AND. JBARbb.LT.NA13)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 200<MBAR<2250
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0D0.AND.MBAR.GE.HATLAS13(1).AND.
     .         MBAR.LE.HATLAS13(NA13).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* ATLAS ggF limit at MBAR:
               LATLASMB=LATLAS13(JBAR-1)+(MBAR-HATLAS13(JBAR-1))/
     .          (HATLAS13(JBAR)-HATLAS13(JBAR-1))*
     .          (LATLAS13(JBAR)-LATLAS13(JBAR-1))
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMB/SIGTOT)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 200<MBARbb<2250
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HATLAS13(1).AND.
     .         MBARbb.LE.HATLAS13(NA13).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* ATLAS bbH limit at MBARbb:
               LATLASMBbb=LATLASbb13(JBARbb-1)+
     .           (MBARbb-HATLAS13(JBARbb-1))/
     .           (HATLAS13(JBARbb)-HATLAS13(JBARbb-1))*
     .           (LATLASbb13(JBARbb)-LATLASbb13(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMBbb/SIGTOTbb)
             ENDIF
            ENDIF
          ENDDO
        ENDIF


* CMS 13 TeV limits
        J=1
        DOWHILE(HCMS13(J).LE.MH(I) .AND. J.LT.NC13)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HCMS13(NC13)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CJ(I)**2*BRLL(I)*ggF13(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS13(J-1)+(MH(I)-HCMS13(J-1))/
     .      (HCMS13(J)-HCMS13(J-1))*(LCMS13(J)-LCMS13(J-1))
          PROB(51)=PROB(51)+DDIM(1D0,LCMSH(I)/SIG(I))

* bbH Signal cross section*BR:
          SIGbb(I)=CB(I)**2*BRLL(I)*bbH13(MH(I))

* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb13(J-1)+(MH(I)-HCMS13(J-1))/
     .      (HCMS13(J)-HCMS13(J-1))*(LCMSbb13(J)-LCMSbb13(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSHbb(I)/SIGbb(I))

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HCMS13(J1).LE.MH(I1) .AND. J1.LT.NC13)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HCMS13(1) .AND. MH(I1).LE.HCMS13(NC13))THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HCMS13(JBAR).LE.MBAR .AND. JBAR.LT.NC13)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HCMS13(JBARbb).LE.MBARbb .AND. JBARbb.LT.NC13)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 90<MBAR<3200
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0D0.AND.MBAR.GE.HCMS13(1).AND.
     .         MBAR.LE.HCMS13(NC13).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* CMS ggF limit at MBAR:
               LCMSMB=LCMS13(JBAR-1)+(MBAR-HCMS13(JBAR-1))/
     .          (HCMS13(JBAR)-HCMS13(JBAR-1))*
     .          (LCMS13(JBAR)-LCMS13(JBAR-1))
               PROB(51)=PROB(51)+DDIM(1D0,LCMSMB/SIGTOT)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<3200
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HCMS13(1).AND.
     .         MBARbb.LE.HCMS13(NC13).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* CMS bbH limit at MBARbb:
               LCMSMBbb=LCMSbb13(JBARbb-1)+
     .           (MBARbb-HCMS13(JBARbb-1))/
     .           (HCMS13(JBARbb)-HCMS13(JBARbb-1))*
     .           (LCMSbb13(JBARbb)-LCMSbb13(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(1D0,LCMSMBbb/SIGTOTbb)
             ENDIF
            ENDIF
          ENDDO
        ENDIF


      ENDDO

      END


      SUBROUTINE ATLAS_H_GAMGAM(PROB)

* Constraints from ggF->H/A->gamgam
* from ATLAS-CONF-2018-025 & CMS-PAS-HIG-17-013
* PROB(53) =/= 0: excluded by ggF->H/A->gamgam

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER I,J,NHGG1,NHGG2,NHGG3

      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HGG1(300,2),HGG2(300,2),HGG3(300,2)
      DOUBLE PRECISION MH,ggF8,ggF13,SHGG,LHGG
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/LHCHGG/HGG1,HGG2,HGG3,NHGG1,NHGG2,NHGG3

      PROB(53)=0d0

* Loop over 5 Higgses
      DO I=1,5
        IF(I.LE.3) MH=SMASS(I)
        IF(I.GE.4) MH=PMASS(I-3)

* CMS 8 TeV limit:
        J=1
        DOWHILE(HGG1(J,1).LE.MH .AND. J.LT.NHGG1)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG1(NHGG1,1)) THEN
          LHGG=HGG1(J-1,2)+(MH-HGG1(J-1,1))/
     .     (HGG1(J,1)-HGG1(J-1,1))*(HGG1(J,2)-HGG1(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF8(MH)
          PROB(53)=PROB(53)+DDIM(1d0,LHGG/SHGG)
        ENDIF

* CMS 13 TeV limit:
        J=1
        DOWHILE(HGG2(J,1).LE.MH .AND. J.LT.NHGG2)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG2(NHGG2,1)) THEN
          LHGG=HGG2(J-1,2)+(MH-HGG2(J-1,1))/
     .     (HGG2(J,1)-HGG2(J-1,1))*(HGG2(J,2)-HGG2(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF13(MH)
          PROB(53)=PROB(53)+DDIM(1d0,LHGG/SHGG)
        ENDIF

* ATLAS 13 TeV limit:
        J=1
        DOWHILE(HGG3(J,1).LE.MH .AND. J.LT.NHGG3)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG3(NHGG3,1)) THEN
          LHGG=HGG3(J-1,2)+(MH-HGG3(J-1,1))/
     .     (HGG3(J,1)-HGG3(J-1,1))*(HGG3(J,2)-HGG3(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF13(MH)
          PROB(53)=PROB(53)+DDIM(1d0,LHGG/SHGG)
        ENDIF

      ENDDO

      RETURN

*   Error catch

 1    WRITE(*,*)"Cannot find the file ",FILENAME
      STOP

 2    WRITE(*,*)"Read error in the file ",FILENAME
      STOP

      END


      SUBROUTINE LHC_HSMAA_LEPTONS(PROB)

*     PROB(52) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0
      INTEGER NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5
      INTEGER NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2
      INTEGER NHAATAUB,NHAAGJ,NHAABS
      INTEGER I,J,K

      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION HAATAUS0(100,2),HAATAUS1(300,2),HAATAUS2(100,2)
      DOUBLE PRECISION HAATAUS3	(100,2),HAAMUS0(100,2),HAAMUS1(100,2)
      DOUBLE PRECISION HAAMUS2(100,2),HAAMUS3(100,2),HAAMUS4(100,4)
      DOUBLE PRECISION HAAMUS5(100,4),HAAMUTAU1(100,4),HAAMUTAU2(100,4)
      DOUBLE PRECISION HAAMUB1(100,2),HAAMUB2(200,2),HAATAUB(100,2)
      DOUBLE PRECISION HAAGJ(100,2),HAABS(100,2)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5)
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION D1(3),D2(3),C2BRHAA,MH,MA
      DOUBLE PRECISION LOWBOUND,HIGHBOUND,LIMIT,ggF8,ggF13

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/LHCHAA/HAATAUS0,HAATAUS1,HAATAUS2,HAATAUS3,
     .      HAAMUS0,HAAMUS1,HAAMUS2,HAAMUS3,HAAMUS4,HAAMUS5,
     .      HAAMUTAU1,HAAMUTAU2,HAAMUB1,HAAMUB2,HAATAUB,
     .      HAAGJ,HAABS,
     .      NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0,
     .      NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5,
     .      NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2,
     .      NHAATAUB,NHAAGJ,NHAABS

      PROB(52)=0d0

* HSM -> a1a1 or h1h1
* Determining the SM-like Higgs and its coupling*BR

      DO I=1,3
       D1(I)=DDIM(SMASS(I)/MHMAX,1d0)-DDIM(1d0,SMASS(I)/MHMIN)
      ENDDO

      DO K=1,2

      MH=0d0
      C2BRHAA=0d0
      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CJ(I)**2*BRHAA(I,1)
         MH=MH+SMASS(I)*CJ(I)**2*BRHAA(I,1)
	ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CJ(I)**2*BRHHH(I-1)
         MH=MH+SMASS(I)*CJ(I)**2*BRHHH(I-1)
	ENDIF
       ENDDO
      ENDIF
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

* ggF->HSM->AA->4mu from CMS-PAS-HIG-13-010, 0.25GeV < M_A < 3.55GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS0(1,1)
     . .AND. HAAMUS0(I,1).LE.MA .AND. I.LT.NHAAMUS0)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS0(NHAAMUS0,1))THEN
       LIMIT=HAAMUS0(I-1,2)
     .    +(MA-HAAMUS0(I-1,1))/(HAAMUS0(I,1)-HAAMUS0(I-1,1))
     .         *(HAAMUS0(I,2)-HAAMUS0(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from CMS-HIG-PAS-15-011, 20GeV < M_A < 62.5GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS0(1,1)
     .  .AND. HAATAUS0(I,1).LE.MA .AND. I.LT.NHAATAUS0)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS0(NHAATAUS0,1))THEN
       LIMIT=HAATAUS0(I-1,2)
     .    +(MA-HAATAUS0(I-1,1))/(HAATAUS0(I,1)-HAATAUS0(I-1,1))
     .         *(HAATAUS0(I,2)-HAATAUS0(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1701.02032 (CMS), 15GeV < M_A < 62.5GeV
* BR(A->2tau) converted in BR(A->2mu), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS1(1,1)
     .  .AND. HAAMUS1(I,1).LE.MA .AND. I.LT.NHAAMUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS1(NHAAMUS1,1))THEN
       LIMIT=HAAMUS1(I-1,2)
     .    +(MA-HAAMUS1(I-1,1))/(HAAMUS1(I,1)-HAAMUS1(I-1,1))
     .         *(HAAMUS1(I,2)-HAAMUS1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1701.02032 (CMS), 15GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU1(1,1)
     .  .AND. HAAMUTAU1(I,1).LE.MA .AND. I.LT.NHAAMUTAU1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU1(NHAAMUTAU1,1))THEN
       LIMIT=HAAMUTAU1(I-1,2)
     .    +(MA-HAAMUTAU1(I-1,1))/(HAAMUTAU1(I,1)-HAAMUTAU1(I-1,1))
     .         *(HAAMUTAU1(I,2)-HAAMUTAU1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRLL(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1805.04865 (CMS), 15GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU2(1,1)
     .  .AND. HAAMUTAU2(I,1).LE.MA .AND. I.LT.NHAAMUTAU2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU2(NHAAMUTAU2,1))THEN
       LIMIT=HAAMUTAU2(I-1,2)
     .    +(MA-HAAMUTAU2(I-1,1))/(HAAMUTAU2(I,1)-HAAMUTAU2(I-1,1))
     .         *(HAAMUTAU2(I,2)-HAAMUTAU2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRLL(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1505.01609 (ATLAS), 3.7GeV < M_A < 50GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS1(1,1)
     .  .AND. HAATAUS1(I,1).LE.MA .AND. I.LT.NHAATAUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS1(NHAATAUS1,1))THEN
       LIMIT=HAATAUS1(I-1,2)
     .    +(MA-HAATAUS1(I-1,1))/(HAATAUS1(I,1)-HAATAUS1(I-1,1))
     .         *(HAATAUS1(I,2)-HAATAUS1(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV

      I=1
      DOWHILE(MA.GE.HAATAUS2(1,1)
     .  .AND. HAATAUS2(I,1).LE.MA .AND. I.LT.NHAATAUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS2(NHAATAUS2,1))THEN
       LIMIT=HAATAUS2(I-1,2)
     .    +(MA-HAATAUS2(I-1,1))/(HAATAUS2(I,1)-HAATAUS2(I-1,1))
     .         *(HAATAUS2(I,2)-HAATAUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(ggF8(MH)*C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV
* BR(A->2tau) converted in BR(A->2mu), Fig. 7 of 1701.02032
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS3(1,1)
     .  .AND. HAAMUS3(I,1).LE.MA .AND. I.LT.NHAAMUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS3(NHAAMUS3,1))THEN
       LIMIT=HAAMUS3(I-1,2)
     .    +(MA-HAAMUS3(I-1,1))/(HAAMUS3(I,1)-HAAMUS3(I-1,1))
     .         *(HAAMUS3(I,2)-HAAMUS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 15GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS3(1,1)
     .  .AND. HAATAUS3(I,1).LE.MA .AND. I.LT.NHAATAUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS3(NHAATAUS3,1))THEN
       LIMIT=HAATAUS3(I-1,2)
     .    +(MA-HAATAUS3(I-1,1))/(HAATAUS3(I,1)-HAATAUS3(I-1,1))
     .         *(HAATAUS3(I,2)-HAATAUS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 15GeV
* BR(A->2tau) converted in BR(A->2mu), Fig. 7 of 1701.02032
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS2(1,1)
     .  .AND. HAAMUS2(I,1).LE.MA .AND. I.LT.NHAAMUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS2(NHAAMUS2,1))THEN
       LIMIT=HAAMUS2(I-1,2)
     .    +(MA-HAAMUS2(I-1,1))/(HAAMUS2(I,1)-HAAMUS2(I-1,1))
     .         *(HAAMUS2(I,2)-HAAMUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from CMS-PAS-HIG-14-041, 25GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB1(1,1)
     .  .AND. HAAMUB1(I,1).LE.MA .AND. I.LT.NHAAMUB1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUB1(NHAAMUB1,1))THEN
       LIMIT=HAAMUB1(I-1,2)
     .    +(MA-HAAMUB1(I-1,1))/(HAAMUB1(I,1)-HAAMUB1(I-1,1))
     .         *(HAAMUB1(I,2)-HAAMUB1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from 1807.00539 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB2(1,1)
     .  .AND. HAAMUB2(I,1).LE.MA .AND. I.LT.NHAAMUB2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUB2(NHAAMUB2,1))THEN
       LIMIT=HAAMUB2(I-1,2)
     .    +(MA-HAAMUB2(I-1,1))/(HAAMUB2(I,1)-HAAMUB2(I-1,1))
     .         *(HAAMUB2(I,2)-HAAMUB2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2tau2b from 1805.10191 (CMS), 15GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUB(1,1)
     .  .AND. HAATAUB(I,1).LE.MA .AND. I.LT.NHAATAUB)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUB(NHAATAUB,1))THEN
       LIMIT=HAATAUB(I-1,2)
     .    +(MA-HAATAUB(I-1,1))/(HAATAUB(I,1)-HAATAUB(I-1,1))
     .         *(HAATAUB(I,2)-HAATAUB(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRLL(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2gamma2jet from 1803.11145 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAGJ(1,1)
     .  .AND. HAAGJ(I,1).LE.MA .AND. I.LT.NHAAGJ)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAGJ(NHAAGJ,1))THEN
       LIMIT=HAAGJ(I-1,2)
     .    +(MA-HAAGJ(I-1,1))/(HAAGJ(I,1)-HAAGJ(I-1,1))
     .         *(HAAGJ(I,2)-HAAGJ(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRGG(7-3*K)*BRJJ(7-3*K)/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* pp->VHSM->AA->4b from 1806.07355 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      DO K=1,2

      MH=0d0
      C2BRHAA=0d0
      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHAA(I,1)
         MH=MH+SMASS(I)*CV(I)**2*BRHAA(I,1)
	ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHHH(I-1)
         MH=MH+SMASS(I)*CV(I)**2*BRHHH(I-1)
	ENDIF
       ENDDO
      ENDIF
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

      I=1
      DOWHILE(MA.GE.HAABS(1,1)
     .  .AND. HAABS(I,1).LE.MA .AND. I.LT.NHAABS)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAABS(NHAABS,1))THEN
       LIMIT=HAABS(I-1,2)
     .    +(MA-HAABS(I-1,1))/(HAABS(I,1)-HAABS(I-1,1))
     .         *(HAABS(I,2)-HAABS(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRBB(7-3*K)**2/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* ggF -> H -> a1a1 or h1h1 -> 4mu from 1506.00424 (CMS 8TeV)
* for 0.25GeV < M_A < 3.55GeV, 86GeV < M_A< 150GeV

      DO I=1,3
       D1(I)=DDIM(SMASS(I)/HAAMUS4(NHAAMUS4,1),1d0)
     .      -DDIM(1d0,SMASS(I)/HAAMUS4(1,1))
      ENDDO

        DO K=1,2

      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
	 D2(I)=CJ(I)**2*BRHAA(I,1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         D2(I)=CJ(I)**2*BRHHH(I-1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ENDIF

       DO J=K,3

      MH=0d0
      C2BRHAA=0d0
      DO I=K,3
       IF(DABS(SMASS(I)-SMASS(J)).LT.3d0)THEN
        C2BRHAA=C2BRHAA+D2(I)
	MH=MH+SMASS(I)*D2(I)
       ENDIF
      ENDDO
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

       I=1
       DOWHILE(MH.GE.HAAMUS4(1,1)
     .   .AND. HAAMUS4(I,1).LE.MH .AND. I.LT.NHAAMUS4)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.HAAMUS4(NHAAMUS4,1))THEN
        IF(MA.GE.0.25d0 .AND. MA.LE.3.55d0)THEN
         IF(MA.LT.2d0)THEN
          LOWBOUND=HAAMUS4(I-1,2)+(MA-0.25d0)/(2d0-0.25d0)
     .            *(HAAMUS4(I-1,3)-HAAMUS4(I-1,2))
          HIGHBOUND=HAAMUS4(I,2)+(MA-0.25d0)/(2d0-0.25d0)
     .             *(HAAMUS4(I,3)-HAAMUS4(I,2))
         ELSE
          LOWBOUND=HAAMUS4(I-1,3)+(MA-2d0)/(3.55d0-2d0)
     .            *(HAAMUS4(I-1,4)-HAAMUS4(I-1,3))
          HIGHBOUND=HAAMUS4(I,3)+(MA-2d0)/(3.55d0-2d0)
     .             *(HAAMUS4(I,4)-HAAMUS4(I,3))
         ENDIF
         LIMIT=(LOWBOUND
     .        +(MH-HAAMUS4(I-1,1))/(HAAMUS4(I,1)-HAAMUS4(I-1,1))
     .        *(HIGHBOUND-LOWBOUND))/1d3 ! translated in pb
         PROB(52)=PROB(52)
     .           +DDIM(ggF8(MH)*C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
        ENDIF
       ENDIF

      ENDIF
       ENDDO
        ENDDO

* ggF -> H -> a1a1 or h1h1 -> 4mu from CMS-PAS-HIG-18-003 (CMS 13TeV)
* for 0.25GeV < M_A < 3.55GeV, 90GeV < M_A< 150GeV

      DO I=1,3
       D1(I)=DDIM(SMASS(I)/150d0,1d0)-DDIM(1d0,SMASS(I)/90d0)
      ENDDO

        DO K=1,2

      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
	 D2(I)=CJ(I)**2*BRHAA(I,1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         D2(I)=CJ(I)**2*BRHHH(I-1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ENDIF

       DO J=K,3

      MH=0d0
      C2BRHAA=0d0
      DO I=K,3
       IF(DABS(SMASS(I)-SMASS(J)).LT.3d0)THEN
        C2BRHAA=C2BRHAA+D2(I)
	MH=MH+SMASS(I)*D2(I)
       ENDIF
      ENDDO
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

       I=1
       DOWHILE(MA.GE.HAAMUS5(1,1)
     .   .AND. HAAMUS5(I,1).LE.MA .AND. I.LT.NHAAMUS5)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MA.LE.HAAMUS5(NHAAMUS5,1))THEN
        IF(MH.GE.90d0 .AND. MH.LE.150d0)THEN
         IF(MA.LT.125d0)THEN
          LOWBOUND=HAAMUS5(I-1,2)+(MH-90d0)/(125d0-90d0)
     .            *(HAAMUS5(I-1,3)-HAAMUS5(I-1,2))
          HIGHBOUND=HAAMUS5(I,2)+(MH-90d0)/(125d0-90d0)
     .             *(HAAMUS5(I,3)-HAAMUS5(I,2))
         ELSE
          LOWBOUND=HAAMUS5(I-1,3)+(MH-125d0)/(150d0-125d0)
     .            *(HAAMUS5(I-1,4)-HAAMUS5(I-1,3))
          HIGHBOUND=HAAMUS5(I,3)+(MH-2d0)/(150d0-125d0)
     .             *(HAAMUS5(I,4)-HAAMUS5(I,3))
         ENDIF
         LIMIT=(LOWBOUND
     .        +(MH-HAAMUS5(I-1,1))/(HAAMUS5(I,1)-HAAMUS5(I-1,1))
     .        *(HIGHBOUND-LOWBOUND))/1d3 ! translated in pb
         PROB(52)=PROB(52)
     .           +DDIM(ggF13(MH)*C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
        ENDIF
       ENDIF

      ENDIF
       ENDDO
        ENDDO

      END


      SUBROUTINE LHC_HSMAA_PHOTONS(PROB)

*     H_125 -> AA -> photons from ATLAS, arXiv:1509.05051
*     PROB(63) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=185)

      DOUBLE PRECISION PROB(*),S,SMAX,M(N),L(N)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz

      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz

      DATA M/1d1,1.0436d1,1.0744d1,1.0913d1,1.1094d1,1.1260d1,1.1594d1,
     .1.2120d1,1.2521d1,1.2750d1,1.3071d1,1.3438d1,1.3974d1,1.4343d1,
     .1.4674d1,1.5156d1,1.5500d1,1.5746d1,1.5950d1,1.6175d1,1.6394d1,
     .1.6727d1,1.7276d1,1.7734d1,1.8192d1,1.8765d1,1.9224d1,1.9637d1,
     .2.0140d1,2.0713d1,2.1101d1,2.1343d1,2.1479d1,2.1660d1,2.1796d1,
     .2.1946d1,2.2139d1,2.2325d1,2.2606d1,2.3004d1,2.3483d1,2.3861d1,
     .2.4098d1,2.4313d1,2.4513d1,2.4769d1,2.5011d1,2.5410d1,2.5945d1,
     .2.6327d1,2.6785d1,2.7051d1,2.7268d1,2.7520d1,2.7759d1,2.8236d1,
     .2.8561d1,2.8788d1,2.8902d1,2.9013d1,2.9125d1,2.9236d1,2.9355d1,
     .2.9568d1,2.9775d1,3.0112d1,3.0555d1,3.0804d1,3.1182d1,3.1697d1,
     .3.2282d1,3.2715d1,3.3201d1,3.3659d1,3.4173d1,3.4498d1,3.4710d1,
     .3.4989d1,3.5329d1,3.5802d1,3.6004d1,3.6211d1,3.6408d1,3.6568d1,
     .3.6797d1,3.7268d1,3.7594d1,3.7949d1,3.8185d1,3.8414d1,3.8643d1,
     .3.8865d1,3.9083d1,3.9331d1,3.9554d1,3.9789d1,4.0114d1,4.0625d1,
     .4.0935d1,4.1164d1,4.1378d1,4.1529d1,4.1729d1,4.1867d1,4.2046d1,
     .4.2188d1,4.2403d1,4.2651d1,4.3054d1,4.3492d1,4.3797d1,4.3930d1,
     .4.4077d1,4.4237d1,4.4360d1,4.4473d1,4.4584d1,4.4707d1,4.4820d1,
     .4.5059d1,4.5460d1,4.5976d1,4.6205d1,4.6320d1,4.6434d1,4.6546d1,
     .4.6648d1,4.6720d1,4.6870d1,4.6964d1,4.7110d1,4.7215d1,4.7444d1,
     .4.7934d1,4.8249d1,4.8470d1,4.8611d1,4.8823d1,4.9179d1,4.9799d1,
     .4.9974d1,5.0207d1,5.0437d1,5.0663d1,5.0900d1,5.1126d1,5.1455d1,
     .5.1876d1,5.2394d1,5.2604d1,5.2839d1,5.3058d1,5.3289d1,5.3413d1,
     .5.3535d1,5.3658d1,5.3881d1,5.4110d1,5.4485d1,5.5027d1,5.5371d1,
     .5.5600d1,5.5810d1,5.5922d1,5.6031d1,5.6139d1,5.6282d1,5.6418d1,
     .5.6613d1,5.6757d1,5.6959d1,5.7204d1,5.7547d1,5.7948d1,5.8521d1,
     .5.8922d1,5.9180d1,5.9467d1,5.9713d1,6.0055d1,6.0399d1,6.0641d1,
     .6.1099d1,6.1638d1,6.2d1/
      DATA L/3.0495d-4,3.3188d-4,3.5175d-4,3.7958d-4,4.1748d-4,
     .4.5015d-4,4.9288d-4,5.0986d-4,4.7180d-4,4.3749d-4,4.0261d-4,
     .3.8321d-4,3.9233d-4,4.1524d-4,4.5016d-4,4.6345d-4,4.3712d-4,
     .4.0954d-4,3.7648d-4,3.4511d-4,3.2145d-4,3.0275d-4,3.0490d-4,
     .3.1516d-4,3.2702d-4,3.3475d-4,3.2452d-4,3.1097d-4,2.9572d-4,
     .2.9321d-4,3.1935d-4,3.4496d-4,3.7242d-4,4.0961d-4,4.4202d-4,
     .4.7699d-4,5.2462d-4,5.7701d-4,6.1537d-4,6.5265d-4,6.3860d-4,
     .6.0308d-4,5.6361d-4,5.2482d-4,4.8616d-4,4.4202d-4,4.1415d-4,
     .3.8329d-4,3.8954d-4,4.1103d-4,4.0633d-4,3.7242d-4,3.4511d-4,
     .3.1981d-4,2.9690d-4,2.8802d-4,3.1981d-4,3.5175d-4,3.7958d-4,
     .4.0961d-4,4.4202d-4,4.7699d-4,5.1473d-4,5.6613d-4,6.2266d-4,
     .6.5711d-4,6.3463d-4,5.8810d-4,5.4809d-4,5.3743d-4,5.5929d-4,
     .5.8265d-4,5.7396d-4,5.6264d-4,5.8253d-4,6.3463d-4,6.8484d-4,
     .7.2922d-4,7.7615d-4,7.3902d-4,6.7193d-4,6.1092d-4,5.6613d-4,
     .5.2462d-4,4.8002d-4,4.5880d-4,4.9458d-4,5.4330d-4,5.8161d-4,
     .6.2988d-4,6.8101d-4,7.3839d-4,7.9749d-4,8.7313d-4,9.6472d-4,
     .1.0645d-3,1.1396d-3,1.1342d-3,1.0634d-3,9.8790d-4,9.1117d-4,
     .8.4436d-4,7.6770d-4,7.1142d-4,6.4682d-4,5.9940d-4,5.4498d-4,
     .5.0833d-4,4.8501d-4,5.0502d-4,5.5516d-4,5.9940d-4,6.4682d-4,
     .7.1142d-4,7.6770d-4,8.2844d-4,8.9399d-4,9.6472d-4,1.0410d-3,
     .1.1140d-3,1.1812d-3,1.1037d-3,1.0042d-3,9.2868d-4,8.6059d-4,
     .7.9749d-4,7.3902d-4,6.8484d-4,6.2266d-4,5.7701d-4,5.2462d-4,
     .4.8616d-4,4.4428d-4,4.3512d-4,4.7699d-4,5.2462d-4,5.6409d-4,
     .6.2243d-4,6.8064d-4,6.8484d-4,6.3463d-4,5.7701d-4,5.2462d-4,
     .4.7699d-4,4.3556d-4,4.0296d-4,3.7598d-4,3.6342d-4,3.8221d-4,
     .4.0961d-4,4.4202d-4,4.8616d-4,5.3411d-4,5.7701d-4,6.2266d-4,
     .6.7193d-4,7.4853d-4,8.1282d-4,8.8536d-4,8.6768d-4,7.9209d-4,
     .7.2098d-4,6.5926d-4,6.1092d-4,5.6613d-4,5.2462d-4,4.7699d-4,
     .4.4202d-4,4.0188d-4,3.7242d-4,3.4511d-4,3.1805d-4,2.9313d-4,
     .2.7738d-4,2.8251d-4,2.9963d-4,3.1981d-4,3.5175d-4,3.7911d-4,
     .4.1656d-4,4.5876d-4,4.9140d-4,5.0985d-4,4.9083d-4,4.6356d-4/

      S=0d0
      DO I=1,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*BRHAA(I,1)*BRGG(4)**2
       ENDIF
      ENDDO

       I=1
       DOWHILE(M(I).LE.PMASS(1) .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. PMASS(1).LE.M(N))THEN
        SMAX=L(I-1)+(PMASS(1)-M(I-1))/(M(I)-M(I-1))*(L(I)-L(I-1))
        PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
       ENDIF

      S=0d0
      DO I=2,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*BRHHH(I-1)*BRGG(1)**2
       ENDIF
      ENDDO

       I=1
       DOWHILE(M(I).LE.SMASS(1) .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. SMASS(1).LE.M(N))THEN
        SMAX=L(I-1)+(SMASS(1)-M(I-1))/(M(I)-M(I-1))*(L(I)-L(I-1))
        PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
       ENDIF

      END


      DOUBLE PRECISION FUNCTION ggF8(MH)

* SM Higgs ggF production cross section at 8 TeV from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt8TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/4.614d3,2.701d3,1.261d3,6.676d2,4.009d2,2.664d2,1.911d2,
     . 1.450d2,1.147d2,9.365d1,7.827d1,6.661d1,5.752d1,5.025d1,4.431d1,
     . 3.940d1,3.526d1,3.175d1,2.873d1,2.611d1,2.383d1,2.184d1,2.008d1,
     . 1.851d1,1.712d1,1.587d1,1.475d1,1.375d1,1.284d1,1.127d1,9.961d0,
     . 8.866d0,7.945d0,7.163d0,6.492d0,5.919d0,5.426d0,5.000d0,4.632d0,
     . 4.314d0,4.040d0,3.807d0,3.608d0,3.444d0,3.314d0,3.222d0,3.175d0,
     . 3.200d0,3.409d0,3.503d0,3.470d0,3.357d0,3.198d0,3.011d0,2.812d0,
     . 2.609d0,2.411d0,2.220d0,2.038d0,1.868d0,1.709d0,1.562d0,1.426d0,
     . 1.301d0,8.236d-1,5.241d-1,3.377d-1,2.207d-1,1.464d-1,9.849d-2,
     . 6.713d-2,4.634d-2,3.237d-2,2.285d-2,1.628d-2,1.172d-2,8.507d-3,
     . 6.222d-3,4.587d-3,3.402d-3,2.540d-3,1.909d-3,1.442d-3,1.095d-3,
     . 8.349d-4,6.400d-4,4.925d-4,3.805d-4,2.951d-4,2.297d-4,1.793d-4,
     . 1.404d-4,1.103d-4,8.685d-5,6.857d-5,5.422d-5,4.311d-5,3.419d-5,
     . 2.722d-5,2.173d-5,1.737d-5,1.390d-5,1.113d-5,8.936d-6,7.190d-6,
     . 5.771d-6,4.665d-6,3.753d-6,3.028d-6,2.455d-6,1.979d-6,1.604d-6,
     . 1.298d-6,1.048d-6/

      ggF8=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        ggF8=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION ggF13(MH)

* SM Higgs ggF production cross section at 13 TeV from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=154)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0,3050d0,3100d0,3150d0,
     . 3200d0,3250d0,3300d0,3350d0,3400d0,3450d0,3500d0,3550d0,3600d0,
     . 3650d0,3700d0,3750d0,3800d0,3850d0,3900d0,3950d0,4000d0,4050d0,
     . 4100d0,4150d0,4200d0,4250d0,4300d0,4350d0,4400d0,4450d0,4500d0,
     . 4550d0,4600d0,4650d0,4700d0,4750d0,4800d0,4850d0,4900d0,4950d0,
     . 5000d0/

      DATA X/1.900d3,1.203d3,8.458d2,6.322d2,4.923d2,3.949d2,3.240d2,
     . 2.706d2,2.294d2,1.968d2,1.706d2,1.492d2,1.315d2,1.166d2,1.041d2,
     . 9.340d1,8.420d1,7.630d1,6.930d1,6.320d1,5.790d1,5.310d1,4.890d1,
     . 4.520d1,4.180d1,3.880d1,3.600d1,3.350d1,3.129d1,2.737d1,2.409d1,
     . 2.132d1,1.896d1,1.694d1,1.520d1,1.369d1,1.237d1,1.122d1,1.020d1,
     . 9.300d0,8.510d0,7.800d0,7.160d0,6.590d0,6.080d0,5.620d0,5.200d0,
     . 4.820d0,4.480d0,4.160d0,3.880d0,3.620d0,3.380d0,3.160d0,2.960d0,
     . 2.770d0,2.600d0,2.440d0,2.300d0,2.160d0,2.030d0,1.920d0,1.810d0,
     . 1.709d0,1.297d0,1.001d0,7.834d-1,6.206d-1,4.969d-1,4.015d-1,
     . 3.271d-1,2.685d-1,2.219d-1,1.845d-1,1.542d-1,1.295d-1,1.093d-1,
     . 9.260d-2,7.880d-2,6.730d-2,5.760d-2,4.950d-2,4.270d-2,3.690d-2,
     . 3.190d-2,2.770d-2,2.410d-2,2.100d-2,1.838d-2,1.609d-2,1.411d-2,
     . 1.239d-2,1.090d-2,9.600d-3,8.470d-3,7.480d-3,6.620d-3,5.860d-3,
     . 5.190d-3,4.610d-3,4.090d-3,3.640d-3,3.240d-3,2.890d-3,2.570d-3,
     . 2.300d-3,2.050d-3,1.830d-3,1.640d-3,1.460d-3,1.310d-3,1.170d-3,
     . 1.050d-3,9.400d-4,8.470d-4,7.600d-4,6.820d-4,6.130d-4,5.510d-4,
     . 4.950d-4,4.450d-4,4.000d-4,3.600d-4,3.240d-4,2.920d-4,2.630d-4,
     . 2.370d-4,2.130d-4,1.920d-4,1.730d-4,1.560d-4,1.410d-4,1.270d-4,
     . 1.140d-4,1.030d-4,9.300d-5,8.400d-5,7.600d-5,6.800d-5,6.100d-5,
     . 5.500d-5,5.000d-5,4.500d-5,4.100d-5,3.700d-5,3.300d-5,3.000d-5,
     . 2.700d-5,2.400d-5,2.200d-5,2.000d-5,1.800d-5,1.600d-5,1.400d-5/

      ggF13=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        ggF13=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION bbH8(MH)

* SM Higgs bbH production cross section at 8 TeV from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt8TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/6.366d1,3.458d1,2.171d1,1.450d1,9.894d0,7.077d0,5.161d0,
     . 3.900d0,2.991d0,2.351d0,1.871d0,1.510d0,1.228d0,1.012d0,8.349d-1,
     . 6.966d-1,5.842d-1,4.955d-1,4.215d-1,3.595d-1,3.094d-1,2.665d-1,
     . 2.315d-1,2.021d-1,1.763d-1,1.550d-1,1.364d-1,1.203d-1,1.066d-1,
     . 8.430d-2,6.734d-2,5.441d-2,4.434d-2,3.636d-2,3.004d-2,2.494d-2,
     . 2.088d-2,1.755d-2,1.484d-2,1.262d-2,1.076d-2,9.229d-3,7.945d-3,
     . 6.846d-3,5.935d-3,5.156d-3,4.494d-3,3.937d-3,3.451d-3,3.033d-3,
     . 2.674d-3,2.361d-3,2.095d-3,1.860d-3,1.653d-3,1.474d-3,1.317d-3,
     . 1.177d-3,1.055d-3,9.463d-4,8.517d-4,7.668d-4,6.913d-4,6.256d-4,
     . 3.832d-4,2.418d-4,1.567d-4,1.037d-4,6.987d-5,4.789d-5,3.326d-5,
     . 2.346d-5,1.672d-5,1.205d-5,8.756d-6,6.423d-6,4.747d-6,3.537d-6,
     . 2.653d-6,2.000d-6,1.517d-6,1.158d-6,8.871d-7,6.828d-7,5.273d-7,
     . 4.096d-7,3.188d-7,2.489d-7,1.949d-7,1.533d-7,1.200d-7,9.539d-8,
     . 7.518d-8,5.992d-8,4.764d-8,3.804d-8,3.028d-8,2.407d-8,1.937d-8,
     . 1.551d-8,1.245d-8,1.001d-8,8.063d-9,6.396d-9,5.217d-9,4.206d-9,
     . 3.421d-9,2.742d-9,2.215d-9,1.798d-9,1.448d-9,1.171d-9,9.427E-10,
     . 7.647E-10/

      bbH8=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        bbH8=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION bbH13(MH)

* SM Higgs bbH production cross section at 13 TeV from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/1.138d2,6.297d1,4.011d1,2.729d1,1.931d1,1.388d1,1.029d1,
     . 7.846d0,6.102d0,4.849d0,3.914d0,3.189d0,2.633d0,2.194d0,1.838d0,
     . 1.549d0,1.317d0,1.126d0,9.669d-1,8.379d-1,7.262d-1,6.325d-1,
     . 5.534d-1,4.879d-1,4.304d-1,3.818d-1,3.383d-1,3.018d-1,2.693d-1,
     . 2.175d-1,1.769d-1,1.451d-1,1.200d-1,1.000d-1,8.397d-2,7.092d-2,
     . 6.021d-2,5.133d-2,4.410d-2,3.799d-2,3.287d-2,2.854d-2,2.491d-2,
     . 2.180d-2,1.915d-2,1.689d-2,1.491d-2,1.321d-2,1.174d-2,1.045d-2,
     . 9.324d-3,8.351d-3,7.492d-3,6.731d-3,6.046d-3,5.470d-3,4.941d-3,
     . 4.472d-3,4.057d-3,3.690d-3,3.352d-3,3.055d-3,2.784d-3,2.547d-3,
     . 1.651d-3,1.101d-3,7.518d-4,5.251d-4,3.723d-4,2.692d-4,1.967d-4,
     . 1.457d-4,1.091d-4,8.258d-5,6.299d-5,4.856d-5,3.773d-5,2.950d-5,
     . 2.322d-5,1.863d-5,1.465d-5,1.172d-5,9.329d-6,7.618d-6,6.179d-6,
     . 5.038d-6,4.120d-6,3.380d-6,2.787d-6,2.303d-6,1.910d-6,1.586d-6,
     . 1.332d-6,1.105d-6,9.255d-7,7.769d-7,6.532d-7,5.511d-7,4.742d-7,
     . 3.941d-7,3.343d-7,2.837d-7,2.415d-7,2.059d-7,1.748d-7,1.466d-7,
     . 1.265d-7,1.090d-7,9.444d-8,8.022d-8,6.968d-8,5.972d-8,5.178d-8,
     . 4.467d-8/

      bbH13=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        bbH13=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END
