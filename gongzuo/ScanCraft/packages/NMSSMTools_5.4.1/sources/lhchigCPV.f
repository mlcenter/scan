      SUBROUTINE LHC_HIGGS_CPV(PAR,PROB)

***********************************************************************
*   Subroutine to check Higgs LHC constraints:
*      PROB(45) =/= 0: Bound on Br(t->bH+)*BR(H+->tau nu)
*      PROB(46) =/= 0:  No Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0:  chi2gam > chi2max
*      PROB(48) =/= 0:  chi2bb > chi2max
*      PROB(49) =/= 0:  chi2zz > chi2max
*      PROB(51) =/= 0: excluded by ggF/bb->H/A->tautau
*      PROB(52) =/= 0: Excluded H_125->AA->4mu
***********************************************************************

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER I,I1,J,J1,NHGG1,NHGG2,NHGG3,NY
      INTEGER JBAR,JBARbb,K,K1,K2,K3,NX,NA13,NC13
      PARAMETER(NX=18,NA13=15,NC13=28)

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION HMAS(NX),LCMS(NX),LCMSbb(NX)
      DOUBLE PRECISION LATLAS(NX),LATLASbb(NX)
      DOUBLE PRECISION HCMS13(NC13),LCMS13(NC13),LCMSbb13(NC13)
      DOUBLE PRECISION HATLAS13(NA13),LATLAS13(NA13),LATLASbb13(NA13)
      DOUBLE PRECISION MH(5),LCMSH(5),SIG1(5),LATLASH(5)
      DOUBLE PRECISION LCMSHbb(5),SIGbb(5),LATLASHbb(5)
      DOUBLE PRECISION DEL,SIGTOT,MBAR,LCMSMB,LATLASMB
      DOUBLE PRECISION SIGTOTbb,MBARbb,LCMSMBbb,LATLASMBbb
      DOUBLE PRECISION LIMIT,LHC_TBH
      DOUBLE PRECISION D1,D2,CJ2,BRHTOAA,MHcen,masstest
      DOUBLE PRECISION agg,bgg,cgg,mugcengg,muvcengg,SSIG(6)
      DOUBLE PRECISION abb,bbb,cbb,mugcenbb,muvcenbb
      DOUBLE PRECISION azz,bzz,czz,mugcenzz,muvcenzz
      DOUBLE PRECISION BRJJSM,BREESM,BRMMSM,BRLLSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION HGG1(300,2),HGG2(300,2),HGG3(300,2)
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION BRNEU(5,5,5),BRCHA(5,3),BRHSQ(5,10),BRHSL(5,7),
     . BRSUSY(5)
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION SIG(5,8),muH2,ggF8,ggF13,bbH8,bbH13,SHGG,LHGG
      DOUBLE PRECISION C(7),CMIN(7),CMAX(7)

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/LHCSIGCPV/SIG
      COMMON/HIGGSMS/muH2
      COMMON/LHCHGG/HGG1,HGG2,HGG3,NHGG1,NHGG2,NHGG3

***********************************************************************
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

***********************************************************************
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

***********************************************************************

c       I- Constraints from ggF/bb->H/A->tautau

      PROB(51)=0d0

* Loop over 5 Higgses
      DO I=1,5
        MH(I)=dsqrt(MH0(I))

* 8 TeV limits
        J=1
        DOWHILE(HMAS(J).LE.MH(I) .AND. J.LT.NX)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HMAS(NX)) THEN

* ggF Signal cross section*BR:
          SIG1(I)=(CJ(I)**2+CJP(I)**2)*BRLL(I)*ggF8(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMS(J)-LCMS(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSH(I)/SIG1(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLAS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLAS(J)-LATLAS(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLAS(J-1)=LATLAS(7)=.96D0 and not .794d0:
          IF(J.EQ.8) THEN
            LATLASH(I)=.96D0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLAS(J)-.96D0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG1(I))

* bbH Signal cross section*BR:
          SIGbb(I)=(CB(I)**2+CBP(I)**2)*BRLL(I)*bbH8(MH(I))

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
             MBAR=(SIG1(I)*MH(I)+SIG1(I1)*MH(I1))/(SIG1(I)+SIG1(I1))
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
             SIGTOT=SIG1(I)+SIG1(I1)
     .         -SIG1(I)*SIG1(I1)*DEL/(SIG1(I)+SIG1(I1))
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
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .           (LCMSbb(JBARbb)-LCMSbb(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(1d0,LCMSMBbb/SIGTOTbb)

* ATLAS bbH limit at MBARbb:
               LATLASMBbb=LATLASbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .           (LATLASbb(JBARbb)-LATLASbb(JBARbb-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBARbb=8, 
* modif. LATLASbb(JBARbb-1)=LATLASbb(7)=.858D0 and not .393d0:
               IF(J.EQ.8) THEN
                 LATLASMBbb=.858D0+(MBARbb-HMAS(JBARbb-1))/
     .             (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .             (LATLASbb(JBARbb)-.858D0)
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
          SIG1(I)=(CJ(I)**2+CJP(I)**2)*BRLL(I)*ggF13(MH(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLAS13(J-1)+(MH(I)-HATLAS13(J-1))/
     .      (HATLAS13(J)-HATLAS13(J-1))*(LATLAS13(J)-LATLAS13(J-1))
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG1(I))

* bbH Signal cross section*BR:
          SIGbb(I)=(CB(I)**2+CBP(I)**2)**2*BRLL(I)*bbH13(MH(I))

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
             MBAR=(SIG1(I)*MH(I)+SIG1(I1)*MH(I1))/(SIG1(I)+SIG1(I1))
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
             SIGTOT=SIG1(I)+SIG1(I1)
     .         -SIG1(I)*SIG1(I1)*DEL/(SIG1(I)+SIG1(I1))
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
          SIG1(I)=(CJ(I)**2+CJP(I)**2)*BRLL(I)*ggF13(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS13(J-1)+(MH(I)-HCMS13(J-1))/
     .      (HCMS13(J)-HCMS13(J-1))*(LCMS13(J)-LCMS13(J-1))
          PROB(51)=PROB(51)+DDIM(1D0,LCMSH(I)/SIG1(I))

* bbH Signal cross section*BR:
          SIGbb(I)=(CB(I)**2+CBP(I)**2)**2*BRLL(I)*bbH13(MH(I))

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
             MBAR=(SIG1(I)*MH(I)+SIG1(I1)*MH(I1))/(SIG1(I)+SIG1(I1))
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
             SIGTOT=SIG1(I)+SIG1(I1)
     .         -SIG1(I)*SIG1(I1)*DEL/(SIG1(I)+SIG1(I1))
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


c        II- Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/LHC_TBH(dsqrt(MHC)),1d0)


c        III- Constraints from ggF->H/A->gamgam

      PROB(53)=0d0

* Loop over 5 Higgses
      DO I=1,5
        MH(I)=DSQRT(MH0(I))

* CMS 8 TeV limit:
        J=1
        DOWHILE(HGG1(J,1).LE.MH(I) .AND. J.LT.NHGG1)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH(I).LE.HGG1(NHGG1,1)) THEN
          LHGG=HGG1(J-1,2)+(MH(I)-HGG1(J-1,1))/
     .     (HGG1(J,1)-HGG1(J-1,1))*(HGG1(J,2)-HGG1(J-1,2))
          SHGG=(CJ(I)**2+CJP(I)**2)*BRGG(I)*ggF8(MH(I))
          PROB(53)=PROB(53)+DDIM(1d0,LHGG/SHGG)
        ENDIF

* CMS 13 TeV limit:
        J=1
        DOWHILE(HGG2(J,1).LE.MH(I) .AND. J.LT.NHGG2)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH(I).LE.HGG2(NHGG2,1)) THEN
          LHGG=HGG2(J-1,2)+(MH(I)-HGG2(J-1,1))/
     .     (HGG2(J,1)-HGG2(J-1,1))*(HGG2(J,2)-HGG2(J-1,2))
          SHGG=(CJ(I)**2+CJP(I)**2)*BRGG(I)*ggF13(MH(I))
          PROB(53)=PROB(53)+DDIM(1d0,LHGG/SHGG)
        ENDIF

* ATLAS 13 TeV limit:
        J=1
        DOWHILE(HGG3(J,1).LE.MH(I) .AND. J.LT.NHGG3)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH(I).LE.HGG3(NHGG3,1)) THEN
          LHGG=HGG3(J-1,2)+(MH(I)-HGG3(J-1,1))/
     .     (HGG3(J,1)-HGG3(J-1,1))*(HGG3(J,2)-HGG3(J-1,2))
          SHGG=(CJ(I)**2+CJP(I)**2)*BRGG(I)*ggF13(MH(I))
          PROB(53)=PROB(53)+DDIM(1d0,LHGG/SHGG)
        ENDIF
      ENDDO


c       IV- Higgs test at ~125 GeV

c      1) Data

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

c      2) Higgs signals / SM

      DO I=1,5

       DO J=1,8
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(dsqrt(MH0(I)),BRJJSM,BREESM,BRMMSM,BRLLSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VBF/VH
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF
       IF(BRLLSM.NE.0d0)SIG(I,2)=(CJ(I)**2+CJP(I)**2)*BRLL(I)/BRLLSM
       
*   H -> bb
* VBF/VH
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH
       IF(BRBBSM.NE.0d0)SIG(I,4)=(CU(I)**2+CUP(I)**2)*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF/VH
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF
       IF(BRZZSM.NE.0d0)SIG(I,6)=(CJ(I)**2+CJP(I)**2)*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF/VH
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF
       IF(BRGGSM.NE.0d0)SIG(I,8)=(CJ(I)**2+CJP(I)**2)*BRGG(I)/BRGGSM

      ENDDO

c      3) Chi^2 test of the couplings

      MHcen=(MHmin+MHmax)/2d0
      D1=(MHmax-MHmin)/2d0
      masstest=dabs(dsqrt(MH0(1))-MHcen)/D1
      K=0
      K1=0
      K2=0
      K3=0
      J=1

      DO I=1,6
       SSIG(I)=0d0
      ENDDO

      DO I=1,5
      masstest=min(masstest,dabs(dsqrt(MH0(I))-MHcen)/D1)
       IF(dabs(dsqrt(MH0(I))-MHcen).le.D1)THEN
        if(K.eq.2)then
         K3=I
         K=3
        endif
        if(K.eq.1)then
         K2=I
         K=2
        endif
        if(K.eq.0)then
         K1=I
         K=1
        endif
       ELSE
        if(K.eq.0.and.I.gt.1)then
         if(dabs(dsqrt(MH0(I))-MHcen)
     .  .lt.dabs(dsqrt(MH0(I-1))-MHcen))J=I
        endif
       ENDIF
      ENDDO

      IF(K.eq.0)THEN
       chi2gam=agg*(SIG(J,1)-mugcengg)**2+cgg*(SIG(J,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(J,1)-mugcengg)*(SIG(J,2)-muvcengg)
       chi2bb=abb*(SIG(J,3)-mugcenbb)**2+cbb*(SIG(J,4)-muvcenbb)**2
     .    +2d0*bbb*(SIG(J,3)-mugcenbb)*(SIG(J,4)-muvcenbb)
       chi2zz=azz*(SIG(J,5)-mugcenzz)**2+czz*(SIG(J,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(J,5)-mugcenzz)*(SIG(J,6)-muvcenzz)
      ENDIF

      IF(K.eq.1)THEN
       chi2gam=agg*(SIG(K1,1)-mugcengg)**2+cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg)
       chi2bb=abb*(SIG(K1,3)-mugcenbb)**2+cbb*(SIG(K1,4)-muvcenbb)**2
     .    +2d0*bbb*(SIG(K1,3)-mugcenbb)*(SIG(K1,4)-muvcenbb)
       chi2zz=azz*(SIG(K1,5)-mugcenzz)**2+czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz)
      ENDIF

      IF(K.eq.2)THEN
       chi2gam=min(agg*(SIG(K1,1)-mugcengg)**2
     .                +cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg),
     .     agg*(SIG(K2,1)-mugcengg)**2+cgg*(SIG(K2,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K2,1)-mugcengg)*(SIG(K2,2)-muvcengg))
       chi2zz=min(azz*(SIG(K1,5)-mugcenzz)**2
     .                +czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz),
     .      azz*(SIG(K2,5)-mugcenzz)**2+czz*(SIG(K2,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K2,5)-mugcenzz)*(SIG(K2,6)-muvcenzz))

        SSIG(1)=SIG(K1,8)+SIG(K2,8)
        SSIG(2)=SIG(K1,7)+SIG(K2,7)
        SSIG(3)=SIG(K1,2)+SIG(K2,2)
        SSIG(4)=SIG(K1,3)+SIG(K2,3)
        SSIG(5)=SIG(K1,6)+SIG(K2,6)
        SSIG(6)=SIG(K1,5)+SIG(K2,5)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).le.3d0)then
       chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
       chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
        endif

       chi2bb=abb*(SSIG(3)-mugcenbb)**2+cbb*(SSIG(4)-muvcenbb)**2
     .    +2d0*bbb*(SSIG(3)-mugcenbb)*(SSIG(4)-muvcenbb)
      ENDIF

      IF(K.eq.3)THEN
       chi2gam=min(agg*(SIG(K1,1)-mugcengg)**2
     .                 +cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg),
     .     agg*(SIG(K2,1)-mugcengg)**2+cgg*(SIG(K2,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K2,1)-mugcengg)*(SIG(K2,2)-muvcengg),
     .     agg*(SIG(K3,1)-mugcengg)**2+cgg*(SIG(K3,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K3,1)-mugcengg)*(SIG(K3,2)-muvcengg))
       chi2zz=min(azz*(SIG(K1,5)-mugcenzz)**2
     .                 +czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz),
     .      azz*(SIG(K2,5)-mugcenzz)**2+czz*(SIG(K2,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K2,5)-mugcenzz)*(SIG(K2,6)-muvcenzz),
     .      azz*(SIG(K3,5)-mugcenzz)**2+czz*(SIG(K3,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K3,5)-mugcenzz)*(SIG(K3,6)-muvcenzz))

        SSIG(1)=SIG(K1,8)+SIG(K2,8)+SIG(K3,8)
        SSIG(2)=SIG(K1,7)+SIG(K2,7)+SIG(K3,7)
        SSIG(3)=SIG(K1,2)+SIG(K2,2)+SIG(K3,2)
        SSIG(4)=SIG(K1,3)+SIG(K2,3)+SIG(K3,3)
        SSIG(5)=SIG(K1,6)+SIG(K2,6)+SIG(K3,6)
        SSIG(6)=SIG(K1,5)+SIG(K2,5)+SIG(K3,5)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K3))).le.3d0)then
       chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
       chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
        endif

       chi2bb=abb*(SSIG(3)-mugcenbb)**2+cbb*(SSIG(4)-muvcenbb)**2
     .    +2d0*bbb*(SSIG(3)-mugcenbb)*(SSIG(4)-muvcenbb)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K3))).gt.3d0)then
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).lt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).gt.3d0)then
         SSIG(1)=SIG(K1,8)+SIG(K2,8)
         SSIG(2)=SIG(K1,7)+SIG(K2,7)
         SSIG(3)=SIG(K1,2)+SIG(K2,2)
         SSIG(4)=SIG(K1,3)+SIG(K2,3)
         SSIG(5)=SIG(K1,6)+SIG(K2,6)
         SSIG(6)=SIG(K1,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                       +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),
     .     agg*(SIG(K3,1)-mugcengg)**2+cgg*(SIG(K3,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K3,1)-mugcengg)*(SIG(K3,2)-muvcengg))
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                       +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),
     .      azz*(SIG(K3,5)-mugcenzz)**2+czz*(SIG(K3,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K3,5)-mugcenzz)*(SIG(K3,6)-muvcenzz))
         endif
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).gt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).lt.3d0)then
         SSIG(1)=SIG(K3,8)+SIG(K2,8)
         SSIG(2)=SIG(K3,7)+SIG(K2,7)
         SSIG(3)=SIG(K3,2)+SIG(K2,2)
         SSIG(4)=SIG(K3,3)+SIG(K2,3)
         SSIG(5)=SIG(K3,6)+SIG(K2,6)
         SSIG(6)=SIG(K3,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                              +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),
     .     agg*(SIG(K1,1)-mugcengg)**2+cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg))
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                              +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),
     .      azz*(SIG(K1,5)-mugcenzz)**2+czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz))
         endif
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).lt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).lt.3d0)then
         SSIG(1)=SIG(K1,8)+SIG(K2,8)
         SSIG(2)=SIG(K1,7)+SIG(K2,7)
         SSIG(3)=SIG(K1,2)+SIG(K2,2)
         SSIG(4)=SIG(K1,3)+SIG(K2,3)
         SSIG(5)=SIG(K1,6)+SIG(K2,6)
         SSIG(6)=SIG(K1,5)+SIG(K2,5)
         chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
         chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
         SSIG(1)=SIG(K3,8)+SIG(K2,8)
         SSIG(2)=SIG(K3,7)+SIG(K2,7)
         SSIG(3)=SIG(K3,2)+SIG(K2,2)
         SSIG(4)=SIG(K3,3)+SIG(K2,3)
         SSIG(5)=SIG(K3,6)+SIG(K2,6)
         SSIG(6)=SIG(K3,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                  +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),chi2gam)
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                  +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),chi2zz)
         endif
        endif

      ENDIF

      IF(masstest.gt.1d0)PROB(46)=masstest-1d0

c      PROB(47)=DDIM(chi2gam/chi2MAX,1d0)
c      PROB(48)=DDIM(chi2bb/chi2MAX,1d0)
c      PROB(49)=DDIM(chi2zz/chi2MAX,1d0)


c      4) kappas(H_SM) and BR(H_SM -> BSM)

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

      IF(K.EQ.0)THEN
       C(1)=CV(J)
       C(2)=DSQRT(CU(J)**2+CUP(J)**2)
       C(3)=DSQRT(CB(J)**2+CBP(J)**2)
       C(4)=DSQRT(CJ(J)**2+CJP(J)**2)
       C(5)=DSQRT(CG(J)**2+CGP(J)**2)
       C(6)=DSQRT(CL(J)**2+CLP(J)**2)
       C(7)=BRHIGGS(J)+BRSUSY(J)
      ELSE
       C(1)=CV(K1)
       C(2)=DSQRT(CU(K1)**2+CUP(K1)**2)
       C(3)=DSQRT(CB(K1)**2+CBP(K1)**2)
       C(4)=DSQRT(CJ(K1)**2+CJP(K1)**2)
       C(6)=DSQRT(CL(K1)**2+CLP(K1)**2)
       C(5)=DSQRT(CG(K1)**2+CGP(K1)**2)
       C(7)=BRHIGGS(K1)+BRSUSY(K1)
       IF(K.GE.2)THEN
        IF(CV(K2)**2.GT.C(1)**2)THEN
         C(1)=CV(K2)
         C(2)=DSQRT(CU(K2)**2+CUP(K2)**2)
         C(3)=DSQRT(CB(K2)**2+CBP(K2)**2)
         C(4)=DSQRT(CJ(K2)**2+CJP(K2)**2)
         C(5)=DSQRT(CG(K2)**2+CGP(K2)**2)
         C(6)=DSQRT(CL(K2)**2+CLP(K2)**2)
         C(7)=BRHIGGS(K2)+BRSUSY(K2)
	ENDIF
        IF(K.EQ.3 .AND. CV(K3)**2.GT.C(1)**2)THEN
         C(1)=CV(K3)
         C(2)=DSQRT(CU(K3)**2+CUP(K3)**2)
         C(3)=DSQRT(CB(K3)**2+CBP(K3)**2)
         C(4)=DSQRT(CJ(K3)**2+CJP(K3)**2)
         C(5)=DSQRT(CG(K3)**2+CGP(K3)**2)
         C(6)=DSQRT(CL(K3)**2+CLP(K3)**2)
         C(7)=BRHIGGS(K3)+BRSUSY(K3)
	ENDIF
       ENDIF
      ENDIF

      PROB(67)=-DDIM(1d0,C(1)/CMIN(1))
      DO I=2,6
       PROB(66+I)=DDIM(C(I)/CMAX(I),1d0)-DDIM(1d0,C(I)/CMIN(I))
      ENDDO
      PROB(73)=DDIM(C(7)/CMAX(7),1d0)


c        V- Constraints from ggF->H->AA->leptons

      CALL LHC_HSMAA_LEPTONS_CPV(PROB)

c        VI- Constraints from ggF->H->AA->photons

      CALL LHC_HSMAA_PHOTONS_CPV(PROB)

      RETURN

c   Error catch

 1    WRITE(*,*)"Cannot find the file ",FILENAME
      STOP

 2    WRITE(*,*)"Read error in the file ",FILENAME
      STOP

      END


      SUBROUTINE LHC_HSMAA_LEPTONS_CPV(PROB)

*     PROB(52) =/= 0: excluded 

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0
      INTEGER NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5
      INTEGER NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2
      INTEGER NHAATAUB,NHAAGJ,NHAABS
      INTEGER I,J,K

      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION D1(5),D2(5),C2BRHAA,MH,MA
      DOUBLE PRECISION LOWBOUND,HIGHBOUND,LIMIT,ggF8,ggF13

      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION HAATAUS0(100,2),HAATAUS1(300,2),HAATAUS2(100,2)
      DOUBLE PRECISION HAATAUS3	(100,2),HAAMUS0(100,2),HAAMUS1(100,2)
      DOUBLE PRECISION HAAMUS2(100,2),HAAMUS3(100,2),HAAMUS4(100,4)
      DOUBLE PRECISION HAAMUS5(100,4),HAAMUTAU1(100,4),HAAMUTAU2(100,4)
      DOUBLE PRECISION HAAMUB1(100,2),HAAMUB2(200,2),HAATAUB(100,2)
      DOUBLE PRECISION HAAGJ(100,2),HAABS(100,2)

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
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

* HSM -> h1h1 or h2h2
* Determining the SM-like Higgs and its coupling*BR

      DO I=2,5
       D1(I)=DDIM(dsqrt(MH0(I))/MHMAX,1d0)-DDIM(1d0,dsqrt(MH0(I))/MHMIN)
      ENDDO

      DO K=1,2

      MH=0d0
      C2BRHAA=0d0
      MA=dsqrt(MH0(K))
      DO I=K+1,5
       IF(D1(I).EQ.0d0)THEN
        C2BRHAA=C2BRHAA+(CJ(I)**2+CJP(I)**2)*BRHHH(I,2*K-1)
        MH=MH+dsqrt(MH0(I))*(CJ(I)**2+CJP(I)**2)*BRHHH(I,2*K-1)
       ENDIF
      ENDDO

      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

* ggF->HSM->AA->4mu from CMS-PAS-HIG-13-010, 0.25GeV < M_A < 3.55GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS0(1,1) .AND. HAAMUS0(I,1).LE.MA
     .        .AND. I.LT.NHAAMUS0)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS0(NHAAMUS0,1))THEN
       LIMIT=HAAMUS0(I-1,2)+(MA-HAAMUS0(I-1,1))
     .         /(HAAMUS0(I,1)-HAAMUS0(I-1,1))
     .         *(HAAMUS0(I,2)-HAAMUS0(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from CMS-HIG-PAS-15-011, 20GeV < M_A < 62.5GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS0(1,1) .AND. HAATAUS0(I,1).LE.MA
     .        .AND. I.LT.NHAATAUS0)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS0(NHAATAUS0,1))THEN
       LIMIT=HAATAUS0(I-1,2)+(MA-HAATAUS0(I-1,1))
     .         /(HAATAUS0(I,1)-HAATAUS0(I-1,1))
     .         *(HAATAUS0(I,2)-HAATAUS0(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1701.02032 (CMS), 15GeV < M_A < 62.5GeV
* BR(A->2tau) converted in BR(A->2mu), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS1(1,1) .AND. HAAMUS1(I,1).LE.MA
     .        .AND. I.LT.NHAAMUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAAMUS1(NHAAMUS1,1))THEN
       LIMIT=HAAMUS1(I-1,2)+(MA-HAAMUS1(I-1,1))
     .         /(HAAMUS1(I,1)-HAAMUS1(I-1,1))
     .         *(HAAMUS1(I,2)-HAAMUS1(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1701.02032 (CMS), 15GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU1(1,1) .AND. HAAMUTAU1(I,1).LE.MA
     .        .AND. I.LT.NHAAMUTAU1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU1(NHAAMUTAU1,1))THEN
       LIMIT=HAAMUTAU1(I-1,2)+(MA-HAAMUTAU1(I-1,1))
     .         /(HAAMUTAU1(I,1)-HAAMUTAU1(I-1,1))
     .         *(HAAMUTAU1(I,2)-HAAMUTAU1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(K)*BRLL(K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1805.04865 (CMS), 15GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU2(1,1) .AND. HAAMUTAU2(I,1).LE.MA
     .        .AND. I.LT.NHAAMUTAU2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU2(NHAAMUTAU2,1))THEN
       LIMIT=HAAMUTAU2(I-1,2)+(MA-HAAMUTAU2(I-1,1))
     .         /(HAAMUTAU2(I,1)-HAAMUTAU2(I-1,1))
     .         *(HAAMUTAU2(I,2)-HAAMUTAU2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(K)*BRLL(K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1505.01609 (ATLAS), 3.7GeV < M_A < 50GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS1(1,1) .AND. HAATAUS1(I,1).LE.MA
     .        .AND. I.LT.NHAATAUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS1(NHAATAUS1,1))THEN
       LIMIT=HAATAUS1(I-1,2) +(MA-HAATAUS1(I-1,1))
     .         /(HAATAUS1(I,1)-HAATAUS1(I-1,1))
     .         *(HAATAUS1(I,2)-HAATAUS1(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV

      I=1
      DOWHILE(MA.GE.HAATAUS2(1,1) .AND. HAATAUS2(I,1).LE.MA
     .        .AND. I.LT.NHAATAUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS2(NHAATAUS2,1))THEN
       LIMIT=HAATAUS2(I-1,2)+(MA-HAATAUS2(I-1,1))
     .         /(HAATAUS2(I,1)-HAATAUS2(I-1,1))
     .         *(HAATAUS2(I,2)-HAATAUS2(I-1,2))
       PROB(52)=PROB(52)+DDIM(ggF8(MH)*C2BRHAA*BRLL(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV
* BR(A->2tau) converted in BR(A->2mu), Fig. 7 of 1701.02032
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS3(1,1) .AND. HAAMUS3(I,1).LE.MA
     .        .AND. I.LT.NHAAMUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAAMUS3(NHAAMUS3,1))THEN
       LIMIT=HAAMUS3(I-1,2)+(MA-HAAMUS3(I-1,1))
     .         /(HAAMUS3(I,1)-HAAMUS3(I-1,1))
     .         *(HAAMUS3(I,2)-HAAMUS3(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 15GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS3(1,1) .AND. HAATAUS3(I,1).LE.MA
     .        .AND. I.LT.NHAATAUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAATAUS3(NHAATAUS3,1))THEN
       LIMIT=HAATAUS3(I-1,2)+(MA-HAATAUS3(I-1,1))
     .         /(HAATAUS3(I,1)-HAATAUS3(I-1,1))
     .         *(HAATAUS3(I,2)-HAATAUS3(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 15GeV
* BR(A->2tau) converted in BR(A->2mu), Fig. 7 of 1701.02032
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS2(1,1) .AND. HAAMUS2(I,1).LE.MA
     .        .AND. I.LT.NHAAMUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAAMUS2(NHAAMUS2,1))THEN
       LIMIT=HAAMUS2(I-1,2)+(MA-HAAMUS2(I-1,1))
     .         /(HAAMUS2(I,1)-HAAMUS2(I-1,1))
     .         *(HAAMUS2(I,2)-HAAMUS2(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from CMS-PAS-HIG-14-041, 25GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB1(1,1) .AND. HAAMUB1(I,1).LE.MA
     .        .AND. I.LT.NHAAMUB1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAAMUB1(NHAAMUB1,1))THEN
       LIMIT=HAAMUB1(I-1,2)+(MA-HAAMUB1(I-1,1))
     .         /(HAAMUB1(I,1)-HAAMUB1(I-1,1))
     .         *(HAAMUB1(I,2)-HAAMUB1(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(K)*BRBB(K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from 1807.00539 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB2(1,1) .AND. HAAMUB2(I,1).LE.MA
     .        .AND. I.LT.NHAAMUB2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAAMUB2(NHAAMUB2,1))THEN
       LIMIT=HAAMUB2(I-1,2)+(MA-HAAMUB2(I-1,1))
     .         /(HAAMUB2(I,1)-HAAMUB2(I-1,1))
     .         *(HAAMUB2(I,2)-HAAMUB2(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(K)*BRBB(K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2tau2b from 1805.10191 (CMS), 15GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUB(1,1) .AND. HAATAUB(I,1).LE.MA
     .        .AND. I.LT.NHAATAUB)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAATAUB(NHAATAUB,1))THEN
       LIMIT=HAATAUB(I-1,2)+(MA-HAATAUB(I-1,1))
     .         /(HAATAUB(I,1)-HAATAUB(I-1,1))
     .         *(HAATAUB(I,2)-HAATAUB(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(K)*BRBB(K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2gamma2jet from 1803.11145 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAGJ(1,1) .AND. HAAGJ(I,1).LE.MA
     .        .AND. I.LT.NHAAGJ)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAAGJ(NHAAGJ,1))THEN
       LIMIT=HAAGJ(I-1,2)+(MA-HAAGJ(I-1,1))
     .         /(HAAGJ(I,1)-HAAGJ(I-1,1))
     .         *(HAAGJ(I,2)-HAAGJ(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRGG(K)*BRJJ(K)/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* pp->VHSM->AA->4b from 1806.07355 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      DO K=1,2

      MH=0d0
      C2BRHAA=0d0
      MA=dsqrt(MH0(K))
      DO I=K+1,5
       IF(D1(I).EQ.0d0)THEN
        C2BRHAA=C2BRHAA+CV(I)**2*BRHHH(I,2*K-1)
        MH=MH+dsqrt(MH0(I))*CV(I)**2*BRHHH(I,2*K-1)
       ENDIF
      ENDDO

      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

      I=1
      DOWHILE(MA.GE.HAABS(1,1) .AND. HAABS(I,1).LE.MA
     .        .AND. I.LT.NHAABS)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LT.HAABS(NHAABS,1))THEN
       LIMIT=HAABS(I-1,2)+(MA-HAABS(I-1,1))
     .         /(HAABS(I,1)-HAABS(I-1,1))
     .         *(HAABS(I,2)-HAABS(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRBB(K)**2/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* ggF -> H -> a1a1 or h1h1 -> 4mu from 1506.00424 (CMS 8TeV)
* for 0.25GeV < M_A < 3.55GeV, 86GeV < M_A< 150GeV

      DO I=2,5
       D1(I)=DDIM(dsqrt(MH0(I))/HAAMUS4(NHAAMUS4,1),1d0)
     .      -DDIM(1d0,dsqrt(MH0(I))/HAAMUS4(1,1))
      ENDDO

        DO K=1,2

      MA=dsqrt(MH0(K))
      DO I=K+1,5
       IF(D1(I).EQ.0d0)THEN
        D2(I)=(CJ(I)**2+CJP(I)**2)*BRHHH(I,2*K-1)
       ELSE
        D2(I)=0d0
       ENDIF
      ENDDO

       DO J=K+1,5

      MH=0d0
      C2BRHAA=0d0
      DO I=K+1,5
       IF(DABS(dsqrt(MH0(I))-dsqrt(MH0(J))).LT.3d0)THEN
        C2BRHAA=C2BRHAA+D2(I)
	MH=MH+dsqrt(MH0(I))*D2(I)
       ENDIF
      ENDDO
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

      I=1
      DOWHILE(MH.GE.HAAMUS4(1,1)
     .        .AND. HAAMUS4(I,1).LE.MH .AND. I.LT.NHAAMUS4)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LE.HAAMUS4(NHAAMUS4,1))THEN
       LOWBOUND=100d0
       HIGHBOUND=100d0
       IF(MA.GE.0.25d0 .and. MA.LT.2d0)THEN
        LOWBOUND=HAAMUS4(I-1,2)+(MA-0.25d0)/(2d0-0.25d0)
     .         *(HAAMUS4(I-1,3)-HAAMUS4(I-1,2))
        HIGHBOUND=HAAMUS4(I,2)+(MA-0.25d0)/(2d0-0.25d0)
     .         *(HAAMUS4(I,3)-HAAMUS4(I,2))
       ELSEIF(MA.GE.2d0 .and. MA.LE.3.55d0)THEN
        LOWBOUND=HAAMUS4(I-1,3)+(MA-2d0)/(3.55d0-2d0)
     .         *(HAAMUS4(I-1,4)-HAAMUS4(I-1,3))
        HIGHBOUND=HAAMUS4(I,3)+(MA-2d0)/(3.55d0-2d0)
     .         *(HAAMUS4(I,4)-HAAMUS4(I,3))
       ENDIF

       LIMIT=(LOWBOUND
     .    +(MH-HAAMUS4(I-1,1))/(HAAMUS4(I,1)-HAAMUS4(I-1,1))
     .         *(HIGHBOUND-LOWBOUND))/1d3 ! translated in pb
       PROB(52)=PROB(52)+DDIM(ggF8(MH)*C2BRHAA*BRMM(K)**2/LIMIT,1d0)
      ENDIF

      ENDIF
       ENDDO
        ENDDO

* ggF -> H -> a1a1 or h1h1 -> 4mu from CMS-PAS-HIG-18-003 (CMS 13TeV)
* for 0.25GeV < M_A < 3.55GeV, 90GeV < M_A< 150GeV

      DO I=2,5
       D1(I)=DDIM(dsqrt(MH0(I))/150d0,1d0)-DDIM(1d0,dsqrt(MH0(I))/90d0)
      ENDDO

        DO K=1,2

      MA=dsqrt(MH0(K))
      DO I=K+1,5
       IF(D1(I).EQ.0d0)THEN
        D2(I)=(CJ(I)**2+CJP(I)**2)*BRHHH(I,2*K-1)
       ELSE
        D2(I)=0d0
       ENDIF
      ENDDO

       DO J=K+1,5

      MH=0d0
      C2BRHAA=0d0
      DO I=K+1,5
       IF(DABS(dsqrt(MH0(I))-dsqrt(MH0(J))).LT.3d0)THEN
        C2BRHAA=C2BRHAA+D2(I)
	MH=MH+dsqrt(MH0(I))*D2(I)
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
       PROB(52)=PROB(52)+DDIM(ggF13(MH)*C2BRHAA*BRMM(K)**2/LIMIT,1d0)
        ENDIF
       ENDIF

      ENDIF
       ENDDO
        ENDDO

      END


      SUBROUTINE LHC_HSMAA_PHOTONS_CPV(PROB)

*     H_125 -> AA -> photons from ATLAS, arXiv:1509.05051
*     PROB(63) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER I,J,K,N
      PARAMETER (N=185)

      DOUBLE PRECISION PROB(*),S,SMAX,M(N),L(N),MA,MH
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
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

      DO I=1,2
       MA=dsqrt(MH0(I))
       IF(MA.GT.M(1) .AND. MA.LE.M(N))THEN
        S=0d0
        J=1
        DOWHILE(M(J).LT.MA)
         J=J+1
        ENDDO
        SMAX=L(J-1)+(MA-M(J-1))/(M(J)-M(J-1))*(L(J)-L(J-1))
        DO K=I+1,5
         MH=dsqrt(MH0(K))
         IF(MH.GE.MHMIN .AND. MH.LE.MHMAX)THEN
          S=S+(CJ(K)**2+CJP(K)**2)*BRHHH(K,2*I-1)*BRGG(I)**2
         ENDIF
        ENDDO
        PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
       ENDIF
      ENDDO

      END
