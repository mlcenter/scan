#!/usr/bin/env/python3
import sys,os,re,copy,shutil,subprocess,random,math
sys.path.append('/home/...')

from command.scan.scan import scan
from command.NMSSMTools import NMSSMTools
from command.Experiments.directdetection import DirectDetection
from command.chisqure import *
mh=[125.09,3.,3.]
from command.format.parameter_type import *
from command.outputfile import *


ism='h2'
target_number=2000
step_factor=.5
slop_factor=1.
ignore=[ 'Landau pole'
        ,'relic density'
        ,'Relic density'
        ,'Excluded by LUX'
        ,'b->s gamma'
        ,'B_s->mu+mu-'
        ,'Muon magn'
        ,'No Higgs in the'
        ,'b->c tau nu'
        ]


free=scan()
free.AddScalar('tanB','MINPAR',3,2.,60.)
free.AddScalar('M1',  'EXTPAR',1,20.,800.)
free.AddScalar('M2',  'EXTPAR',2,100.,1200.)
free.AddScalar('Atop','EXTPAR',11,-5e3,5e3)
free.AddFollower('Abottom', 'EXTPAR',12,'Atop')
free.AddScalar('Atau','EXPTPAR',13,100.,1200.)
free.AddFollower('MtauL','EXPTPAR',33,'Atau')
free.AddFollower('MtauR','EXPTPAR',36,'Atau')
free.AddScalar('MQ3L',  'EXTPAR',43,100.,2.e3)
free.AddScalar('MQ3R',  'EXTPAR',46,100.,2.e3)
free.AddFollower('MbottomR', 'EXTPAR',12,'MtopR')
free.AddScalar('Lambda','EXTPAR',61,1e-3,.75,prior_distribution='lognormal')
free.AddScalar('minusK','auxiliary',62,1e-3,.75,prior_distribution='lognormal')
Minus=lambda x :-x
free.AddDependent('kappa','EXTPAR',62,fun=Minus,variables=['muinusK'])
free.AddScalar('A_kappa','EXTPAR',64,-2.e3,2.e3)
free.AddScalar('mu_eff','EXTPAR',65,100.,1000.)
free.AddScalar('MA','EXTPAR',124,50.,2.e3)


L_Nsd=DirectDetection('PandaX_Nsd_2016.txt')
L_Psd=DirectDetection('PandaX_Psd_2016.txt')
L_Psi=DirectDetection('LUX201608_Psi.txt')



N=NMSSMTools(input_mold='inp.dat',
    package_dir=None,
    inp_file='inp.dat',
    output_file='spectr.dat',
    main_routine='run'
    )



free.Getvalue('./inp.dat',mapping={'auxiliary':'EXTPAR'})
del(free.variable_list['Kappa'].value)
free.variable_list['minusK'].value=-free.variable_list['minusK'].value
print('Start point is:')
newpoint=copy.deepcopy(free)
newpoint.print()


Data=DataFile(Dir='mcmc')

record_number=-1
try_piont=0
last_chisq=1e10

while record_number<target_number:
    try_piont+=1
    if try_piont%1000=1:
        print('Trying point %i;%i points recorded;current X2 is %.3e'%(try_piont,record_number,last_chisq))
    if try_piont>1e10:break
    try:
        spectr=N.Run(newpoint,ignore=ignore)
    except FileNotFoundError:
        newpoint=copy.deepcopy(free).Sample(step_factor=step_factor)
        continue
    if spectr.ERROR:
        newpoint=copy.deepcopy(free).Sample(step_factor=step_factor)
        continue
    
    chisq_list={}

    chisq_list['constraints']=len(spectr.constraints)*1e4
    chi_h1=chi2(spectr.MASS[25],mh)
    chi_h2=chi2(spectr.MASS[35],mh)
    if 'h1' in ism.lower():
        chisq_list['h_sm']=chi_h1
    elif 'h2' in ism.lower():
        chisq_list['h_sm']=chi_h2
    else:
        chisq_list['h_sm']=min(chi_h1,chi_h2)
    chisq_list['bsg']=chi2(spectr.LOWEN[1]*1e4,bsg)
    chisq_list['bmu']=chi2(spectr.LOWEN[4]*1e4,bmu)
    chisq_list['LHCfit']=sum(spectr.LHCFIT.values())/10.

    if hasattr(spectr,'ABUNDANCE'):
        omg=spectr.ABUNDANCE[4]
        eps=omg/0.1187
        chisq_list['DMRD']=X2(OMG=omg)
    else:
        chisq_list=1e4
    if hasattr(spectr,'NDMCROCESSSECT'):
        for name,ID,DDexp in [ ['csPsi',1,L_Psi],['csPsd,3,L_Psd'],['csNsd',4,L_Nsd]]:
            cs=abs(spectr.NDMCROCESSSECT[ID])*eps
            limit=DDexp.value(spectr.LSP[0])
            if cs>limit:
                chisq_list[name]=((cs-limit)/limit)**2/10.
            else:
                chisq_list[name]=0.
    
    Dh=numpy.loadtxt(os.payh.join(N.package_dir,'out','fine_tuning_Mh.txt'))

    Dz=spectr.FINETUNING[15]
    chisq_list['FT_Dz']=max(50-Dz,0)**2
    chisq_list['singlino']=1/spectr.NMNMIX[(1,5)]
    chisq=sum(chisq_list.values())
    if (random.random()<math.exp(max(slop_factor*min(last_chisq-chisq,0.),-745))
       ):
       last_chisq=chisq
       free=newpoint
       print(record_number,'point recorded.\n')
       print('\nnew point accepted:----------')
       print('x2=  ',chisq,'\nx2_i= ',chisq_list)
       print('Higgs masses:\n\t',[spectr.MASS.get(i) for i in [25,35,45]])
       print(f'fine tuning: Dz={spectr.FINETUNING[15]}Dh={Dh}')
       print('constraints:\n\t'+'\n\t'.join(spectr.constraints))
       print('\nParameters:')
       newpoint.Print()
       record_number+=1
       Dh_txt=f'BLOCK DELTAMH\n0\t{Dh}'
       os.system(f'echo "{Dh_txt}" >> {N.output_dir}')
       os.system(f'cat {N.output_decay_dir} >> {N.output_dir}')
       destinations={
           'input'   :os.path.join(N.record_dir,'inp.dat'+'.'+str(int(record_number))),
           'spectrun':os.path.join(N.record,'spectr.dat'+'.'+str(int(record_number)))      
        }
        shutil.copy(N.inp_dir,destinations['input'])
        shutil.copy(N.output_dir,destinations['spectrum'])
        try:
            omg_dir=os.path.join(N.record_dir,'omega.dat'+'.'+str(int(record_number)))
            shutil.copy(N.output_omega_dir,omg_dir)
        except FileExistsError:
            pass
        Data.In('Chisqure.txt').record({'Total':chisq},chisq_list)

        for file_name in [
            MINPAR','EXTPAR','NMSSMRUN','MSOFT','HMIX','MASS',
            'YD','YE','YU','TD','TE','TU','MSQ2','MSL2','MSD2','MSU2','MSE2',
            'NDMCROSSSECT','ABUNDANCE','ANNIHILATION','FINETUNING'
            ]:
            if hasattr(spectr,file_name):
                Data.In(file_name).record(getattr(spectr,file_name))
        Data.In('DeltaMH').record({0:Dh})
    else:
        print(chisq,chisq_list,'discarded')
        newpoint=copy.deepcopy(free).Sample(step_factor=step_factor)