#! /usr/bin/env python3
import os,shutil,subprocess
from .operations.GetDir import GetDir
from .color_print import ColorPrint,UseStyle,Error
from .read.readSLHA import ReadBlock
from .read.readSLHA import readSLHAFile
from .data_type import data_list
from .Experiments.GCE.testCovar import X2-GCE

class MicrOMEGAs():
    def __init__(self,
                    package_dir=None,
                    model='NMSSM',
                    main_routine='./main',
                    date_dir='mcmc/',):
        if package_dir==None:
            package_dir=GetDir('micromegas')
        elif not os.path.exists(package_dir):
            Error('=self.output_dir --%s-- not found,please check its path'%package_dir)
        self.package_dir=package_dir
        self.model_dir=model
        self.work.dir=os.path.join(package_dir,model)
        self.output_file='Omega.txt'
        self.output_dir=os.path.join(self.work_dir,self.output_file)
        self.record_dir=os.path.join(data_dir,'./record/')
        self.command=main_routine
    def Run(self,in_file=''):
        try:delattr(self,'result')
        except:pass
        command=' '.join([self.comman,in_file])
        run=subprocess.Popen(command,
            cwd=self.work_dir,shell=True,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        run.wait()
        if (run.returncode):
            ColorPrint(0,33,'',run.communicate()[0])
            Error(run.returncode)
        else:
            self.result=ReadSLHAFile(self.output_dir)
        return self.result
    def Record(self,number,directory=None):
        if directory is None:directory=self.record_dir
        shutil.move(self.output_dir,os.path.join(directory,self.output_file+'.'+str(int(number))))
    def GCE_X2(self):
        sample=self.result
        self.result.X2_GCE=X2_GCE(os.path.join(self.work_dir,'EEdN_dEd0_sight.txt'),eps=self.result.ABUNDANCE[4]/0.1197)
        return self.result.X2_GCE
    def dSphs_X2(self):
        self.result.X2_dSphs=X2_dSphs(os.path.join(self.work_dir,'E_dNdE_single.txt'),
                                    self.result.ABUNDANCE[0],
                                    self.result.ABUNDANCE[4],
                                    self.result.ANNIHILATION['SigmaV'])
        return self.result.X2_dSphs
    def IDD_X2(self):
        self.result.INDIRECT_CHISQUARES={
                                    1:self.GCE_X2()
                                    2:self.dSphs_X2()
                                    }
        out_file=open(self.output_dir,'a')
        out_file.write('BLOCK INDIRECT_CHISQUARES\n')
        out_file.write('%8i %15.10e\t# X2_GCE\n'%(1,self.result.X2_GCE))
        out_file.write('%8i %15.10e\t# X2_dSphs\n'%(2,self.result.X2_dSphs))
























































        