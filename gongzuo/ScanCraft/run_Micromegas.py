#! /usr/bin/env python3
import sys,os,re,copy,shutil
# sys.path.append('/home/heyangle/Desktop/ScanCommando/ScanCommando')

from command.MicrOMEGAs import MicrOMEGAs
from command.file_operations.GetSamples import GetSamples

spectr_list=GetSamples('./mcmc/record/',patterns='spectr')
MOmega=MicrOMEGAs()

number=0
#print(os.getcwd())
#print('***:',spectr_list)
for spectr in spectr_list:
    sample=MOmega.Run('spectr.dat')
    print(sample.__dict__)

# print(sample.ABUNDANCE[4],MOmega.IDD_X2(sample),MOmega.GCE_X2())

# sample=MOmega.Run('./spectr.20928')
# print(sample.ABUNDANCE)
# print(MOmega.GCE_X2())
# print(MOmega.IDD_X2(sample))