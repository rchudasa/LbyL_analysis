# coding: utf-8
# A simple hack to use python to submit jobs to bsub. It generates .csh scripts that are then sent to the scheduler.
# It creates file lists for all the samples, as well as output directories on eos
# Author: Stepan Obraztsov

import os, re
import commands
import math, time
import sys
import datetime
import subprocess

print 
print 'START'
print 

#Sample list
#files = ['_01','_02','_03','_04','_05']
files = ['_01','_02','_03','_04','_05','_06','_07','_08','_09','_10','_11','_12','_13','_14','_15','_16','_17','_18','_19','_20','_21','_22','_23','_24','_25','_26','_27','_28','_29','_30','_31','_32','_33','_34','_35','_36','_37','_38','_39','_40','_41','_42','_43','_44','_45','_46','_47','_48','_49','_50','_51','_52','_53','_54','_55','_56','_57','_58','_59','_60','_61','_62','_63','_64','_65','_66','_67','_68','_69','_70','_71','_72','_73','_74','_75','_76','_77','_78','_79','_80','_81','_82','_83','_84','_85','_86','_87','_88','_89','_90','_91','_92','_93','_94','_95','_96','_97','_98','_99','_100','_101','_102','_103','_104','_105','_106','_107','_108','_109','_110','_111','_112','_113','_114','_115','_116','_117','_118','_119','_120','_121','_122','_123','_124','_125','_126','_127','_128','_129','_130','_131','_132','_133','_134','_135','_136','_137','_138','_139','_140','_141','_142','_143','_144','_145','_146','_147','_148','_149','_150','_151','_152','_153','_154','_155','_156','_157','_158','_159','_160','_161','_162','_163','_164','_165','_166','_167','_168','_169','_170','_171','_172','_173','_174','_175','_176','_177','_178','_179','_180','_181','_182','_183','_184','_185','_186','_187','_188','_189','_190','_191','_192','_193','_194','_195']
#,'_196','_197','_198','_199']


queue = "2nd" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
#interval = 4 # number files to be processed in a single job, take care to split your file so that you run on all files.


#Everything for data samples here

macro = {}


#Obsolete. Single macro for all type of samples. 
#macro = 'data_ele_gsf_tree'
macro = 'check_data_driven_reco_eff'

# generate .csh scripts
for x in files:
    with open('closure'+x+'.csh', 'w') as fout:
        fout.write("#!/bin/tcsh\n")
        fout.write("setenv HOME /afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch4/src/HeavyIonsAnalysis/PhotonAnalysis/test/aug_reco_data_check_for_lumi/hiforest/reco_eff/ged\n")
        #fout.write("setenv HOME $PWD\n")
        fout.write("setenv WORK $PWD\n")
        fout.write("cd ${HOME}\n")
        fout.write("cmsenv\n")
        fout.write("cp "+macro+".C ${WORK}\n")
        fout.write("cp filelist"+x+".txt ${WORK}\n")
        fout.write("cd ${WORK}\n")
        fout.write("g++ "+macro+".C `root-config --cflags --libs` -O2 -o "+macro+".exe\n")
        fout.write("./"+macro+".exe " "\"filelist"+x+".txt\" ""\n")
        fout.write("cp test.root /afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch4/src/HeavyIonsAnalysis/PhotonAnalysis/test/aug_reco_data_check_for_lumi/hiforest/reco_eff/ged"+"/ntuple"+x+".root\n")
        #fout.write("cp test.root $PWD"+"/ntuple"+x+".root\n")
        fout.write("cp "+macro+".C /afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch4/src/HeavyIonsAnalysis/PhotonAnalysis/test/aug_reco_data_check_for_lumi/hiforest/reco_eff/ged"+"/"+macro+".C\n")      
        #fout.write("cp "+macro+".C $PWD"+"/"+macro+".C\n")      
    os.system('chmod 755 closure'+x+'.csh')
    os.system('bsub -q '+queue+' closure'+x+'.csh ')
    print "job "+x+" "+" submitted"
   
print
#print "your jobs:"
#os.system("bjobs")
print
print 'END'
print



