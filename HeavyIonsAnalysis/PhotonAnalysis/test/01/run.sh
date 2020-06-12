#!/bin/bash                                          
                                                     
function command(){                                  
     echo ">>>" $@                                   
     eval $@                                         
}

workDir=$PWD  
                                                                       
logDir=${workDir}"/testBATCHJOBS/"
command "mkdir -p $logDir"

castorDirOut="/eos/cms/store/group/phys_diffraction/diphoton/aug_reco_data_check_for_lumi/hiforest/01/"
#command "mkdir -p $castorDirOut"


j=0
jobNb=""

for fileList in `ls filelist_??`
do
    echo $fileList;
    jobNb=${j};
    let j=${j}+1;
    name="mycheck_hiforest_${jobNb}"
    outfilename="${name}.root"


#Start to write the script
    cat > job_Onia2MuMu_${jobNb}.sh << EOF
#!/bin/bash

function command(){                                  
     echo ">>>" \$@                                   
     eval \$@                                         
}


cd ${CMSSW_BASE}/src
eval \`scramv1 runtime -sh\`
cd -
cp ${workDir}/runForestAOD_mod_pbpb_DATA.py tmp.py
cp ${workDir}/${fileList} .
echo "before running cmsRun"
ls -l

command "cmsRun tmp.py inputFiles_clear inputFiles_load=${fileList} print 2>&1 | tee -a $logDir${name}.log"

echo "after running cmsRun"
ls -l

command "cp hiforest.root ${castorDirOut}${outfilename}"
#command "cp onia2MuMuPAT.root ${castorDirOut}${outfilename}"
command "gzip -f $logDir${name}.log"
#command "rm -f tmp.py"
#command "rm -f ${fileList}"

EOF

	chmod 755 job_Onia2MuMu_${jobNb}.sh
	
	#bsub -q cmscaf1nd -J $name ${workDir}/job_Onia2MuMu_EOS_${jobNb}.sh
	bsub  -q 2nd -J $name ${workDir}/job_Onia2MuMu_${jobNb}.sh

done
