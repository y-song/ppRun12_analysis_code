#!/bin/bash

#source runimage.sh
source setup.sh
make
#./bin/RunppTestAna -i /gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/$1 -intype pico -c JetTree -trig ppJP2 -o /gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Results/data/out_$1 -N -1 -pj 15 2000 -mj 1 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 0 -hadcorr 1 -towunc 0 -fakeeff 1
 ./bin/RunppTestAna -i /gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/$1 -intype pico -c JetTree -trig ppJP2 -o /gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Results/data/out_$1 -N -1 -pj 15 2000 -mj 1 -pc 0.2 2000 -lja antikt -ec 1 -R 0.4 -geantnum 0 -hadcorr 1 -towunc 0 -fakeeff 1 
