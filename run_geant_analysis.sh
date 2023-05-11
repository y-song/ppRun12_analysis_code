#!/bin/bash

for file in `cat geant.list`
	do
	#./bin/RunppTestAna -i '/gpfs01/star/pwg/imooney/P12id/picos/pp2012embed_2022sample_2535/out/out/trigger_set_0/'$file -intype pico -c JetTree -trig ppJP2 -o 'Results/geant_526_2535/Geant'$file -N -1 -pj 5 2000 -mj 1 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 1 -geantnum 1
	./bin/RunppTestAna -i '/gpfs01/star/pwg/elayavalli/ppRun12Embpicos/'$file -intype pico -c JetTree -trig ppJP2 -o 'Results/geant/Geant'$file -N -1 -pj 10 2000 -mj 1 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -hadcorr 0.5 -towunc 0 -fakeeff 1 
	#./bin/RunppTestAna -i '/gpfs01/star/pwg/elayavalli/ppRun12Embpicos/'$file -intype pico -c JetTree -trig ppJP2 -o 'Results/geant/Geant'$file -N -1 -pj 0.2 2000 -mj 0.2 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -hadcorr 1 -towunc 0 -fakeeff 1 
done
