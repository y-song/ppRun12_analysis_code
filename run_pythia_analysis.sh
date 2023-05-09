#!/bin/bash

for file in `cat geant.list`
	do
	./bin/RunppTestAna -i '/gpfs01/star/pwg/elayavalli/ppRun12Embpicos/'$file -intype mcpico -c JetTreeMc -trig all -o 'Results/pythia/Pythia'$file -N -1 -pj 5 2000 -mj 0.0 -pc 2.0 30 -lja antikt -ec 1 -R 0.4 -geantnum 1 -jetnef 1
	#./bin/RunppTestAna -i '/gpfs01/star/pwg/elayavalli/ppRun12Embpicos/'$file -intype mcpico -c JetTreeMc -trig all -o 'Results/pythia/Pythia'$file -N -1 -pj 0.2 2000 -mj 0.0 -pc 0.2 1000 -lja antikt -ec 1 -R 0.4 -geantnum 1 -jetnef 1
	
done
