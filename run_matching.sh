#!/bin/bash

#root -l -q 'Test.cxx+("geant0901_allevt.root", "geant0901.root")'
#root -l -q 'GetMcevi.cxx+("pythia0901.root", "geant0901_allevt.root")'
root -l -q 'MatchGeantToPythia.cxx+("pythia0901.root", "geant0901_eff096.root", "all0901_eff096.root", "mcevi.list")'
#root -l -q 'MatchGeantToPythia.cxx+("pythia0901.root", "geant0901_trkcut_evtptcut_jp2.root", "all0901_trkcut_evtptcut_jp2.root", "mcevi.list")'
