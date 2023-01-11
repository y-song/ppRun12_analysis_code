//! Code to read in geant + pythia output trees and match them
//! Raghav Kunnawalkam Elayavalli & Kolja Kauder
//! contact - raghavke@wayne.edu
//! HAS to be compiled,
//! root -l macros/MatchGeantToPythia.cxx+

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLine.h>

#include <TROOT.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

#include <TStarJetVector.h>
#include <TStarJetVectorJet.h>
#include <TStarJetPicoReader.h>

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

using namespace std;

//! Load helper macro
#include "NewGeantWeightReject.hh"

class RootResultStruct
{
public:
    TStarJetVectorJet orig;
    double pt;
    double eta;
    double y;
    double phi;
    double m;
    int n;
    int nch;
    double q;
    double rg;
    double zg;
    double mg;
    int evid;
    double weight;

    RootResultStruct(TStarJetVectorJet orig, double pt, double eta, double y, double phi, double m, int n, int nch, double q, double rg, double zg, double mg, int evid, double weight) : orig(orig), pt(pt), eta(eta), y(y), phi(phi), m(m), n(n), nch(nch), q(q), rg(rg), zg(zg), mg(mg), evid(evid), weight(weight){};
    ClassDef(RootResultStruct, 1)
};

//! ----------------------------------------------------
int GetMcevi(string McFile, string PpFile)
{
    TString dir = "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Results/";

    TFile *Mcf = new TFile(dir + "pythia/" + McFile);
    TTree *McChain = (TTree *)Mcf->Get("ResultTree");
    McChain->BuildIndex("runid", "eventid");

    int mcrunid;
    McChain->SetBranchAddress("runid", &mcrunid);
    int mceventid;
    McChain->SetBranchAddress("eventid", &mceventid);
    
    TFile *Ppf = new TFile(dir + "geant/" + PpFile);
    TTree *PpChain = (TTree *)Ppf->Get("ResultTree");
    PpChain->BuildIndex("runid", "eventid");

    int ppeventid;
    PpChain->SetBranchAddress("eventid", &ppeventid);
    int pprunid1;
    PpChain->SetBranchAddress("runid1", &pprunid1);
    int pprunid;
    PpChain->SetBranchAddress("runid", &pprunid);
    double ppvz;
    PpChain->SetBranchAddress("vz", &ppvz);
    int pphighw;
    PpChain->SetBranchAddress("highw", &pphighw);

    cout << "Number of Pythia events: " << McChain->GetEntries() << endl;
    cout << "Number of Geant events:  " << PpChain->GetEntries() << endl;

    // read bad run list
    vector<int> runid1_list;
    FILE *file = fopen("pp200Y12_badrun_isaac.list", "r");
    char file_content[256];
    runid1_list.push_back(0);
    while (fgets(file_content, 256, file) != NULL)
    {
	file_content[strcspn(file_content, "\n")] = 0;
	runid1_list.push_back(atoi(file_content));
    }

    // find pairs of (ppeventid,pprunid) that are in the bad run list / have high vz / have high weight
    vector<int> ppeventid_list;
    vector<int> pprunid_list;
    for (Long64_t ppevi = 0; ppevi < PpChain->GetEntries(); ++ppevi)
    {
	PpChain->GetEntry(ppevi);
	// check if pprunid1 is in bad runid1 list
	if (find(runid1_list.begin(), runid1_list.end(), pprunid1) != runid1_list.end() || abs(ppvz)>30 || pphighw>0)
	{
	    ppeventid_list.push_back(ppeventid);
	    pprunid_list.push_back(pprunid);
        }
    }
    cout << "Number of bad Geant events:  " << ppeventid_list.size() << endl;

    // find mcevi that corresponds to (ppeventid,pprunid) that are in the bad run list / have high vz / have high weight
    vector<int> mcevi_list;
    for (unsigned int i = 0; i < ppeventid_list.size(); i++)
    {
	int mcevi = McChain->GetEntryNumberWithIndex(pprunid_list.at(i), ppeventid_list.at(i));
	if (mcevi >= 0)
	{
		mcevi_list.push_back(mcevi);
	}
    }
    cout << "Number of Pythia events matched to bad Geant events: " << mcevi_list.size() << endl;

    // find mcevi that do not have a matching ppevi
    vector<int> mcevi_list2;
    for (Long64_t mcevi = 0; mcevi < McChain->GetEntries(); ++mcevi)
    {
	McChain->GetEntry(mcevi);
	int ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
	if (ppevi < 0) // event not embedded
	{
	    mcevi_list2.push_back(mcevi);
	}
    }
    cout << "Number of Pythia events that are not embedded: " << mcevi_list2.size() << endl;

    ofstream out_file("./mcevi.list");
    ostream_iterator<int> output_iterator(out_file, "\n");
    copy(mcevi_list.begin(), mcevi_list.end(), output_iterator);

    ofstream out_file2("./mcevi2.list");
    ostream_iterator<int> output_iterator2(out_file2, "\n");
    copy(mcevi_list2.begin(), mcevi_list2.end(), output_iterator2);

    return 0;
}
