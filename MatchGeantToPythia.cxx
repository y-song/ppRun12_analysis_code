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

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//! ----------------------------------------------------
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
    double q0;
    double q2;
    double rg;
    double zg;
    double mg;
    int evid;
    double weight;
    int reject;
    int mcevi;
    int mult;
    int b;
    double dr;
    RootResultStruct(TStarJetVectorJet orig, double pt, double eta, double y, double phi, double m, int n, int nch, double q0, double q2, double rg, double zg, double mg, int evid, double weight, int reject, int mcevi, int mult, int b, double dr) : orig(orig), pt(pt), eta(eta), y(y), phi(phi), m(m), n(n), nch(nch), q0(q0), q2(q2), rg(rg), zg(zg), mg(mg), evid(evid), weight(weight), reject(reject), mcevi(mcevi), mult(mult), b(b), dr(dr){};
    ClassDef(RootResultStruct, 1)
};

typedef pair<RootResultStruct, RootResultStruct> MatchedRootResultStruct;

double findweight(double val = 10.0, TH1D *hpy8 = 0, TH1D *hhw7 = 0)
{
    double val_r_py8 = hpy8->GetBinContent(hpy8->FindBin(val));
    double val_r_hw7 = hhw7->GetBinContent(hhw7->FindBin(val));
    return std::max(val_r_py8, val_r_hw7);
}

double findweight_max(double val = 10.0, TH1D *hpy8 = 0, TH1D *hhw7 = 0)
{
    double val_r_py8 = hpy8->GetBinContent(hpy8->FindBin(val));
    if (val_r_py8 > 0 && val_r_py8 < 1)
        val_r_py8 = 1. / val_r_py8;
    double val_r_hw7 = hhw7->GetBinContent(hhw7->FindBin(val));
    if (val_r_hw7 > 0 && val_r_hw7 < 1)
        val_r_hw7 = 1. / val_r_hw7;
    return std::max(val_r_py8, val_r_hw7);
}

double findweight_avg(double val = 10.0, TH1D *hpy8 = 0, TH1D *hhw7 = 0)
{
    double val_r_py8 = hpy8->GetBinContent(hpy8->FindBin(val));
    if (val_r_py8 > 0 && val_r_py8 < 1)
        val_r_py8 = 1. / val_r_py8;
    double val_r_hw7 = hhw7->GetBinContent(hhw7->FindBin(val));
    if (val_r_hw7 > 0 && val_r_hw7 < 1)
        val_r_hw7 = 1. / val_r_hw7;
    return (val_r_py8 + val_r_hw7) / 2;
}

double findweight_v2(double val = 10.0, TH1D *hweight = 0)
{
    return hweight->GetBinContent(hweight->FindBin(val));
}

int MatchGeantToPythia(string McFile, string PpFile, string OutFile = "test.root", string badrunlist="mcevi.list", int RADIUS = 4, int mode = 0)
{
    bool RejectHiweights = true;
    float RCut = (float)RADIUS / 10;
    float EtaCut = 1.0 - RCut;
    int nj = 0;

    TString dir = "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Results/";

    TFile *Mcf = new TFile(dir + "pythia/" + McFile);
    TTree *McChain = (TTree *)Mcf->Get("ResultTree");
    McChain->BuildIndex("runid", "eventid");

    TClonesArray *McJets = new TClonesArray("TStarJetVectorJet");
    McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
    McChain->SetBranchAddress("Jets", &McJets);

    int mcrunid;
    McChain->SetBranchAddress("runid", &mcrunid);
    int mceventid;
    McChain->SetBranchAddress("eventid", &mceventid);
    double mcweight;
    McChain->SetBranchAddress("weight", &mcweight);
    int mcnjets = 0;
    McChain->SetBranchAddress("njets", &mcnjets);
    double mcpt[1000];
    McChain->SetBranchAddress("pt", mcpt);
    double mcm[1000];
    McChain->SetBranchAddress("m", mcm);
    double mcq0[1000];
    McChain->SetBranchAddress("q0", mcq0);
    double mcq2[1000];
    McChain->SetBranchAddress("q2", mcq2);
    int mcb[1000];
    McChain->SetBranchAddress("b", mcb);
    double mcdr[1000];
    McChain->SetBranchAddress("dr", mcdr);
    double mcrg[1000];
    McChain->SetBranchAddress("rg", mcrg);
    double mczg[1000];
    McChain->SetBranchAddress("zg", mczg);
    double mcmg[1000];
    McChain->SetBranchAddress("mg", mcmg);
    int mcn[1000];
    McChain->SetBranchAddress("n", mcn);
    int mcnch[1000];
    McChain->SetBranchAddress("nch", mcnch);
    int mcreject[1000];
    McChain->SetBranchAddress("reject", mcreject);
    int mcmult;   
    McChain->SetBranchAddress("mult", &mcmult);

    TFile *Ppf = new TFile(dir + "geant/" + PpFile);
    TTree *PpChain = (TTree *)Ppf->Get("ResultTree");
    PpChain->BuildIndex("runid", "eventid");

    TClonesArray *PpJets = new TClonesArray("TStarJetVectorJet");
    PpChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
    PpChain->SetBranchAddress("Jets", &PpJets);

    int ppeventid;
    PpChain->SetBranchAddress("eventid", &ppeventid);
    int pprunid;
    PpChain->SetBranchAddress("runid", &pprunid);
    double ppweight;
    PpChain->SetBranchAddress("weight", &ppweight);
    int ppnjets = 0;
    PpChain->SetBranchAddress("njets", &ppnjets);
    double rcpt[1000];
    PpChain->SetBranchAddress("pt", rcpt);
    double rcm[1000];
    PpChain->SetBranchAddress("m", rcm);
    double rcq0[1000];
    PpChain->SetBranchAddress("q0", rcq0);
    double rcq2[1000];
    PpChain->SetBranchAddress("q2", rcq2);
    int rcb[1000];
    PpChain->SetBranchAddress("b", rcb);
    double rcdr[1000];
    PpChain->SetBranchAddress("dr", rcdr);
    double rcrg[1000];
    PpChain->SetBranchAddress("rg", rcrg);
    double rczg[1000];
    PpChain->SetBranchAddress("zg", rczg);
    double rcmg[1000];
    PpChain->SetBranchAddress("mg", rcmg);
    int rcn[1000];
    PpChain->SetBranchAddress("n", rcn);
    int rcnch[1000];
    PpChain->SetBranchAddress("nch", rcnch);
    int rcreject[1000];
    PpChain->SetBranchAddress("reject", rcreject);
    int rcmult;
    PpChain->SetBranchAddress("mult", &rcmult);

    //! Output and histograms
    TFile *fout = new TFile(dir + "matched/" + OutFile, "RECREATE");
    TTree *MatchedTree = new TTree("MatchedTree", "Matched Jets");
    double pt;
    MatchedTree->Branch("pt", &pt, "pt/D");
    double pts;
    MatchedTree->Branch("pts", &pts, "pts/D");
    double m;
    MatchedTree->Branch("m", &m, "m/D");
    double ms;
    MatchedTree->Branch("ms", &ms, "ms/D");
    double q0;
    MatchedTree->Branch("q0", &q0, "q0/D");
    double q0s;
    MatchedTree->Branch("q0s", &q0s, "q0s/D");
    double q2;
    MatchedTree->Branch("q2", &q2, "q2/D");
    double q2s;
    MatchedTree->Branch("q2s", &q2s, "q2s/D");
    int b;
    MatchedTree->Branch("b", &b, "b/I");
    int bs;
    MatchedTree->Branch("bs", &bs, "bs/I");
    double dr;
    MatchedTree->Branch("dr", &dr, "dr/D");
    double drs;
    MatchedTree->Branch("drs", &drs, "drs/D");
    double rg;
    MatchedTree->Branch("rg", &rg, "rg/D");
    double rgs;
    MatchedTree->Branch("rgs", &rgs, "rgs/D");
    double zg;
    MatchedTree->Branch("zg", &zg, "zg/D");
    double zgs;
    MatchedTree->Branch("zgs", &zgs, "zgs/D");
    double mg;
    MatchedTree->Branch("mg", &mg, "mg/D");
    double mgs;
    MatchedTree->Branch("mgs", &mgs, "mgs/D");
    int n;
    MatchedTree->Branch("n", &n, "n/I");
    int ns;
    MatchedTree->Branch("ns", &ns, "ns/I");
    int nch;
    MatchedTree->Branch("nch", &nch, "nch/I");
    int nchs;
    MatchedTree->Branch("nchs", &nchs, "nchs/I");
    double w;
    MatchedTree->Branch("w", &w, "w/D");
    double ws;
    MatchedTree->Branch("ws", &ws, "ws/D");
    int evid;
    MatchedTree->Branch("evid", &evid, "evid/I");
    int evids;
    MatchedTree->Branch("evids", &evids, "evids/I");
    int reject;
    MatchedTree->Branch("reject", &reject, "reject/I");
    int rejects;
    MatchedTree->Branch("rejects", &rejects, "rejects/I");
    int mcevi;
    MatchedTree->Branch("mcevi", &mcevi, "mcevi/I");
    int mult;
    MatchedTree->Branch("mult", &mult, "mult/I");
    int mults;
    MatchedTree->Branch("mults", &mults, "mults/I");

    //! Loop over particle level
    int missed = 0;
    int N = McChain->GetEntries();
    cout << "Number of Pythia events: " << N << endl;
    cout << "Number of Geant events:  " << PpChain->GetEntries() << endl;
    McChain->GetEntry(4698); //this event contains a truth jet
    TStarJetVectorJet *dummyjet = (TStarJetVectorJet *)McJets->At(0);
    
    // read lists of pythia events to skip
    vector<int> badmcevi_list;
    FILE *file = fopen("./mcevi.list", "r");
    char file_content[256];
    badmcevi_list.push_back(0);
    while (fgets(file_content, 256, file) != NULL)
    {
	file_content[strcspn(file_content, "\n")] = 0;
	badmcevi_list.push_back(atoi(file_content));
    }
    
    vector<int> badmcevi_list2;
    FILE *file2 = fopen("./mcevi2.list", "r");
    char file_content2[256];
    badmcevi_list2.push_back(0);
    while (fgets(file_content2, 256, file2) != NULL)
    {
	file_content2[strcspn(file_content2, "\n")] = 0;
	badmcevi_list2.push_back(atoi(file_content2));
    }

    for (Long64_t mcEvi = 0; mcEvi < N; ++mcEvi) //event loop
    {
        if (!(mcEvi % 500000))
            cout << "Working on " << mcEvi << " / " << N << endl;
        
        if (binary_search(badmcevi_list.begin(), badmcevi_list.end(), mcEvi))
		continue;

	if (binary_search(badmcevi_list2.begin(), badmcevi_list2.end(), mcEvi))
		continue;

	McChain->GetEntry(mcEvi);
       
	int ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
	//! Fill results in vectors for easier manipulation
        //! Also check whether there's something true in the acceptance
        //bool TruthInAcceptance = false;
        vector<RootResultStruct> mcresult;
        for (int j = 0; j < mcnjets; ++j)
        {
 	    TStarJetVectorJet *mcjet = (TStarJetVectorJet *)McJets->At(j);

	    // THE FOLLOWING IS NOT THE DEFAULT
	    if (fabs(mcjet->Eta()) < EtaCut)  
	    {
		mcresult.push_back(RootResultStruct(*mcjet, mcjet->Pt(), mcjet->Eta(), mcjet->Rapidity(), mcjet->Phi(), mcjet->M(), mcn[j], mcnch[j], mcq0[j], mcq2[j], mcb[j], mcdr[j], mcrg[j], mczg[j], mcmg[j], mceventid, mcweight, mcreject[j], mcEvi, mcmult));
		nj++;
	    }
        }

        //! Get geant
        PpChain->GetEntry(ppevi);
        vector<RootResultStruct> ppresult;
        for (int j = 0; j < ppnjets; ++j)
        {
            TStarJetVectorJet *ppjet = (TStarJetVectorJet *)PpJets->At(j);

           if (rcreject[j]>0)
	    {
		cout << "mcevi,ppevi,rerejct[j]=" << mcEvi << ", " << ppevi << ", " << rcreject[j] << endl;
	    }
	    if (fabs(ppjet->Eta()) < EtaCut && ppevi>=0)
            {
                ppresult.push_back(RootResultStruct(*ppjet, ppjet->Pt(), ppjet->Eta(), ppjet->Rapidity(), ppjet->Phi(), ppjet->M(), rcn[j], rcnch[j], rcq0[j], rcq2[j], rcb[j], rcdr[j], rcrg[j], rczg[j], rcmg[j], ppeventid, ppweight, rcreject[j], mcEvi, rcmult));
            }
        }
 
        //! Sort them together
        vector<MatchedRootResultStruct> MatchedResult;
	if (mcresult.size() > 0)
        {
	for (vector<RootResultStruct>::iterator mcit = mcresult.begin(); mcit != mcresult.end();)
        {
            bool matched = false;
            if (ppresult.size() > 0)
	    {
	    	for (vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();)
           	{
               		 Double_t deta = mcit->y - ppit->y;
               		 Double_t dphi = TVector2::Phi_mpi_pi(mcit->phi - ppit->phi);
               		 Double_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
               		 if (dr < RCut)
               		 {
                   		 MatchedResult.push_back(MatchedRootResultStruct(*mcit, *ppit));
                   		 ppit = ppresult.erase(ppit);
                 		 matched = true;
                   		 break;
               		 }
                	 else
               		 {
                 		 ++ppit;
               		 }
           	 }
	    }
            if (matched)
            {
                mcit = mcresult.erase(mcit);
            }
            else
            {
        	MatchedResult.push_back(MatchedRootResultStruct(*mcit, RootResultStruct(*dummyjet, -9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9, ppeventid, ppweight, -9, mcEvi, rcmult)));
                ++mcit;
            }
        }
	}
	for (vector<RootResultStruct>::iterator ppit = ppresult.begin(); ppit != ppresult.end();)
	{
		MatchedResult.push_back(MatchedRootResultStruct(RootResultStruct(*dummyjet, -9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9, mceventid, mcweight, -9, mcEvi, mcmult), *ppit));
		++ppit;
	}
	for (vector<MatchedRootResultStruct>::iterator res = MatchedResult.begin(); res != MatchedResult.end(); ++res)
        {
            pt = res->first.pt;
            pts = res->second.pt;
            m = res->first.m;
            ms = res->second.m;
            n = res->first.n;
            ns = res->second.n;
            nch = res->first.nch;
            nchs = res->second.nch;
            q0 = res->first.q0;
            q0s = res->second.q0;
            q2 = res->first.q2;
            q2s = res->second.q2;
            b = res->first.b;
	    bs = res->second.b;
	    dr = res->first.dr;
	    drs = res->second.dr;
	    rg = res->first.rg;
            rgs = res->second.rg;
            zg = res->first.zg;
            zgs = res->second.zg;
            mg = res->first.mg;
            mgs = res->second.mg;
            evid = res->first.evid;
            evids = res->second.evid;
            w = res->first.weight;
            ws = res->second.weight;
	    reject = res->first.reject;
	    rejects = res->second.reject;
	    mult = res->first.mult;
	    mults = res->second.mult;
	    if (res->second.mcevi > res->first.mcevi)
	    {
		mcevi = res->second.mcevi;
	    }
	    else
	    {
		mcevi = res->first.mcevi;
	    }
            MatchedTree->Fill();
        }
    }
    fout->Write();
    return missed;
}
