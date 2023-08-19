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
    int pid1;
    int pid2;
    double pt1;
    double pt2;
    double z;
    double epair;
    double qlead;
    double qsublead;
    double y1;
    double y2;
    double phi1;
    double phi2;
    double nef;
    RootResultStruct(TStarJetVectorJet orig, double pt, double eta, double y, double phi, double m, int n, int nch, double q0, double q2, double rg, double zg, double mg, int b, double dr, double pt1, double pt2, int pid1, int pid2, double z, double epair, double qlead, double qsublead, double y1, double y2, double phi1, double phi2, double nef, int evid, double weight, int reject, int mcevi, int mult) : orig(orig), pt(pt), eta(eta), y(y), phi(phi), m(m), n(n), nch(nch), q0(q0), q2(q2), rg(rg), zg(zg), mg(mg), b(b), dr(dr), pid1(pid1), pid2(pid2), pt1(pt1), pt2(pt2), z(z), epair(epair), qlead(qlead), qsublead(qsublead), y1(y1), y2(y2), phi1(phi1), phi2(phi2), nef(nef), evid(evid), weight(weight), reject(reject), mcevi(mcevi), mult(mult){};
    ClassDef(RootResultStruct, 1)
};

int PlotVz()
{
    TString dir = "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Results/";

    TFile *Mcf = new TFile(dir + "pythia/" + "pythia0817.root");
    TTree *McChain = (TTree *)Mcf->Get("ResultTree");
    McChain->BuildIndex("runid", "eventid");

    TClonesArray *McJets = new TClonesArray("TStarJetVectorJet");
    McChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
    McChain->SetBranchAddress("Jets", &McJets);

    int mcrunid;
    McChain->SetBranchAddress("runid", &mcrunid);
    int mceventid;
    McChain->SetBranchAddress("eventid", &mceventid);
    double mcvz = 0;
    McChain->SetBranchAddress("vz", &mcvz);

    TFile *Ppf = new TFile(dir + "geant/" + "geant0817.root");
    TTree *PpChain = (TTree *)Ppf->Get("ResultTree");
    PpChain->BuildIndex("runid", "eventid");

    TClonesArray *PpJets = new TClonesArray("TStarJetVectorJet");
    PpChain->GetBranch("Jets")->SetAutoDelete(kFALSE);
    PpChain->SetBranchAddress("Jets", &PpJets);

    int ppeventid;
    PpChain->SetBranchAddress("eventid", &ppeventid);
    int pprunid;
    PpChain->SetBranchAddress("runid", &pprunid);
    double ppvz = 0;
    PpChain->SetBranchAddress("vz", &ppvz);
    TH1D *h_vz = new TH1D("h", ";vz", 120, -30, 30);
    gStyle->SetOptStat(0);

    //! Loop over particle level
    int missed = 0;
    int N = McChain->GetEntries();
    cout << "Number of Pythia events: " << N << endl;
    cout << "Number of Geant events:  " << PpChain->GetEntries() << endl;
 
    // read lists of pythia events to skip
    vector<int> badmcevi_list;
    FILE *file = fopen("./mcevi_new.list", "r");
    char file_content[256];
    badmcevi_list.push_back(0);

    while (fgets(file_content, 256, file) != NULL)
    {
        file_content[strcspn(file_content, "\n")] = 0;
        badmcevi_list.push_back(atoi(file_content));
    }

    vector<int> badmcevi_list2;
    FILE *file2 = fopen("./mcevi2_new.list", "r");
    char file_content2[256];
    badmcevi_list2.push_back(0);
    while (fgets(file_content2, 256, file2) != NULL)
    {
        file_content2[strcspn(file_content2, "\n")] = 0;
        badmcevi_list2.push_back(atoi(file_content2));
    }

    vector<int> badmcevi_list3;
    FILE *file3 = fopen("./mcevi3_new.list", "r");
    char file_content3[256];
    badmcevi_list3.push_back(0);
    while (fgets(file_content3, 256, file3) != NULL)
    {
        file_content3[strcspn(file_content3, "\n")] = 0;
        badmcevi_list3.push_back(atoi(file_content3));
    }

    vector<int> badmcevi_list4;
    FILE *file4 = fopen("./mcevi4_new.list", "r");
    char file_content4[256];
    badmcevi_list4.push_back(0);
    while (fgets(file_content4, 256, file4) != NULL)
    {
        file_content4[strcspn(file_content4, "\n")] = 0;
        badmcevi_list4.push_back(atoi(file_content4));
    }

    for (Long64_t mcEvi = 0; mcEvi < N; ++mcEvi) // event loop
    {
        if (!(mcEvi % 500000))
            cout << "Working on " << mcEvi << " / " << N << endl;

        if (binary_search(badmcevi_list.begin(), badmcevi_list.end(), mcEvi))
            continue;

        if (binary_search(badmcevi_list2.begin(), badmcevi_list2.end(), mcEvi))
            continue;

        if (binary_search(badmcevi_list3.begin(), badmcevi_list3.end(), mcEvi))
            continue;

        if (binary_search(badmcevi_list4.begin(), badmcevi_list4.end(), mcEvi))
            continue;

        McChain->GetEntry(mcEvi);
        vector<RootResultStruct> mcresult;

        //if (abs(mcvz) > 15)
        //     continue;

        int ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
        PpChain->GetEntry(ppevi);
        vector<RootResultStruct> ppresult;

        if (ppevi >= 0)
        {
            h_vz->Fill(ppvz);
        }
    }

    TFile *Mcf2 = new TFile(dir + "pythia/" + "pythia0901.root");
    TTree *McChain2 = (TTree *)Mcf2->Get("ResultTree");
    McChain2->BuildIndex("runid", "eventid");

    mcrunid = 0;
    McChain2->SetBranchAddress("runid", &mcrunid);
    mceventid = 0;
    McChain2->SetBranchAddress("eventid", &mceventid);
    mcvz = 0;
    McChain2->SetBranchAddress("vz", &mcvz);

    TFile *Ppf2 = new TFile(dir + "geant/" + "geant0901.root");
    TTree *PpChain2 = (TTree *)Ppf2->Get("ResultTree");
    PpChain2->BuildIndex("runid", "eventid");

    ppeventid = 0;
    PpChain2->SetBranchAddress("eventid", &ppeventid);
    pprunid = 0;
    PpChain2->SetBranchAddress("runid", &pprunid);
    ppvz = 0;
    PpChain2->SetBranchAddress("vz", &ppvz);

    TH1D *h_vz_2 = new TH1D("h", ";vz", 120, -30, 30);

    //! Loop over particle level
    N = McChain2->GetEntries();
    cout << "Number of Pythia events: " << N << endl;
    cout << "Number of Geant events:  " << PpChain2->GetEntries() << endl;

    // read lists of pythia events to skip
    vector<int> badmcevi_list22;
    FILE *file22 = fopen("./mcevi.list", "r");
    char file_content22[256];
    badmcevi_list22.push_back(0);
    while (fgets(file_content22, 256, file22) != NULL)
    {
        file_content22[strcspn(file_content22, "\n")] = 0;
        badmcevi_list22.push_back(atoi(file_content22));
    }

    vector<int> badmcevi_list222;
    FILE *file222 = fopen("./mcevi2.list", "r");
    char file_content222[256];
    badmcevi_list222.push_back(0);
    while (fgets(file_content222, 256, file222) != NULL)
    {
        file_content222[strcspn(file_content222, "\n")] = 0;
        badmcevi_list222.push_back(atoi(file_content222));
    }

    for (Long64_t mcEvi = 0; mcEvi < N; ++mcEvi) // event loop
    {
        if (!(mcEvi % 500000))
            cout << "Working on " << mcEvi << " / " << N << endl;

        if (binary_search(badmcevi_list22.begin(), badmcevi_list22.end(), mcEvi))
            continue;

        if (binary_search(badmcevi_list222.begin(), badmcevi_list222.end(), mcEvi))
            continue;

        McChain2->GetEntry(mcEvi);

        int ppevi = PpChain2->GetEntryNumberWithIndex(mcrunid, mceventid);
        PpChain2->GetEntry(ppevi);

        if (ppevi >= 0)
        {
            h_vz_2->Fill(ppvz);
        }
    }

TCanvas *c = new TCanvas("c", "c", 1000, 600);
h_vz_2->SetMinimum(0);
h_vz_2->SetLineColor(2);
h_vz_2->DrawNormalized();
h_vz->DrawNormalized("same");
c->SaveAs("plots/vz.png");
return missed;
}