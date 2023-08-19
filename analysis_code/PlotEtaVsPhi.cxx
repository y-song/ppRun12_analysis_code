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

int PlotEtaVsPhi(string McFile, string PpFile, string OutFile = "test.root")
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
    /*double rcpt[1000];
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
    double rcpt1[1000];
    PpChain->SetBranchAddress("pt1", rcpt1);
    double rcpt2[1000];
    PpChain->SetBranchAddress("pt2", rcpt2);
    double rcqlead[1000];
    PpChain->SetBranchAddress("qlead", rcqlead);
    double rcqsublead[1000];
    PpChain->SetBranchAddress("qsublead", rcqsublead);
    double rcy1[1000];
    PpChain->SetBranchAddress("y1", rcy1);
    double rcy2[1000];
    PpChain->SetBranchAddress("y2", rcy2);
    double rcphi1[1000];
    PpChain->SetBranchAddress("phi1", rcphi1);
    double rcphi2[1000];
    PpChain->SetBranchAddress("phi2", rcphi2);
    int rcpid1[1000];
    PpChain->SetBranchAddress("pid1", rcpid1);
    int rcpid2[1000];
    PpChain->SetBranchAddress("pid2", rcpid2);
    double rcz[1000];
    PpChain->SetBranchAddress("z", rcz);
    double rcepair[1000];
    PpChain->SetBranchAddress("epair", rcepair);
    double rcnef[1000];
    PpChain->SetBranchAddress("nef", rcnef);*/

    //! Output and histograms
    TH2D *h_eta_vs_phi = new TH2D("h", ";phi;eta", 100, -3.14, 3.14, 100, -0.6, 0.6);
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

        int ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
        if (ppevi < 0)
            continue;

        PpChain->GetEntry(ppevi);
        for (int j = 0; j < ppnjets; ++j)
        {
            TStarJetVectorJet *ppjet = (TStarJetVectorJet *)PpJets->At(j);
            h_eta_vs_phi->Fill(ppjet->Phi(), ppjet->Eta());//, 212563/3284499);
            //cout << j << ", " << ppjet->Phi() << ", " << ppjet->Eta() << endl;
        }
    }

    TCanvas *c = new TCanvas("c", "c", 1000, 600);
    //c->SetLogz(1);
    h_eta_vs_phi->GetZaxis()->SetRangeUser(0, 400);
    h_eta_vs_phi->Draw("colz");
    c->SaveAs("plots/test.png");
    return missed;
}