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

//! ----------------------------------------------------
int GetMcevi_new(string McFile, string PpFile)
{
    TString dir = "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/";

    TFile *Mcf = new TFile(dir + "Results/pythia/" + McFile);
    TTree *McChain = (TTree *)Mcf->Get("ResultTree");
    McChain->BuildIndex("runid", "eventid");

    int mcrunid;
    McChain->SetBranchAddress("runid", &mcrunid);
    int mceventid;
    McChain->SetBranchAddress("eventid", &mceventid);
    float mcpttot;
    McChain->SetBranchAddress("pttot", &mcpttot);

    TFile *Ppf = new TFile(dir + "Results/geant/" + PpFile);
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
    //float pppttot;
    //PpChain->SetBranchAddress("pttot", &ppptot);

    cout << "Number of Pythia events: " << McChain->GetEntries() << endl;
    cout << "Number of Geant events:  " << PpChain->GetEntries() << endl;

    // read bad run list
    vector<int> runid1_list;
    FILE *file = fopen(dir + "pp200Y12_badrun_isaac.list", "r");
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
    vector<float> pppt_list;
    for (Long64_t ppevi = 0; ppevi < PpChain->GetEntries(); ++ppevi)
    {
        PpChain->GetEntry(ppevi);
        // check if pprunid1 is in bad runid1 list
        if (find(runid1_list.begin(), runid1_list.end(), pprunid1) != runid1_list.end() || abs(ppvz) > 30 || pphighw > 0)
        {
            ppeventid_list.push_back(ppeventid);
            pprunid_list.push_back(pprunid);
        }
        /*if (find(pppt_list.begin(), pppt_list.end(), pppttot))
        {
            ppeventid_list.push_back(ppeventid);
            pprunid_list.push_back(pprunid);
        }
        else
        {
            pppt_list.push_back(pppttot);
        }*/
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
    vector<int> mcevi_list3;
    set<float> mcpt_list;
    //vector<float> mceta_list;
    //vector<float> mcphi_list;

    for (Long64_t mcevi = 0; mcevi < McChain->GetEntries(); ++mcevi)
    {
        McChain->GetEntry(mcevi);
        int ppevi = PpChain->GetEntryNumberWithIndex(mcrunid, mceventid);
        if (ppevi < 0) // event not embedded
        {
            mcevi_list2.push_back(mcevi);
        }
        if (mcpt_list.find(mcpttot) != mcpt_list.end())// && count(mceta_list.begin(), mceta_list.end(), eta) && count(mcphi_list.begin(), mcphi_list.end(), phi))
        {
            mcevi_list3.push_back(mcevi); // some events have identical total particle pT
        }
        else
        {
            mcpt_list.insert(mcpttot);
            //mceta_list.push_back(eta);
            //mcphi_list.push_back(phi);
        }
    }
    cout << "Number of Pythia events that are not embedded: " << mcevi_list2.size() << endl;
    cout << "Number of Pythia events that are repeated: " << mcevi_list3.size() << endl;

    ofstream out_file(dir + "mcevi_new_new.list");
    ostream_iterator<int> output_iterator(out_file, "\n");
    copy(mcevi_list.begin(), mcevi_list.end(), output_iterator);

    ofstream out_file2(dir + "mcevi2_new_new.list");
    ostream_iterator<int> output_iterator2(out_file2, "\n");
    copy(mcevi_list2.begin(), mcevi_list2.end(), output_iterator2);

    ofstream out_file3(dir + "mcevi3_new_new.list");
    ostream_iterator<int> output_iterator3(out_file3, "\n");
    copy(mcevi_list3.begin(), mcevi_list3.end(), output_iterator3);

    return 0;
}