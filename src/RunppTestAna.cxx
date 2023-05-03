/* @file RunppTestAna.cxx
   @author Raghav Kunnawalkam Elayavalli
   @version Revision 1.0
   @brief Test Analysis for Run12 jets in both data and PYTHIA 6 embedding
   @date March 16, 2022
*/

#include "ppTestParameters.hh"
#include "ppTestAnalysis.hh"
#include "TStarJetVectorJet.h"

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TFile.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TParameter.h>
#include "TString.h"
#include "TObjString.h"

#include <set>
#include <vector>
#include <algorithm>

#include <cmath>
#include <climits>

#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <exception>

using namespace std;
using namespace fastjet;
using namespace contrib;

// Mostly for run 12
bool readinbadrunlist(vector<int> &badrun, TString csvfile);

int main(int argc, const char **argv)
{

  ppTestAnalysis *ppana = nullptr;
  try
  {
    ppana = new ppTestAnalysis(argc, argv);
  }
  catch (std::exception &e)
  {
    cerr << "Initialization failed with exception " << e.what() << endl;
    return -1;
  }

  if (ppana->InitChains() == false)
  {
    cerr << "Chain initialization failed" << endl;
    return -1;
  }

  // Get parameters we used
  // ----------------------
  const ppTestParameters pars = ppana->GetPars();

  // Explicitly choose bad tower list here
  // -------------------------------------
  shared_ptr<TStarJetPicoReader> pReader = ppana->GetpReader();
  if (pReader)
  {
    TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
    towerCuts->AddBadTowers("/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Combined_pp200Y12_badtower_isaac.list");
  }

  // Explicitly add bad run list here
  // --------------------------------
  if (pReader)
  {
    TString csvfile = "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/pp200Y12_badrun_isaac.list";
    // TString csvfile="/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/pp200Y12_badrun.list";
    vector<int> badruns;
    if (readinbadrunlist(badruns, csvfile) == false)
    {
      cerr << "Problems reading bad run list" << endl;
      return -1;
    }
    pReader->AddMaskedRuns(badruns);
  }

  // Files and histograms
  // --------------------
  TFile *fout = new TFile(pars.OutFileName, "RECREATE");
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  TH3::SetDefaultSumw2(true);

  // List of miscellaneous info
  // --------------------------
  TTree *info = new TTree("info", "Information");
  info->Branch("InputName", (void *)pars.InputName.Data(), "InputName/C");
  info->Branch("ChainName", (void *)pars.ChainName.Data(), "ChainName/C");
  info->Branch("R", (void *)&pars.R, "R/D");
  info->Branch("LargeJetAlgorithm", (UInt_t *)&pars.LargeJetAlgorithm, "LargeJetAlgorithm/i");
  info->Branch("PtJetMin", (void *)&pars.PtJetMin, "PtJetMin/D");
  info->Branch("PtJetMax", (void *)&pars.PtJetMax, "PtJetMax/D");
  info->Branch("EtaConsCut", (void *)&pars.EtaConsCut, "EtaConsCut/D");
  info->Branch("PtConsMin", (void *)&pars.PtConsMin, "PtConsMin/D");
  info->Branch("PtConsMax", (void *)&pars.PtConsMax, "PtConsMax/D");

  // Save results
  // ------------
  TTree *ResultTree = new TTree("ResultTree", "Result Jets");

  // Give each event a unique ID to compare event by event with different runs
  int runid;
  ResultTree->Branch("runid", &runid, "runid/I");
  int runid1;
  ResultTree->Branch("runid1", &runid1, "runid1/I");
  int eventid;
  ResultTree->Branch("eventid", &eventid, "eventid/I");
  double weight = 1;
  ResultTree->Branch("weight", &weight, "weight/D");
  double refmult;
  ResultTree->Branch("refmult", &refmult, "refmult/D");
  int njets = 0;
  ResultTree->Branch("njets", &njets, "njets/I");
  double vz;
  ResultTree->Branch("vz", &vz, "vz/D");
  int highw;
  ResultTree->Branch("highw", &highw, "highw/I");
  int mult;
  ResultTree->Branch("mult", &mult, "mult/I");

  TClonesArray Jets("TStarJetVectorJet");
  ResultTree->Branch("Jets", &Jets);
  double nef[1000];
  ResultTree->Branch("nef", nef, "nef[njets]/D");
  double pt[1000];
  ResultTree->Branch("pt", pt, "pt[njets]/D");
  double m[1000];
  ResultTree->Branch("m", m, "m[njets]/D");
  int b[1000]; // -1 if the 2 leading tracks are of opposite charge, 0 if any is neutral
  ResultTree->Branch("b", b, "b[njets]/I");
  double dr[1000];
  ResultTree->Branch("dr", dr, "dr[njets]/D");
  double q0[1000];
  ResultTree->Branch("q0", q0, "q0[njets]/D");
  double q2[1000];
  ResultTree->Branch("q2", q2, "q2[njets]/D");
  double epair[1000];
  ResultTree->Branch("epair", epair, "epair[njets]/D");
  double z[1000];
  ResultTree->Branch("z", z, "z[njets]/D");
  double rg[1000];
  ResultTree->Branch("rg", rg, "rg[njets]/D");
  double zg[1000];
  ResultTree->Branch("zg", zg, "zg[njets]/D");
  double mg[1000];
  ResultTree->Branch("mg", mg, "mg[njets]/D");
  int n[1000];
  ResultTree->Branch("n", n, "n[njets]/I");
  int nch[1000];
  ResultTree->Branch("nch", nch, "nch[njets]/I");
  int index[1000];
  ResultTree->Branch("index", index, "index[njets]/I");
  int reject[1000];
  ResultTree->Branch("reject", reject, "reject[njets]/I");
  int pid1[1000];
  ResultTree->Branch("pid1", pid1, "pid1[njets]/I");
  int pid2[1000];
  ResultTree->Branch("pid2", pid2, "pid2[njets]/I");
  double pt1[1000];
  ResultTree->Branch("pt1", pt1, "pt1[njets]/D");
  double pt2[1000];
  ResultTree->Branch("pt2", pt2, "pt2[njets]/D");

  // Helpers
  TStarJetVector *sv;
  TObjString *tobjs;

  // Go through events
  // -----------------
  Long64_t Naccepted = 0;
  cout << "Running analysis" << endl;
  try
  {
    bool ContinueReading = true;

    while (ContinueReading)
    {

      Jets.Clear();
      refmult = 0;
      runid = -(INT_MAX - 1);
      eventid = -(INT_MAX - 1);

      EVENTRESULT ret = ppana->RunEvent();

      // Understand what happened in the event
      switch (ret)
      {
      case EVENTRESULT::PROBLEM:
        cerr << "Encountered a serious issue" << endl;
        return -1;
        break;
      case EVENTRESULT::ENDOFINPUT:
        cout << "End of Input" << endl;
        ContinueReading = false;
        continue;
        break;
      case EVENTRESULT::NOTACCEPTED:
        cout << "Event rejected" << endl;
        // continue;
        break;
      case EVENTRESULT::NOCONSTS:
        // cout << "Event empty." << endl;
        // continue;
        break;
      case EVENTRESULT::NOJETS:
        // cout << "No jets found." << endl;
        // continue;
        break;
      case EVENTRESULT::JETSFOUND:
        // The only way not to break out or go back to the top
        // Do Something
        // cout << "Jets found." << endl;
        Naccepted++;
        break;
      default:
        cerr << "Unknown return value." << endl;
        return -1;
        // I understand that a break after continue or return is silly...
        // But it's necessary in nested switches in root and I don't want to lose the habit
        break;
      }

      // Now we can pull out details and results
      // ---------------------------------------
      weight = ppana->GetEventWeight();
      refmult = ppana->GetRefmult();
      runid = ppana->GetRunid();
      runid1 = ppana->GetRunid1();
      eventid = ppana->GetEventid();
      vz = ppana->GetVz();
      highw = ppana->CheckHighW();
      mult = ppana->GetEventMult();

      vector<ResultStruct> Result = ppana->GetResult();

      njets = Result.size();
      int ijet = 0;
      if (njets == 0)
      {
        ResultTree->Fill();
        continue;
      }
      for (auto &gr : Result)
      {

        TStarJetVector sv = TStarJetVector(MakeTLorentzVector(gr.orig));
        sv.SetCharge(gr.orig.user_info<JetAnalysisUserInfo>().GetQuarkCharge() / 3);
        new (Jets[ijet]) TStarJetVectorJet(sv);
        nef[ijet] = gr.orig.user_info<JetAnalysisUserInfo>().GetNumber();
        pt[ijet] = gr.orig.perp();
        m[ijet] = gr.orig.m();
        q0[ijet] = gr.q0;
        q2[ijet] = gr.q2;
        b[ijet] = gr.b;
        dr[ijet] = gr.dr;
        epair[ijet] = gr.epair;
        z[ijet] = gr.z;
        rg[ijet] = gr.sd.structure_of<fastjet::contrib::SoftDrop>().delta_R();
        zg[ijet] = gr.sd.structure_of<fastjet::contrib::SoftDrop>().symmetry();
        mg[ijet] = gr.sd.m();
        n[ijet] = gr.orig.constituents().size();
        nch[ijet] = gr.nch;
        index[ijet] = ijet;
        reject[ijet] = gr.reject;
        pid1[ijet] = gr.pid1;
        pid2[ijet] = gr.pid2;
        pt1[ijet] = gr.pt1;
        pt2[ijet] = gr.pt2;

        ijet++;
      }

      ResultTree->Fill();
    }
  }
  catch (std::string &s)
  {
    cerr << "RunEvent failed with string " << s << endl;
    return -1;
  }
  catch (std::exception &e)
  {
    cerr << "RunEvent failed with exception " << e.what() << endl;
    return -1;
  }

  info->Branch("Naccepted", &Naccepted, "Naccepted/L");
  info->Fill();

  fout->Write();

  cout << "Done." << endl;

  delete ppana;
  return 0;
}

//----------------------------------------------------------------------
bool readinbadrunlist(vector<int> &badrun, TString csvfile)
{

  // open infile
  std::string line;
  std::ifstream inFile(csvfile);

  std::cout << "Loading bad run id from " << csvfile.Data() << std::endl;
  ;

  if (!inFile.good())
  {
    std::cout << "Can't open " << csvfile.Data() << std::endl;
    return false;
  }

  while (std::getline(inFile, line))
  {
    if (line.size() == 0)
      continue; // skip empty lines
    if (line[0] == '#')
      continue; // skip comments

    std::istringstream ss(line);
    while (ss)
    {
      std::string entry;
      std::getline(ss, entry, ',');
      int ientry = atoi(entry.c_str());
      if (ientry)
      {
        badrun.push_back(ientry);
        // std::cout<<"Added bad runid "<<ientry<<std::endl;
      }
    }
  }

  return true;
}