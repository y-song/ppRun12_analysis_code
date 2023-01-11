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
bool readinbadrunlist(vector<int> & badrun, TString csvfile);
  
// DEBUG
void decluster (PseudoJet j);

/*
    - Parse parameters
    - Set up input tree
    - Set up output histos and tree
    - Initialize SubjetAnalysis object
    - If needed, combine input from two sources
    - Loop through events
    \arg argv: flags.
    Display options with
    <BR><tt>% UnifiedSubjetWrapper -h </tt> 
    <BR>Note that wildcarded file patterns should be in single quotes.
*/

int main( int argc, const char** argv ){

  ppTestAnalysis* ppana = nullptr; 
  try {
    ppana = new ppTestAnalysis(argc, argv );
  } catch ( std::exception& e ){
    cerr << "Initialization failed with exception " << e.what() << endl;
    return -1;
  }
  
  if ( ppana->InitChains() == false ){
    cerr << "Chain initialization failed" << endl;
    return -1;
  }
  
  // Get parameters we used
  // ----------------------
  const ppTestParameters pars  = ppana->GetPars();

  // Explicitly choose bad tower list here
  // -------------------------------------
  // Otherwise too easy to hide somewhere and forget...
  shared_ptr<TStarJetPicoReader> pReader = ppana->GetpReader();
  
  if ( pReader ){
    TStarJetPicoTowerCuts* towerCuts = pReader->GetTowerCuts();
    towerCuts->AddBadTowers( "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/Combined_pp200Y12_badtower_isaac.list");
    //towerCuts->AddBadTowers("/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/pp200Y12_badtower.list");
  }
  
  // Explicitly add bad run list here
  // --------------------------------
  if ( pReader ){
      TString csvfile= "/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/pp200Y12_badrun_isaac.list";
      //TString csvfile="/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/pp200Y12_badrun.list";
      vector<int> badruns;
      if ( readinbadrunlist( badruns, csvfile) == false ){
	cerr << "Problems reading bad run list" << endl;
	return -1;
      }
      pReader->AddMaskedRuns (badruns);
  }
  

  // Files and histograms
  // --------------------
  TFile* fout = new TFile( pars.OutFileName, "RECREATE");
  
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  TH3::SetDefaultSumw2(true);
    
  TH2D* cphieta = new TH2D("cphieta","", 100, 0, TMath::TwoPi(), 100, -1, 1);
  TH2D* nphieta = new TH2D("nphieta","", 100, 0, TMath::TwoPi(), 100, -1, 1);
  
  // List of miscellaneous info
  // --------------------------
  TTree* info = new TTree("info", "Information");
  info->Branch("InputName"         , (void*)pars.InputName.Data()     , "InputName/C" );
  info->Branch("ChainName"         , (void*)pars.ChainName.Data()     , "ChainName/C" );
  info->Branch("R"                 , (void*) &pars.R                          , "R/D" );
  info->Branch("LargeJetAlgorithm" , (UInt_t*)&pars.LargeJetAlgorithm , "LargeJetAlgorithm/i" );
  info->Branch("PtJetMin"          , (void*) &pars.PtJetMin                   , "PtJetMin/D" );
  info->Branch("PtJetMax"          , (void*) &pars.PtJetMax                   , "PtJetMax/D" );
  info->Branch("EtaConsCut"        , (void*) &pars.EtaConsCut                 , "EtaConsCut/D" );
  info->Branch("PtConsMin"         , (void*) &pars.PtConsMin                  , "PtConsMin/D" );
  info->Branch("PtConsMax"         , (void*) &pars.PtConsMax                  , "PtConsMax/D" );


  // Save results
  // ------------
  TTree* ResultTree=new TTree("ResultTree","Result Jets");

  // Give each event a unique ID to compare event by event with different runs
  int runid;
  ResultTree->Branch("runid",   &runid, "runid/I");
  int runid1;
  ResultTree->Branch("runid1", &runid1, "runid1/I");
  int eventid;
  ResultTree->Branch("eventid", &eventid, "eventid/I");
  double weight=1;
  ResultTree->Branch("weight",  &weight, "weight/D" );      
  double refmult;
  ResultTree->Branch("refmult",&refmult, "refmult/D");  
  int njets=0;
  ResultTree->Branch("njets",   &njets, "njets/I" );
  double vz;
  ResultTree->Branch("vz", &vz, "vz/D");
  int highw;
  ResultTree->Branch("highw", &highw, "highw/I");
  int mult;
  ResultTree->Branch("mult", &mult, "mult/I");

  TClonesArray Jets( "TStarJetVectorJet" );
  ResultTree->Branch("Jets", &Jets );
  double nef[1000];
  ResultTree->Branch("nef",  nef, "nef[njets]/D" );
  double pt[1000];
  ResultTree->Branch("pt", pt, "pt[njets]/D");
  double m[1000];
  ResultTree->Branch("m", m, "m[njets]/D");
  double q[1000];
  ResultTree->Branch("q", q, "q[njets]/D");
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
  int nue[1000];
  ResultTree->Branch("nue", nue, "nue[njets]/I");
  double ptue[1000];
  ResultTree->Branch("ptue", ptue, "ptue[njets]/D");
  int reject[1000];
  ResultTree->Branch("reject", reject, "reject[njets]/I");

  // Helpers
  TStarJetVector* sv;
  TObjString* tobjs;
  
  // Go through events
  // -----------------
  Long64_t Naccepted=0;
  cout << "Running analysis" << endl;
  try {
    bool ContinueReading = true;

    while ( ContinueReading ){

      Jets.Clear();
      refmult=0;
      runid   =-(INT_MAX-1);
      eventid =-(INT_MAX-1);

      EVENTRESULT ret=ppana->RunEvent();

      // Understand what happened in the event
      switch (ret){
      case  EVENTRESULT::PROBLEM :
	cerr << "Encountered a serious issue" << endl;
	return -1;
	break;	
      case  EVENTRESULT::ENDOFINPUT :
	cout << "End of Input" << endl;
	ContinueReading=false;
	continue;
	break;
      case  EVENTRESULT::NOTACCEPTED :
	cout << "Event rejected" << endl;
	//continue;
	break;
      case  EVENTRESULT::NOCONSTS :
	//cout << "Event empty." << endl;
	//continue;
	break;
      case  EVENTRESULT::NOJETS :
	//cout << "No jets found." << endl;
	//continue;
	break;
      case  EVENTRESULT::JETSFOUND:
	// The only way not to break out or go back to the top
	// Do Something
	//cout << "Jets found." << endl;
	Naccepted++;
	break;
      default :
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
      // softdrop
      // fastjet::contrib::SoftDrop sd(0, 0.1);
      // sdjet = sd(Result);	
	     
      njets=Result.size();
      int ijet=0;
      if (njets == 0)
      {
	ResultTree->Fill();
	continue;
      }
      for ( auto& gr : Result ){	

	TStarJetVector sv = TStarJetVector( MakeTLorentzVector( gr.orig) );
	sv.SetCharge( gr.orig.user_info<JetAnalysisUserInfo>().GetQuarkCharge() / 3 );      
	new ( Jets[ijet] ) TStarJetVectorJet ( sv );
	nef[ijet] = gr.orig.user_info<JetAnalysisUserInfo>().GetNumber() ;
	pt[ijet] = gr.orig.perp();
	m[ijet] = gr.orig.m();
	q[ijet] = gr.q;
	rg[ijet] = gr.sd.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	zg[ijet] = gr.sd.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	mg[ijet] = gr.sd.m();
	n[ijet] = gr.orig.constituents().size();
	nch[ijet] = gr.nch;
	index[ijet] = ijet;
	nue[ijet] = gr.nue;
	ptue[ijet] = gr.ptue;
	reject[ijet] = gr.reject;

	ijet++;
      }
      
      ResultTree->Fill();

      const vector<PseudoJet>& particles = ppana->GetParticles();
      for (auto& p : particles ){
      if ( p.has_user_info<JetAnalysisUserInfo>() ){
	          if ( abs( p.user_info<JetAnalysisUserInfo>().GetQuarkCharge() )>0 ){
	            cphieta->Fill(p.phi(),p.eta());
//	            track_pT.push_back(p.pt());
//	            track_phi.push_back(p.phi());
//	            track_eta.push_back(p.eta());
//	            track_id.push_back(0);
     		  }
     		  else {
	            nphieta->Fill(p.phi(),p.eta());
//	            tower_pT.push_back(p.pt());
//	            tower_phi.push_back(p.phi());
//	            tower_eta.push_back(p.eta());
//	            tower_id.push_back(0);
     		  }
      }
      else {
		  cout << "constitnuent: no charge info" << endl;
      }
      }
 } } catch (std::string& s) {
    cerr << "RunEvent failed with string " << s << endl;
    return -1;
  } catch ( std::exception& e ){
    cerr << "RunEvent failed with exception " << e.what() << endl;
    return -1;
  }

    
  info->Branch("Naccepted" , &Naccepted , "Naccepted/L" );
  info->Fill();
  
  fout->Write();

  cout << "Done." << endl;

  delete ppana;
  return 0;
}
//----------------------------------------------------------------------

// DEBUG
void decluster (PseudoJet j){
  cout << " ### Declustering ### " << endl;
  cout << j.pt() << endl;

  PseudoJet piece1, piece2;

  static int level=0;
  if ( j.has_parents(piece1, piece2) ){
    cout << "level = " << level++ << endl;
    cout << "piece1.pt() = " << piece1.pt() << endl;
    cout << "piece2.pt() = " << piece2.pt() << endl;

    if (! piece1.is_pure_ghost() ) decluster(piece1);
    if (! piece2.is_pure_ghost() ) decluster(piece2);

  } else cout << " Done with this branch" << endl;
  return;
}

//----------------------------------------------------------------------
bool readinbadrunlist(vector<int> & badrun, TString csvfile) {
	
  // open infile
  std::string line;
  std::ifstream inFile (csvfile );
	
  std::cout<<"Loading bad run id from "<< csvfile.Data()<<std::endl;;
	        
  if ( !inFile.good() ) {
    std::cout<<"Can't open "<<csvfile.Data()<<std::endl;
    return false;
  }
	
  while (std::getline (inFile, line) ){
    if ( line.size()==0 ) continue; // skip empty lines
    if ( line[0] == '#' ) continue; // skip comments
	
    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if (ientry) {
	badrun.push_back( ientry );
	//std::cout<<"Added bad runid "<<ientry<<std::endl;
      }
    }
  }
	
  return true;
}

