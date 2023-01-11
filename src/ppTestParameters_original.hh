/* @file parameters file
   @author Kolja Kauder
   @brief Common parameters
   @details Used to quickly include the same parameters into different macros.
   @date Mar 23, 2017   
 */

#ifndef PPTESTPARAMETERS_HH
#define PPTESTPARAMETERS_HH

#include "JetAnalyzer.hh"

using fastjet::JetAlgorithm;
using fastjet::antikt_algorithm;
using fastjet::kt_algorithm;
using fastjet::cambridge_algorithm;

enum INTYPE{ MCTREE, INTREE, INPICO, MCPICO, HERWIGTREE };

// Return values for the main routine
enum class EVENTRESULT{
  PROBLEM,
  ENDOFINPUT,
  NOTACCEPTED,
  NOCONSTS,
  NOJETS,
  JETSFOUND
};

class ppTestParameters{

public :
  double R = 0.4;            ///< Resolution parameter ("radius").

  /// Jet algorithm for the original jets
  JetAlgorithm LargeJetAlgorithm=fastjet::antikt_algorithm;
  // JetAlgorithm LargeJetAlgorithm = fastjet::cambridge_algorithm;
  
  /// Alternative reclustering to Cambridge/Aachen (at our own risk)
  bool CustomRecluster=false;
  JetAlgorithm CustomReclusterJetAlgorithm;
  // bool CustomRecluster=true;
  // JetAlgorithm CustomReclusterJetAlgorithm = fastjet::antikt_algorithm;
  bool Recursive=false;  ///< Repeat on subjets?

  /// Repetitions in the background. Anything other than 1 WILL NOT WORK because
  /// a) we're using explicit ghosts (though we don't have to)
  /// b) more importantly, the background subtractor contains fastjet::SelectorNHardest(2)
  ///    which doesn't work jet-by-jet and throws an error
  int GhostRepeat = 1;
  float GhostArea = 0.005;    ///< ghost area

  // int ghost_repeat = 1;
  // double ghost_area = 0.01;    ///< ghost area
  // const double ghost_area = 0.0005;    ///< ghost area

  // const double PtJetMin = 20.0;    ///< Min jet pT
  double PtJetMin = 5.0;      ///< Min jet pT
  double PtJetMax = 1000.0;   ///< Max jet pT
  double MJetMin = 0.0;
  double LeadPtMin=5.0;                 ///< leading jet minimum p<SUB>T</SUB>
    
  double MaxJetNEF=0.9;       ///< Max neutral energy fraction

  double EtaConsCut = 1.0;    ///< Constituent |&eta;| acceptance
  double PtConsMin=0.2;       ///< Constituent pT minimum
  // const double PtConsMin=2.0;       ///< Constituent pT minimum
  double PtConsMax=30;        ///< Constituent pT maximum
  
  double RefMultCut=0;        ///< Reference multiplicity. Needs to be rethought to accomodate pp and AuAu
  
  double VzCut=30;            ///< Vertex z 
  // const double VzDiffCut=6;         ///< |Vz(TPC) - Vz(VPD)| <-- NOT WORKING in older data (no VPD)
  double VzDiffCut=99999;      ///< |Vz(TPC) - Vz(VPD)|
  
  double DcaCut=1.0;          ///< track dca
  double NMinFit=20;             ///< minimum number of fit points for tracks
  double FitOverMaxPointsCut=0.52; ///< NFit / NFitPossible

  double HadronicCorr = 0.9999; ///< Fraction of hadronic correction

  double FakeEff = 1.0; ///< fake efficiency for systematics. 0.95 is a reasonable example.

  Int_t IntTowScale=0;
  /// Tower GAIN: 4.8%
  Float_t fTowUnc=0.038;
  /// Tower scale for uncertainty;
  float fTowScale=1.0;
  
  // ************************************
  // Do NOT cut high tracks and towers!
  // Instead, reject the whole event when
  // of these is found
  // ************************************
  double MaxEtCut=1000;       ///< tower ET cut
  double MaxTrackPt=1000;     ///< track pT cut

  
  // EVENT rejection cuts
  double MaxEventPtCut=30;       ///< max track pT cut for event
  double MaxEventEtCut=30;       ///< max tower ET cut for event
  double MinEventEtCut=0;        ///< min event ET cut for event
  double ManualHtCut=0.0;        ///< necessary for some embedding picos. Should always equal MinEventEtCut


  // Geant files have messed up runid and event id.
  // Switch for fixing. Should be turned on by default for Geant files, off otherwise.
  bool UseGeantNumbering=false;


  // GEANT12 example
  // ---------------
  TString InputName = "Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt25_35_g0.root"; ///< Input name
  INTYPE intype = INPICO;             ///< Input type (can be a pico dst, a result tree, an MC tree)
  TString ChainName = "JetTree";      ///< Name of the input chain
  TString TriggerName = "ppJP2";        ///< Trigger type (All, MB, HT, pp, ppHT, ppJP2, ppJP(careful!))

  // Only putting it here so that it can be initialized in the analysis class and then used
  TString OutFileName = "TmpResult.root";     ///< Output file


};
#endif // PPTESTPARAMETERS_HH
