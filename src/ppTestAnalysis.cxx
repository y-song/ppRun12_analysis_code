/* @file ppTestAnalysis.cxx
    @author Raghav Kunnawalkam Elayavalli
    @version Revision 1.0
    @brief test pp analysis for Run12 data and embedding
    @details Uses JetAnalyzer objects
    @date March 16, 2022
*/

#include "ppTestAnalysis.hh"
#include <stdlib.h> // for getenv, atof, atoi
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

using std::cerr;
using std::cout;
using std::endl;

// Standard ctor
ppTestAnalysis::ppTestAnalysis(const int argc, const char **const argv)
{
  // Parse arguments
  // ---------------
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay = true;
  bool forcedgeantnum = false;
  NEvents = -1;
  for (auto parg = arguments.begin(); parg != arguments.end(); ++parg)
  {
    string arg = *parg;
    if (arg == "-R")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.R = atof(parg->data());
    }
    else if (arg == "-lja")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.LargeJetAlgorithm = AlgoFromString(*parg);
    }
    else if (arg == "-pj")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.PtJetMin = atof((parg)->data());
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.PtJetMax = atof((parg)->data());
    }
    else if (arg == "-mj")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.MJetMin = atof((parg)->data());
    }
    else if (arg == "-ec")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.EtaConsCut = atof((parg)->data());
    }
    else if (arg == "-pc")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.PtConsMin = atof((parg)->data());
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.PtConsMax = atof((parg)->data());
    }
    else if (arg == "-hadcorr")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.HadronicCorr = atof(parg->data());
    }
    else if (arg == "-o")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.OutFileName = *parg;
    }
    else if (arg == "-i")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.InputName = *parg;
    }
    else if (arg == "-c")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.ChainName = *parg;
    }
    else if (arg == "-trig")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.TriggerName = *parg;
    }
    else if (arg == "-intype")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      if (*parg == "pico")
      {
        pars.intype = INPICO;
        continue;
      }
      if (*parg == "mcpico")
      {
        pars.intype = MCPICO;
        continue;
      }
      if (*parg == "mctree")
      {
        pars.intype = MCTREE;
        continue;
      }
      argsokay = false;
      break;
    }
    else if (arg == "-N")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      NEvents = atoi(parg->data());
    }
    else if (arg == "-tracksmear")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      switch (atoi(parg->data()))
      {
      case 0:
        cout << "Track Smearing disabled." << endl;
        break;
      case 1:
        cout << "Track Smearing set to pp primaries. DeltapT/pT = 0.01+0.005 pT " << endl;
        SigmaPt = new TF1("SigmaPt", "[0] + [1]*x", 0, 100);
        SigmaPt->FixParameter(0, 0.01);
        SigmaPt->FixParameter(1, 0.005);
        break;
      case 2:
        cout << "Track Smearing set to AA primaries. DeltapT/pT = 0.005+0.0025 pT " << endl;
        SigmaPt = new TF1("SigmaPt", "[0] + [1]*x", 0, 100);
        SigmaPt->FixParameter(0, 0.005);
        SigmaPt->FixParameter(1, 0.0025);

        break;
      case 3:
        cout << "Track Smearing set to AA globals not implemented. DeltapT/pT = 0.01 * pT^2" << endl;
        SigmaPt = new TF1("SigmaPt", "[0]*x*x", 0, 100);
        SigmaPt->FixParameter(0, 0.01);
        break;
      case 4:
        cout << "Unrecognized pT smearing option. " << endl;
        argsokay = false;
        break;
      }
    }
    else if (arg == "-fakeeff")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.FakeEff = atof(parg->data());
      cout << "Setting fake efficiency to " << pars.FakeEff << endl;
      if (pars.FakeEff < 0 || pars.FakeEff > 1)
      {
        argsokay = false;
        break;
      }
    }
    else if (arg == "-towunc")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.IntTowScale = atoi(parg->data());
      pars.fTowScale = 1.0 + pars.IntTowScale * pars.fTowUnc;
      cout << "Setting tower scale to " << pars.fTowScale << endl;
      if (pars.IntTowScale < -1 || pars.FakeEff > 1)
      {
        argsokay = false;
        break;
      }
    }
    else if (arg == "-geantnum")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      if (*parg != "0" && *parg != "1")
      {
        argsokay = false;
        break;
      }
      forcedgeantnum = true;
      pars.UseGeantNumbering = bool(atoi(parg->data()));
    }
    else if (arg == "-jetnef")
    {
      if (++parg == arguments.end())
      {
        argsokay = false;
        break;
      }
      pars.MaxJetNEF = atof(parg->data());
      cout << "Setting Max Jet NEF to " << pars.MaxJetNEF << endl;
      if (pars.MaxJetNEF < 0 || pars.MaxJetNEF > 1)
      {
        argsokay = false;
        break;
      }
    }
    else
    {
      argsokay = false;
      break;
    }
  }

  if (!argsokay)
  {
    cerr << "usage: " << argv[0] << endl
         << " [-o OutFileName]" << endl
         << " [-N Nevents (<0 for all)]" << endl
         << " [-R radius]" << endl
         << " [-lja LargeJetAlgorithm]" << endl
         << " [-i infilepattern]" << endl
         << " [-c chainname]" << endl
         << " [-intype pico|mcpico(embedding)|tree|mctree(pythia)|herwigtree]" << endl
         << " [-trig trigger name (e.g. HT)]" << endl
         << " [-pj PtJetMin PtJetMax]" << endl
         << " [-ec EtaConsCut]" << endl
         << " [-pc PtConsMin PtConsMax]" << endl
         << " [-hadcorr HadronicCorrection]  -- Set to a negative value for MIP correction." << endl
         << " [-psc PtSubConsMin PtSubConsMax]" << endl
         << " [-tracksmear number] -- enable track pT smearing. " << endl
         << "                      -- 1: pp primaries, 2: AuAu primaries, 3: AuAu (maybe also pp) globals." << endl
         << " [-fakeeff (0..1)] -- enable fake efficiency for systematics. 0.95 is a reasonable example." << endl
         << " [-towunc -1|0|1 ] -- Shift tower energy by this times " << pars.fTowUnc << endl
         << " [-geantnum true|false] (Force Geant run event id hack)" << endl
         << endl
         << endl
         << "NOTE: Wildcarded file patterns should be in single quotes." << endl
         << endl;
    throw std::runtime_error("Not a valid list of options");
  }

  // Consistency checks
  // ------------------

  if (pars.PtJetMin <= 0)
  {
    throw std::runtime_error("PtJetMin needs to be positive (0.001 will work).");
  }

  if (pars.intype == MCPICO)
  {
    // This refers to the JetTreeMc branch
    if (pars.ChainName != "JetTreeMc")
    {
      throw std::runtime_error("Unsuitable chain name for pars.intype==MCPICO");
    }
  }
  if (pars.ChainName == "JetTreeMc")
  {
    if (pars.intype != MCPICO)
    {
      throw std::runtime_error("Unsuitable chain name for pars.intype==MCPICO");
    }
  }

  // Derived rapidity cuts
  // ---------------------
  EtaJetCut = pars.EtaConsCut - pars.R;

  // Jet candidate selectors
  // -----------------------
  // select_jet_eta = SelectorAbsRapMax( EtaJetCut );
  select_jet_eta = SelectorAbsEtaMax(EtaJetCut);
  select_jet_pt = SelectorPtRange(pars.PtJetMin, pars.PtJetMax);
  select_jet_m = SelectorMassMin(pars.MJetMin);
  if (pars.intype == MCPICO)
  {
    select_jet = select_jet_eta * select_jet_pt;
  }
  if (pars.intype == INPICO)
  {
    select_jet = select_jet_eta * select_jet_pt; // * select_jet_m;
  }

  // Repeat on subjets?
  // ------------------
  pars.Recursive = pars.InputName.Contains("Pythia") && false;

  // Initialize jet finding
  // ----------------------
  JetDef = JetDefinition(pars.LargeJetAlgorithm, pars.R);
  cout << " R = " << pars.R << endl;
  cout << " Original jet algorithm : " << pars.LargeJetAlgorithm << endl;
  cout << " PtJetMin = " << pars.PtJetMin << endl;
  cout << " PtJetMax = " << pars.PtJetMax << endl;
  cout << " MJetMin = " << pars.MJetMin << endl;
  cout << " PtConsMin = " << pars.PtConsMin << endl;
  cout << " PtConsMax = " << pars.PtConsMax << endl;
  cout << " Constituent eta cut = " << pars.EtaConsCut << endl;
  cout << " Jet eta cut = " << EtaJetCut << endl;
  cout << " Ghosts out to eta = " << EtaGhostCut << endl;
  cout << " Reading tree named \"" << pars.ChainName << "\" from " << pars.InputName << endl;
  cout << " intype = " << pars.intype << endl;
  cout << " Writing to " << pars.OutFileName << endl;
  cout << " ----------------------------" << endl;

  // Provide a gaussian for track pt smearing
  // ----------------------------------------
  SmearPt = new TF1("SmearPt", "gaus(0)", -1, 1);
}
//----------------------------------------------------------------------
ppTestAnalysis::~ppTestAnalysis()
{
  if (pJA)
  {
    delete pJA;
    pJA = 0;
  }
}
//----------------------------------------------------------------------
bool ppTestAnalysis::InitChains()
{

  // For trees of TStarJetVector
  // (like a previous result)
  // ------------------------
  Events = new TChain(pars.ChainName);
  Events->Add(pars.InputName);
  if (NEvents < 0)
    NEvents = INT_MAX;

  // For picoDSTs
  // -------------
  if (pars.intype == INPICO || pars.intype == MCPICO)
  {
    pReader = SetupReader(Events, pars);

    InitializeReader(pReader, pars.InputName, NEvents, PicoDebugLevel, pars.HadronicCorr);
    if (pars.intype == MCPICO)
    {
      TurnOffCuts(pReader);
    }

    cout << "Don't forget to set tower cuts for pReader!!" << endl;
  }

  cout << "N = " << NEvents << endl;

  cout << "Done initializing chains. " << endl;
  return true;
}
//----------------------------------------------------------------------
// Main routine for one event.
EVENTRESULT ppTestAnalysis::RunEvent()
{
  TStarJetVector *sv;

  TStarJetPicoEventHeader *header = 0;

  UInt_t filehash = 0;
  TString cname = "";

  // Reset results (from last event)
  // -------------------------------
  Result.clear();
  weight = 1;
  particles.clear();
  // cout << " test " << endl;

  switch (pars.intype)
  {
    // =====================================================
  case INPICO:
  case MCPICO:
    if (!pReader->NextEvent())
    {
      pReader->PrintStatus();
      return EVENTRESULT::ENDOFINPUT;
      break;
    }
    pReader->PrintStatus(10);

    pFullEvent = pReader->GetOutputContainer()->GetArray(); // https://github.com/wsu-yale-rhig/TStarJetPicoReader/blob/master/TStarJetPicoReader/TStarJetPicoReaderBase.h#LL24C1-L24C1

    header = pReader->GetEvent()->GetHeader();

    refmult = header->GetProperReferenceMultiplicity();
    eventid1 = header->GetEventId();
    runid1 = header->GetRunId();
    vz = header->GetPrimaryVertexZ();
    pttot = 0;
    
    // For GEANT: Need to devise a runid that's unique but also
    // reproducible to match Geant and GeantMc data.
    if (pars.UseGeantNumbering)
    {
      TString cname = gSystem->BaseName(pReader->GetInputChain()->GetCurrentFile()->GetName());
      UInt_t filehash = cname.Hash();
      while (filehash > INT_MAX - 100000)
        filehash -= INT_MAX / 4; // some random large number
      if (filehash < 1000000)
        filehash += 1000001;
      runid = filehash;
      eventid = pReader->GetNOfCurrentEvent();
    }

    break;
    // =====================================================
  case INTREE:
  case MCTREE:
    if (evi >= NEvents)
    {
      return EVENTRESULT::ENDOFINPUT;
      break;
    }

    if (!(evi % 200))
      cout << "Working on " << evi << " / " << NEvents << endl;
    Events->GetEntry(evi);

    cname = Events->GetCurrentFile()->GetName();
    eventid = Events->GetLeaf("eventid")->GetValue();
    runid = Events->GetLeaf("runid")->GetValue();

    if (pars.intype == MCTREE)
    {
      filehash = cname.Hash();
      while (filehash > INT_MAX - 100000)
        filehash /= 10;
      if (filehash < 1000000)
        filehash += 1000001;
      runid += filehash;
    }

    ++evi;
    break;
    // =====================================================
  default:
    cerr << "Unknown intype " << pars.intype << endl;
    return EVENTRESULT::PROBLEM;
  }

  // Fill particle container
  // -----------------------
  for (int i = 0; i < pFullEvent->GetEntries(); ++i)
  {
    sv = (TStarJetVector *)pFullEvent->At(i);

    // Ensure kinematic similarity
    if (sv->Pt() < pars.PtConsMin || sv->Pt() > pars.PtConsMax)
      continue;
    if (fabs(sv->Eta()) > pars.EtaConsCut)
      continue;

    //! setting pid mass for particles from pythia
    if (pars.intype == MCPICO)
    {
      TParticlePDG *pid = (TParticlePDG *)PDGdb.GetParticle((Int_t)sv->mc_pdg_pid());
      double sv_mass = pid->Mass();
      double E = TMath::Sqrt(sv->P() * sv->P() + sv_mass * sv_mass);
      sv->SetPxPyPzE(sv->Px(), sv->Py(), sv->Pz(), E);
      pttot += sv->Pt();
    }

    //! setting pion mass for charged particles and zero mass for towers
    if (pars.intype == INPICO)
    {
      TParticlePDG *pid = (TParticlePDG *)PDGdb.GetParticle(211);
      double sv_mass;
      if (sv->GetCharge() != 0)
        sv_mass = pid->Mass(); //! tracks get pion mass
      else
        sv_mass = 0.0; //! towers get photon mass ~ 0
      double E = TMath::Sqrt(sv->P() * sv->P() + sv_mass * sv_mass);
      sv->SetPxPyPzE(sv->Px(), sv->Py(), sv->Pz(), E);
      pttot += sv->Pt();
    }

    // TRACKS
    // ------
    if (sv->GetCharge() != 0)
    {
      // EFFICIENCY uncertainty
      // ----------------------
      Double_t mran = gRandom->Uniform(0, 1);
      if (mran > pars.FakeEff)
      {
        continue;
      }
      if (SigmaPt)
      {
        // cout << sv->Pt() << "  " << SigmaPt->Eval( sv->Pt() ) << endl;
        SmearPt->SetParameters(1, 0, SigmaPt->Eval(sv->Pt()));
        double ptscaler = 1 + SmearPt->GetRandom();
        (*sv) *= ptscaler;
      }
    }

    // TOWERS
    // ------
    // Shift gain
    if (!sv->GetCharge())
    {
      (*sv) *= pars.fTowScale; // for systematics
    }

    particles.push_back(PseudoJet(*sv));
    particles.back().set_user_info(new JetAnalysisUserInfo(3 * sv->GetCharge(), sv->mc_pdg_pid(), "", sv->GetTowerID()));
  }

  mult = particles.size();
  if (particles.size() == 0)
    return EVENTRESULT::NOCONSTS;

  // For pythia, use cross section as weight
  // ---------------------------------------
  if (TParameter<double> *sigmaGen = (TParameter<double> *)Events->GetCurrentFile()->Get("sigmaGen"))
  {
    weight = sigmaGen->GetVal();
  }

  if (pars.InputName.Contains("Clean") || pars.InputName.Contains("pt-hat"))
  {
    TString currentfile = pReader->GetInputChain()->GetCurrentFile()->GetName();
    weight = LookupRun12Xsec(currentfile);
    if (fabs(weight - 1) < 1e-4)
    {
      throw std::runtime_error("mcweight unchanged!");
    }
  }

  // Run analysis
  // ------------
  if (pJA)
  {
    delete pJA;
    pJA = 0;
  }
  pJA = new JetAnalyzer(particles, JetDef);

  JetAnalyzer &JA = *pJA;
  vector<PseudoJet> JAResult = sorted_by_pt(select_jet(JA.inclusive_jets()));

  if (JAResult.size() == 0)
  {
    return EVENTRESULT::NOJETS;
  }

  highw = 0;
  // record event with a jet pT > 2*pTHat upper bound
  // if ((pars.InputName.Contains("2_3") && JAResult[0].perp() > 6.0) || (pars.InputName.Contains("3_4") && JAResult[0].perp() > 8.0) || (pars.InputName.Contains("4_5") && JAResult[0].perp() > 10.0) || (pars.InputName.Contains("5_7") && JAResult[0].perp() > 14.0) || (pars.InputName.Contains("7_9") && JAResult[0].perp() > 18.0) || (pars.InputName.Contains("9_11") && JAResult[0].perp() > 22.0) || (pars.InputName.Contains("11_15") && JAResult[0].perp() > 30.0) || (pars.InputName.Contains("15_20") && JAResult[0].perp() > 40.0) || (pars.InputName.Contains("20_25") && JAResult[0].perp() > 50.0) || (pars.InputName.Contains("25_35") && JAResult[0].perp() > 70.0) || (pars.InputName.Contains("35_-1") && JAResult[0].perp() > 1000.0))
  if ((pars.InputName.Contains("hat23_") && JAResult[0].perp() > 6.0) || (pars.InputName.Contains("hat34_") && JAResult[0].perp() > 8.0) || (pars.InputName.Contains("hat45_") && JAResult[0].perp() > 10.0) || (pars.InputName.Contains("hat57_") && JAResult[0].perp() > 14.0) || (pars.InputName.Contains("hat79_") && JAResult[0].perp() > 18.0) || (pars.InputName.Contains("hat911_") && JAResult[0].perp() > 22.0) || (pars.InputName.Contains("hat1115_") && JAResult[0].perp() > 30.0) || (pars.InputName.Contains("hat1520_") && JAResult[0].perp() > 40.0) || (pars.InputName.Contains("hat2025_") && JAResult[0].perp() > 50.0) || (pars.InputName.Contains("hat2535_") && JAResult[0].perp() > 70.0) || (pars.InputName.Contains("hat3545_") && JAResult[0].perp() > 90.0) || (pars.InputName.Contains("hat4555_") && JAResult[0].perp() > 110.0) || (pars.InputName.Contains("hat55999_") && JAResult[0].perp() > 1000.0))
  {
    highw = 1;
  }

  int njets = JAResult.size();

  fastjet::contrib::SoftDrop sd(0, 0.1);
  vector<fastjet::PseudoJet> sdjets = sd(JAResult);

  for (unsigned ijet = 0; ijet < JAResult.size(); ijet++)
  {
    PseudoJet &CurrentJet = JAResult[ijet];
    PseudoJet &SDJet = sdjets[ijet];
    vector<PseudoJet> IncPart = sorted_by_pt(CurrentJet.constituents());
    PseudoJet NeutralPart = join(OnlyNeutral(CurrentJet.constituents()));
    PseudoJet ChargedPart = join(OnlyCharged(CurrentJet.constituents()));
    int nch = ChargedPart.constituents().size();

    // leading track signs
    int b = -9;
    double dr = -9.0;
    double epair = -9.0;
    double z = -9.0;
    int pid1 = -9999;
    int pid2 = -9999;
    double pt1 = -9;
    double pt2 = -9;
    double qlead = -9;
    double qsublead = -9;
    double y1 = -9;
    double y2 = -9;
    double phi1 = -9;
    double phi2 = -9;

    vector<PseudoJet> ChargedPartVec = sorted_by_pt(ChargedPart.constituents());

    if (ChargedPartVec.size() >= 2)
    {
      qlead = ChargedPartVec.at(0).user_info<JetAnalysisUserInfo>().GetQuarkCharge();
      qsublead = ChargedPartVec.at(1).user_info<JetAnalysisUserInfo>().GetQuarkCharge();
      double elead = ChargedPartVec.at(0).e();
      double esublead = ChargedPartVec.at(1).e();
      pid1 = ChargedPartVec.at(0).user_info<JetAnalysisUserInfo>().GetPID();
      pid2 = ChargedPartVec.at(1).user_info<JetAnalysisUserInfo>().GetPID();
      pt1 = ChargedPartVec.at(0).perp();
      pt2 = ChargedPartVec.at(1).perp();
      y1 = ChargedPartVec.at(0).rap();
      y2 = ChargedPartVec.at(1).rap();
      phi1 = ChargedPartVec.at(0).phi();
      phi2 = ChargedPartVec.at(1).phi();

      epair = elead + esublead;
      z = elead / (elead + esublead);
      dr = ChargedPartVec.at(0).delta_R(ChargedPartVec.at(1));
      if (qlead * qsublead > 0)
        b = 1;
      else if (qlead * qsublead < 0)
        b = -1;
      else
        b = 0;
    }

    double ptne = 0.0;
    double pttot = 0.0;

    for (PseudoJet &n : NeutralPart.constituents())
    {
      ptne += n.perp();
    }

    for (PseudoJet &part : CurrentJet.constituents())
    {
      pttot += part.perp();
    }

    // jet charges
    double numerator0 = 0.0;
    double numerator2 = 0.0;

    for (PseudoJet &c : ChargedPart.constituents())
    {
      double charge = c.user_info<JetAnalysisUserInfo>().GetQuarkCharge() / 3.0;
      numerator0 += charge;
      numerator2 += pow(c.perp(), 2) * charge;
    }

    double q0 = numerator0;
    double q2 = numerator2 / pow(CurrentJet.perp(), 2);
    double nef = ptne / pttot;

    if (pars.MaxJetNEF < 1.0 && nef > pars.MaxJetNEF)
      continue;

    // check if the event has high weight
    int reject = 0;
    // if ((pars.InputName.Contains("2_3") && JAResult[0].perp() > 6.0) || (pars.InputName.Contains("3_4") && JAResult[0].perp() > 8.0) || (pars.InputName.Contains("4_5") && JAResult[0].perp() > 10.0) || (pars.InputName.Contains("5_7") && JAResult[0].perp() > 14.0) || (pars.InputName.Contains("7_9") && JAResult[0].perp() > 18.0) || (pars.InputName.Contains("9_11") && JAResult[0].perp() > 22.0) || (pars.InputName.Contains("11_15") && JAResult[0].perp() > 30.0) || (pars.InputName.Contains("15_20") && JAResult[0].perp() > 40.0) || (pars.InputName.Contains("20_25") && JAResult[0].perp() > 50.0) || (pars.InputName.Contains("25_35") && JAResult[0].perp() > 70.0) || (pars.InputName.Contains("35_-1") && JAResult[0].perp() > 1000.0))
    if ((pars.InputName.Contains("hat23_") && JAResult[0].perp() > 6.0) || (pars.InputName.Contains("hat34_") && JAResult[0].perp() > 8.0) || (pars.InputName.Contains("hat45_") && JAResult[0].perp() > 10.0) || (pars.InputName.Contains("hat57_") && JAResult[0].perp() > 14.0) || (pars.InputName.Contains("hat79_") && JAResult[0].perp() > 18.0) || (pars.InputName.Contains("hat911_") && JAResult[0].perp() > 22.0) || (pars.InputName.Contains("hat1115_") && JAResult[0].perp() > 30.0) || (pars.InputName.Contains("hat1520_") && JAResult[0].perp() > 40.0) || (pars.InputName.Contains("hat2025_") && JAResult[0].perp() > 50.0) || (pars.InputName.Contains("hat2535_") && JAResult[0].perp() > 70.0) || (pars.InputName.Contains("hat3545_") && JAResult[0].perp() > 90.0) || (pars.InputName.Contains("hat4555_") && JAResult[0].perp() > 110.0) || (pars.InputName.Contains("hat55999_") && JAResult[0].perp() > 1000.0))
    {
      reject = 1;
    }
    if (abs(vz) > pars.VzCut)
    {
      reject = 2;
    }
    if (eventid1 <= 0)
    {
      reject = 3;
    }
    Result.push_back(ResultStruct(CurrentJet, SDJet, nch, b, dr, q0, q2, epair, z, reject, pid1, pid2, pt1, pt2, qlead, qsublead, y1, y2, phi1, phi2, nef));
  }
  // By default, sort for original jet pt
  sort(Result.begin(), Result.end(), ResultStruct::origptgreater);

  return EVENTRESULT::JETSFOUND;
}
//----------------------------------------------------------------------
void InitializeReader(std::shared_ptr<TStarJetPicoReader> pReader, const TString InputName, const Long64_t NEvents,
                      const int PicoDebugLevel, const double HadronicCorr)
{

  TStarJetPicoReader &reader = *pReader;

  if (HadronicCorr < 0)
  {
    reader.SetApplyFractionHadronicCorrection(kFALSE);
    reader.SetApplyMIPCorrection(kTRUE);
    reader.SetRejectTowerElectrons(kTRUE);
  }
  else
  {
    reader.SetApplyFractionHadronicCorrection(kTRUE);
    reader.SetFractionHadronicCorrection(HadronicCorr);
    reader.SetApplyMIPCorrection(kFALSE);
    reader.SetRejectTowerElectrons(kFALSE);
  }

  reader.Init(NEvents);
  TStarJetPicoDefinitions::SetDebugLevel(PicoDebugLevel);
}
//----------------------------------------------------------------------
// Helper to deal with repetitive stuff
shared_ptr<TStarJetPicoReader> SetupReader(TChain *chain, const ppTestParameters &pars)
{
  TStarJetPicoDefinitions::SetDebugLevel(0); // 10 for more output

  shared_ptr<TStarJetPicoReader> pReader = make_shared<TStarJetPicoReader>();
  TStarJetPicoReader &reader = *pReader;
  reader.SetInputChain(chain);

  // Event and track selection
  // -------------------------
  TStarJetPicoEventCuts *evCuts = reader.GetEventCuts();
  evCuts->SetTriggerSelection(pars.TriggerName); // All, MB, HT, pp, ppHT, ppJP
  // Additional cuts
  evCuts->SetVertexZCut(pars.VzCut);
  evCuts->SetRefMultCut(pars.RefMultCut);
  evCuts->SetVertexZDiffCut(pars.VzDiffCut);
  evCuts->SetMaxEventPtCut(pars.MaxEventPtCut);
  evCuts->SetMaxEventEtCut(pars.MaxEventEtCut);

  evCuts->SetMinEventEtCut(pars.MinEventEtCut);

  std::cout << "Using these event cuts:" << std::endl;
  std::cout << " Vz: " << evCuts->GetVertexZCut() << std::endl;
  std::cout << " Refmult: " << evCuts->GetRefMultCutMin() << " -- " << evCuts->GetRefMultCutMax() << std::endl;
  std::cout << " Delta Vz:  " << evCuts->GetVertexZDiffCut() << std::endl;
  std::cout << " MaxEventPt:  " << evCuts->GetMaxEventPtCut() << std::endl;
  std::cout << " MaxEventEt:  " << evCuts->GetMaxEventEtCut() << std::endl;

  // This method does NOT WORK for GEANT MC trees because everything is in the tracks...
  // Do it by hand later on, using pars.ManualHtCut;
  // Also doesn't work for general trees, but there it can't be fixed

  // Tracks cuts
  TStarJetPicoTrackCuts *trackCuts = reader.GetTrackCuts();
  trackCuts->SetDCACut(pars.DcaCut);
  trackCuts->SetMinNFitPointsCut(pars.NMinFit);
  trackCuts->SetFitOverMaxPointsCut(pars.FitOverMaxPointsCut);
  trackCuts->SetMaxPtCut(pars.MaxTrackPt);

  std::cout << "Using these track cuts:" << std::endl;
  std::cout << " dca : " << trackCuts->GetDCACut() << std::endl;
  std::cout << " nfit : " << trackCuts->GetMinNFitPointsCut() << std::endl;
  std::cout << " nfitratio : " << trackCuts->GetFitOverMaxPointsCut() << std::endl;
  std::cout << " maxpt : " << trackCuts->GetMaxPtCut() << std::endl;

  // Towers
  TStarJetPicoTowerCuts *towerCuts = reader.GetTowerCuts();
  towerCuts->SetMaxEtCut(pars.MaxEtCut);

  std::cout << "Using these tower cuts:" << std::endl;
  std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
  std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;

  // V0s: Turn off
  reader.SetProcessV0s(false);

  // cout << "towerCuts=" << towerCuts << endl;
  // cout << "trackCuts=" << trackCuts << endl;
  return pReader;
}
//----------------------------------------------------------------------
void TurnOffCuts(std::shared_ptr<TStarJetPicoReader> pReader)
{

  pReader->SetProcessTowers(false);
  TStarJetPicoEventCuts *evCuts = pReader->GetEventCuts();
  evCuts->SetTriggerSelection("All"); // All, MB, HT, pp, ppHT, ppJP
  evCuts->SetVertexZCut(30);
  evCuts->SetRefMultCut(0);
  evCuts->SetVertexZDiffCut(999999);

  evCuts->SetMaxEventPtCut(99999);
  evCuts->SetMaxEventEtCut(99999);
  // evCuts->SetMinEventPtCut (-1);
  evCuts->SetMinEventEtCut(-1);

  evCuts->SetPVRankingCutOff(); //  Use SetPVRankingCutOff() to turn off vertex ranking cut.  default is OFF

  // Tracks cuts
  TStarJetPicoTrackCuts *trackCuts = pReader->GetTrackCuts();
  trackCuts->SetDCACut(99999);
  trackCuts->SetMinNFitPointsCut(-1);
  trackCuts->SetFitOverMaxPointsCut(-1);
  trackCuts->SetMaxPtCut(99999);

  // Towers: should be no tower in MC. All (charged or neutral) are handled in track
  TStarJetPicoTowerCuts *towerCuts = pReader->GetTowerCuts();
  towerCuts->SetMaxEtCut(99999);

  cout << " TURNED OFF ALL CUTS" << endl;
}

//----------------------------------------------------------------------
double LookupRun12Xsec(TString filename)
{

  const int NUMBEROFPT = 13;
  // const char *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};
  const static float XSEC[NUMBEROFPT] = {9.0012, 1.46253, 0.354566, 0.151622, 0.0249062, 0.00584527, 0.00230158, 0.000342755, 4.57002e-05, 9.72535e-06, 4.69889e-07, 2.69202e-08, 1.43453e-09};
  const static float NUMBEROFEVENT[NUMBEROFPT] = {3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 3e6, 2e6, 2e6, 1e6, 1e6};

  const static vector<string> vptbins = {"pt-hat23_", "pt-hat34_", "pt-hat45_", "pt-hat57_", "pt-hat79_", "pt-hat911_", "pt-hat1115_", "pt-hat1520_", "pt-hat2025_", "pt-hat2535_", "pt-hat3545_", "pt-hat4555_", "pt-hat55999_"};
  for (int i = 0; i < vptbins.size(); ++i)
  {
    if (filename.Contains(vptbins.at(i).data()))
      return XSEC[i] / NUMBEROFEVENT[i];
  }

  throw std::runtime_error("Not a valid filename");
  return -1;
}

//----------------------------------------------------------------------
/*
    convenient output
*/
ostream &operator<<(ostream &ostr, const PseudoJet &jet)
{
  if (jet == 0)
  {
    ostr << " 0 ";
  }
  else
  {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi()
         << " ClusSeq = " << (jet.has_associated_cs() ? "yes" : "no");
  }
  return ostr;
}