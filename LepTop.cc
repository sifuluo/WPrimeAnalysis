#include "Utilities/Analyzer.cc"
#include "Utilities/JESTools.cc"
#include "Utilities/GenTools.cc"
#include "Utilities/ROOTMini.cc"

#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TEfficiency.h>
#include <TString.h>
#include <TTree.h>

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>

using namespace std;



void LepTop(int SampleType = 0, int irun = 1, int OptionCode = 0, int debug = -1) {
  cout << "start" << endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "LepTopHad";
  double dRToMatch = 0.4;
  cout << "Running " << savename << endl;

  a->SetupROOTMini();
  a->SetOutput(savepath, savename);
  a->DebugMode(debug);

  a->CDOut();
  a->AddPlot(new TH1F("LepTopMatch","Leptonic Top Matching; Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepTopMatchScaled","Leptonic Top Matching After Scaling;Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepTopdR","Leptonic Top DeltaR with Gen",50,0.,5.));
  a->AddPlot(new TH1F("LepTopdRScaled","Leptonic Top DeltaR with Gen After Scaling",50,0.,5.));
  a->AddPlot(new TH1F("LepbMatch","Leptonic b Matching; Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepbMatchScaled","Leptonic b Matching After Scaling;Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepbdR","Leptonic b DeltaR with Gen",50,0.,5.));
  a->AddPlot(new TH1F("LepbdRScaled","Leptonic b DeltaR with Gen After Scaling",50,0.,5.));

  a->AddPlot(new TH2F("NeutrinodRVsLepbdR","Neutrino dR Vs Lepb dR;LepbdR;Neutrino dR",50,0.,5.,50,0.,5.));
  a->AddPlot(new TH2F("NeutrinodRVsLepbdRScaled","Neutrino dR Vs Lepb dR After Scale;LepbdR;Neutrino dR",50,0.,5.,50,0.,5.));
  a->AddPlot(new TH2F("LepbPtRatioVsdR","LepbPtRatio Vs Lepb dR;LepbdR;Pt Ratio",50,0.,5.,80,0.,4.));
  a->AddPlot(new TH2F("LepbPtRatioVsdR_PtLimited","LepbPtRatio Vs Pt Limited Lepb dR;LepbdR;Pt Ratio",50,0.,5.,40,0.,2.));
  a->AddPlot(new TH2F("WbdPhiVsMatching","dPhi of Lep W and b Vs matching;Matched, Not Matched, Not Found, Matched Scaled, Not Matched Scaled, Not Found Scaled",6,0,6,40,0.,4.));
  a->AddPlot(new TH2F("PtopVsMatching","Ptop Vs Matching;Matched, Not Matched, Not Found, Matched Scaled, Not Matched Scaled, Not Found Scaled",6,0,6,201,-20.,0.1));

  a->AddPlot(new TH1F("LepTopRes","Leptonic Top Mass Using Matched Leptonic b",600,0.,300.));
  a->AddPlot(new TH1F("LepTopResScaled","Leptonic Top Mass Using Matched Leptonic b",600,0.,300.));
  a->AddPlot(new TH1F("NJetCandidate","Number of Jets within dR < 0.4 of Lepb",5,-0.5,4.5));
  a->AddPlot(new TH2F("GenLepBHasNoMatch","Property of GenLepb without a possible match; Pt ; Eta",500,0,500,100,0,10));
  a->AddPlot(new TH2F("GenLepBHasMatch","Property of GenLepb with a possible match; Pt ; Eta",500,0,500,100,0,10));
  a->AddPlot(new TH1F("WrongLepbdPhi","dPhi of Wrong Lepb Selected with Gen", 40,0.,4.));
  a->AddPlot(new TH1F("WrongLepbdR","dR of Wrong Lepb Selected with Gen", 60,0.,6.));
  a->AddPlot(new TH1F("WrongLepbPtRatio","Pt Ratio of Wrong Lepb Selected with Gen", 80,0.,4.));
  a->AddPlot(new TH2F("LepbEtaCompare","Right Eta Vs Wrong Eta; Right Eta; Wrong Eta",80,-4.,4.,80,-4.,4.));
  a->AddPlot(new TH1F("WrongLepbdPhiScaled","dPhi of scaled Wrong Lepb Selected with Gen", 40,0.,4.));
  a->AddPlot(new TH1F("WrongLepbdRScaled","dR of scaled Wrong Lepb Selected with Gen", 60,0.,6.));
  a->AddPlot(new TH1F("WrongLepbPtRatioScaled","Pt Ratio of scaled Wrong Lepb Selected with Gen", 80,0.,4.));
  a->AddPlot(new TH2F("LepbEtaCompareScaled","Right Eta Vs scaled Wrong Eta; Right Eta; Wrong Eta",80,-4.,4.,80,-4.,4.));

  a->AddPlot(new TH2F("TopMassVsPLepMatched","TopMassVsPLepMatched; PLep; Top Mass",101,-10,0.1, 600,0.,300.));
  a->AddPlot(new TH2F("TopMassVsPLepNoMatch","TopMassVsPLepNoMatch; PLep; Top Mass",101,-10,0.1, 600,0.,300.));
  a->AddPlot(new TH2F("TopMassVsPLepMatchedScaled","TopMassVsPLepMatchedScaled; PLep; Top Mass",101,-10,0.1, 600,0.,300.));
  a->AddPlot(new TH2F("TopMassVsPLepNoMatchScaled","TopMassVsPLepNoMatchScaled; PLep; Top Mass",101,-10,0.1, 600,0.,300.));

  a->AddPlot(new TH1F("MatchWithLepbMatch","Match Efficiency with Lepb matched; Total LF0 LF1 Hadb WPb",8,-0.5,7.5));
  a->AddPlot(new TH1F("MatchWithLepbNoMatch","Match Efficiency with Lepb unmatched; Total LF0 LF1 Hadb WPb",8,-0.5,7.5));
  a->AddPlot(new TH2F("TopMassVsPHad","TopMassVsPHad; PHad; Top Mass",101,-10,0.1, 600,0.,300.));
  a->AddPlot(new TH2F("WPMassVsP","Correct SampleType Type WPMassVsP; P; W\' Mass",101,-10,0.1, 1000,0,1000));
  a->AddPlot(new TH2F("WPMass2D","WPMass2D;FL;LL",100,0,1000,100,0,1000));

  // a->Tree_Init();
  cout <<"Start Loop" <<endl;
  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue;

    map<string, TH1F*> p1d = a->Plots1D;
    map<string, TH2F*> p2d = a->Plots2D;
    a->CountEvent("PassPreSelection");
    vector<TLorentzVector> Jets = a->LVJets;
    vector<bool> BTags = a->BTags;
    if (Jets.size() != BTags.size()) {
      cout << "Jets size is different from BTags size" <<endl;
      continue;
    }
    TLorentzVector Lepton = a->LVLeptons.at(0);
    TLorentzVector MET = a->LVMET;
    TLorentzVector GenLepT = a->LVGenLepT;
    TLorentzVector GenLepB = a->LVGenLepB;
    TLorentzVector GenNeu = a->LVGenNeu;
    a->RM->SetLep(Lepton,MET);
    a->RM->SetGenNeu(GenNeu);
    a->JT->SetfdPhi();
    double PLep = 0;
    int LepTopBIndex = 0;
    TLorentzVector LepTopB = TLorentzVector();
    TLorentzVector LepTop = TLorentzVector();
    TLorentzVector Neutrino = TLorentzVector();

    double PLepScaled = 0;
    TLorentzVector LepTopBScaled = TLorentzVector();
    TLorentzVector LepTopScaled = TLorentzVector();
    TLorentzVector NeutrinoScaled = TLorentzVector();

    TLorentzVector MatchGenb = TLorentzVector();
    int TrueMatchb = 0;
    TLorentzVector MatchGenbLimited = TLorentzVector();
    int TrueMatchbLimited = 0;
    double MatchGenbdR = 6.;
    double MatchGenbdRLimited = 6.;

    int NJetCandidate = 0;

    for (unsigned ij = 0; ij < Jets.size(); ++ij) {
      double PBtag;
      if (BTags.at(ij)) PBtag = 1;
      else PBtag = 0.3 / 0.7;
      TLorentzVector Jet = Jets.at(ij);
      TLorentzVector Neutrinotemp = TLorentzVector();
      double JetPLep = a->JT->CalcPLep(Jet, Lepton, MET, GenNeu, Neutrinotemp);
      JetPLep *= PBtag;
      if (JetPLep > PLep) {
        PLep = JetPLep;
        LepTopBIndex = ij;
        LepTopB = Jet;
        LepTop = Jet + Lepton + Neutrinotemp;
        Neutrino = Neutrinotemp;
      }
      Neutrinotemp = TLorentzVector();

      double LepBScale = 0;
      double JetPLepScaled = a->RM->MinimizePLep(Jet,LepBScale);
      JetPLepScaled *= PBtag;
      if (JetPLepScaled > PLepScaled) {
        PLepScaled = JetPLepScaled;
        LepTopBScaled = Jet * LepBScale;
        TLorentzVector ScaledMET = MET + Jet - LepTopBScaled;
        a->JT->CalcPLep(LepTopBScaled, Lepton, ScaledMET, GenNeu, Neutrinotemp);
        LepTopScaled = LepTopBScaled + Lepton + Neutrinotemp;
        NeutrinoScaled = Neutrinotemp;
      }

      double GendR = Jet.DeltaR(GenLepB);
      double GenPtRatio =  Jet.Pt() / GenLepB.Pt();
      bool PtInRange = (GenPtRatio <= 1.5 && GenPtRatio >= 0.5);
      if (GendR < MatchGenbdR) {
        MatchGenbdR = GendR;
        MatchGenb = Jet;
      }
      if (GendR < MatchGenbdRLimited && PtInRange) {
        MatchGenbdRLimited = GendR;
        MatchGenbLimited = Jet;
      }
      if (GendR < dRToMatch && PtInRange) NJetCandidate++;
    }

    if (MatchGenbdR < dRToMatch)  {
      TrueMatchb = 1;
      p2d["LepbPtRatioVsdR"]->Fill(MatchGenbdR, MatchGenb.Pt() / GenLepB.Pt());
    }
    if (MatchGenbdRLimited < dRToMatch) {
      TrueMatchbLimited = 1;
      p2d["LepbPtRatioVsdR_PtLimited"]->Fill(MatchGenbdRLimited, MatchGenbLimited.Pt() / GenLepB.Pt());
    }

    if (PLep > 0) {
      double tdr = LepTop.DeltaR(GenLepT);
      // if (tdr < dRToMatch) p1d["LepTopMatch"]->Fill(0);
      // else p1d["LepTopMatch"]->Fill(1);
      p1d["LepTopMatch"]->Fill(tdr < dRToMatch ? 0 : 1);
      double bdr = LepTopB.DeltaR(GenLepB);
      int bmatch = (bdr < dRToMatch ? 0 : 1); // 0 for matched, 1 for not matched
      p1d["LepbMatch"]->Fill(bmatch);
      p1d["LepTopdR"]->Fill(tdr);
      p1d["LepbdR"]->Fill(bdr);
      double neudr = Neutrino.DeltaR(GenNeu);
      p2d["NeutrinodRVsLepbdR"]->Fill(bdr, neudr);
      double wbdphi = fabs((Lepton + Neutrino).DeltaPhi(LepTopB));
      p2d["WbdPhiVsMatching"]->Fill(bmatch,wbdphi);
      p2d["PtopVsMatching"]->Fill(bmatch, log10(PLep));
      p1d["LepTopRes"]->Fill( (MatchGenb + Neutrino + Lepton).M() );
      if (bmatch) {
        p1d["WrongLepbdPhi"]->Fill(fabs(LepTopB.DeltaPhi(GenLepB)));
        p1d["WrongLepbdR"]->Fill(bdr);
        p1d["WrongLepbPtRatio"]->Fill(LepTopB.Pt() / GenLepB.Pt());
        if (TrueMatchbLimited) p2d["LepbEtaCompare"]->Fill(MatchGenbLimited.Eta(),LepTopB.Eta());
        p2d["TopMassVsPLepNoMatch"]->Fill(log10(PLep),LepTop.M());
      }
      else p2d["TopMassVsPLepMatched"]->Fill(log10(PLep),LepTop.M());
    }
    else p1d["LepTopMatch"]->Fill(2);

    if (PLepScaled > 0) {
      double tdr = LepTopScaled.DeltaR(GenLepT);
      if (tdr < dRToMatch) p1d["LepTopMatchScaled"]->Fill(0);
      else p1d["LepTopMatchScaled"]->Fill(1);
      double bdr = LepTopBScaled.DeltaR(GenLepB);
      int bmatch = (bdr < dRToMatch ? 0 : 1); // 0 for matched, 1 for not matched
      p1d["LepbMatchScaled"]->Fill(bmatch);
      p1d["LepTopdRScaled"]->Fill(tdr);
      p1d["LepbdRScaled"]->Fill(bdr);
      double neudr = NeutrinoScaled.DeltaR(GenNeu);
      p2d["NeutrinodRVsLepbdRScaled"]->Fill(bdr, neudr);
      double wbdphi = fabs((Lepton + NeutrinoScaled).DeltaPhi(LepTopBScaled));
      p2d["WbdPhiVsMatching"]->Fill(bmatch + 3, wbdphi);
      p2d["PtopVsMatching"]->Fill(bmatch+3,log10(PLepScaled));

      if (bmatch) {
        p1d["WrongLepbdPhiScaled"]->Fill(fabs(LepTopBScaled.DeltaPhi(GenLepB)));
        p1d["WrongLepbdRScaled"]->Fill(bdr);
        p1d["WrongLepbPtRatioScaled"]->Fill(LepTopBScaled.Pt() / GenLepB.Pt());
        if (TrueMatchbLimited) p2d["LepbEtaCompareScaled"]->Fill(MatchGenbLimited.Eta(),LepTopBScaled.Eta());
        p2d["TopMassVsPLepNoMatchScaled"]->Fill(log10(PLepScaled),LepTopScaled.M());
      }
      else p2d["TopMassVsPLepMatchedScaled"]->Fill(log10(PLepScaled),LepTopScaled.M());
    }
    else p1d["LepTopMatchScaled"]->Fill(2);

    p1d["NJetCandidate"]->Fill(NJetCandidate);
    (NJetCandidate ? p2d["GenLepBHasMatch"] : p2d["GenLepBHasNoMatch"])->Fill(GenLepB.Pt(),GenLepB.Eta());

    if (PLep <= 0) continue;

    // Starting the Hadronic Part
    bool LepbMatch = (LepTopB.DeltaR(GenLepB) < dRToMatch);
    int jetsize = Jets.size();
    vector<vector<int> > Permutations;
    for (int lf0 = 0;lf0 < jetsize - 1; ++lf0) {
      if (lf0 == LepTopBIndex) continue;
      for (int lf1 = lf0 + 1; lf1 < jetsize; ++lf1) {
        if (lf1 == LepTopBIndex) continue;
        for (int hadb = 0; hadb < jetsize; ++hadb) {
          if (hadb == lf0 || hadb == lf1 || hadb == LepTopBIndex) continue;
          vector<int> perm{lf0,lf1,hadb};
          Permutations.push_back(perm);
        }
      }
    }

    double PHad = 0;
    vector<int> HadPerm;
    vector<TLorentzVector> HadJets;
    for (unsigned iperm = 0; iperm < Permutations.size(); ++iperm) {
      double phad = 1;
      vector<int> perm = Permutations.at(iperm);
      double PBtag = 1;
      if (BTags.at(perm.at(0))) PBtag*= 0.01 / 0.99;
      if (BTags.at(perm.at(1))) PBtag*= 0.01 / 0.99;
      if (!BTags.at(perm.at(2))) PBtag*= 0.3 / 0.7;
      vector<TLorentzVector> permjets;
      for (unsigned ijet = 0; ijet < perm.size(); ++ijet) permjets.push_back(Jets.at(perm.at(ijet)));
      phad = a->JT->CalcPHad(permjets) * PBtag;
      if (phad > PHad) {
        PHad = phad;
        HadPerm = perm;
        HadJets = permjets;
      }
    }
    int LFMatch = 0;
    int HadbMatch = 0;
    TLorentzVector HadTop = HadJets.at(0) + HadJets.at(1) + HadJets.at(2);
    if (HadJets.at(0).DeltaR(a->LVGenLF0) < dRToMatch || HadJets.at(0).DeltaR(a->LVGenLF1) < dRToMatch) LFMatch++;
    if (HadJets.at(1).DeltaR(a->LVGenLF0) < dRToMatch || HadJets.at(1).DeltaR(a->LVGenLF1) < dRToMatch) LFMatch++;
    if (HadJets.at(2).DeltaR(a->LVGenHadB) < dRToMatch) HadbMatch = 1;
    TH1F* HadMatchEff = (LepbMatch ? p1d["MatchWithLepbMatch"] : p1d["MatchWithLepbNoMatch"]);
    HadMatchEff->Fill(0);
    if (LFMatch > 0) HadMatchEff->Fill(1);
    if (LFMatch > 1) HadMatchEff->Fill(2);
    if (HadbMatch > 0) HadMatchEff->Fill(3);
    p2d["TopMassVsPHad"]->Fill(PHad,HadTop.M());

    // Staring WPb selection
    double WPbPt = 0;
    TLorentzVector WPb = TLorentzVector();
    for (int ib = 0; ib < Jets.size(); ++ib) {
      if (ib == LepTopBIndex || ib == HadPerm.at(0) || ib == HadPerm.at(1) || ib == HadPerm.at(2) ) continue;
      TLorentzVector wpb = Jets.at(ib);
      if (wpb.Pt() > WPbPt) {
        WPb = wpb;
      }
    }
    if (WPb.DeltaR(a->LVGenWPB) < dRToMatch) HadMatchEff->Fill(4);
    double PTotal = PLep * PHad;
    double FLTypeWPMass = (HadTop + WPb).M();
    double LLTypeWPMass = (LepTop + WPb).M();
    if (SampleType == 0) p2d["WPMassVsP"]->Fill(PTotal, FLTypeWPMass);
    if (SampleType == 1) p2d["WPMassVsP"]->Fill(PTotal, LLTypeWPMass);
    if (SampleType == 2) {
      double WPbdRHadT = WPb.DeltaR(HadTop);
      double WPbdRLepT = WPb.DeltaR(LepTop);
      if (WPbdRHadT > WPbdRLepT) p2d["WPMassVsP"]->Fill(PTotal, FLTypeWPMass);
      else p2d["WPMassVsP"]->Fill(PTotal, LLTypeWPMass);
    }
    p2d["WPMass2D"]->Fill(FLTypeWPMass, LLTypeWPMass);

  }

  a->SaveOutput();
}
