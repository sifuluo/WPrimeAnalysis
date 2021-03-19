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

void FitTTBar(int SampleType = 0, int irun = 1, int OptionCode = 0, int debug = -1) {
  cout << "start" << endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "FitTTBarWPB";
  double dRToMatch = 0.4;
  cout << "Running " << savename << endl;
  a->SetupROOTMini();
  a->SetOutput(savepath, savename);
  a->DebugMode(debug);
  a->CDOut();
  TH1F* METPt = new TH1F("METPt","METPt",400,0,400);
  TH2F* WPMass = new TH2F("WPMass","WPMass; FL-like ; LL-like",1000,0,1000,1000,0,1000);
  TH2F* WPMassMatched = new TH2F("WPMassMatched","WPMassMatched; FL-like ; LL-like",1000,0,1000,1000,0,1000);
  TH1F* WPBPt = new TH1F("WPBPt","WPBPt",500,0,500);
  TH1F* OtherRemainPt = new TH1F("OtherRemainPt","OtherRemainPt",500,0,500);
  TH1F* WPBdPhi = new TH1F("WPBdPhi","WPBdPhi",40,0,4);
  TH1F* OtherRemaindPhi = new TH1F("OtherRemaindPhi","OtherRemaindPhi",40,0,4);
  a->Tree_Init();
  // TTree* t = new TTree("t","EventTree");
  vector<double>* PScales = new vector<double>; // LF0 LF1 HadB LepB
  vector<double>* PPreMass = new vector<double>; // HadW HadT LepW LepT
  vector<double>* PPostMass = new vector<double>; // HadW HadT LepW LepT
  vector<double>* PBTags = new vector<double>;
  // vector<double>* POthers = new vector<double>;
  vector<double>* Scales = new vector<double>;

  vector<double>* PScales_Match = new vector<double>; // LF0 LF1 HadB LepB
  vector<double>* PPreMass_Match = new vector<double>; // HadW HadT LepW LepT
  vector<double>* PPostMass_Match = new vector<double>; // HadW HadT LepW LepT
  vector<double>* PBTags_Match = new vector<double>;
  // vector<double>* POthers_Match = new vector<double>;
  vector<double>* Scales_Match = new vector<double>;
  Hypothesis Pre, Post, Match;

  a->t->Branch("PScales",&PScales);
  a->t->Branch("PPreMass",&PPreMass);
  a->t->Branch("PPostMass",&PPostMass);
  a->t->Branch("PBTags",&PBTags);
  // a->t->Branch("POthers",&POthers);
  a->t->Branch("Scales",&Scales);

  a->t->Branch("PScales_Match",&PScales_Match);
  a->t->Branch("PPreMass_Match",&PPreMass_Match);
  a->t->Branch("PPostMass_Match",&PPostMass_Match);
  a->t->Branch("PBTags_Match",&PBTags_Match);
  // a->t->Branch("POthers_Match",&POthers_Match);
  a->t->Branch("Scales_Match",&Scales_Match);
  Pre.BookBranches(a->t, "Pre", true);
  Post.BookBranches(a->t, "Post", true);
  Match.BookBranches(a->t, "Match", true);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue;
    a->CountEvent("PassPreSelection");
    if (a->nBJets < 2) continue;
    a->CountEvent("Pass2BJets");
    vector<int> TruePerm = a->Tree_Reco();
    a->Tree_FitReco();
    bool AllFilled = a->Reco.AllFilled();
    PScales->clear();
    PPreMass->clear();
    PPostMass->clear();
    PBTags->clear();
    Scales->clear();

    PScales_Match->clear();
    PPreMass_Match->clear();
    PPostMass_Match->clear();
    PBTags_Match->clear();
    Scales_Match->clear();

    vector<TLorentzVector> Jets = a->LVJets;
    TLorentzVector hadb;
    vector<bool> BTags = a->BTags;
    a->RM->SetLep(a->LVLeptons[0],a->LVMET);
    METPt->Fill(a->LVMET.Pt());
    vector<vector<int> > perms = a->JT->MakePermutations5(Jets.size());
    double BestP = 0;
    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      bool IsTruePerm;
      vector<int> thisperm = perms.at(iperm);
      if (thisperm == TruePerm && AllFilled) IsTruePerm = true;
      double thispbtags = a->JT->CalcPFlavor(thisperm, BTags);
      vector<TLorentzVector> LVPerm = a->JT->GetPermutationLV(thisperm, Jets);
      vector<TLorentzVector> LVPerm4 = LVPerm;
      LVPerm4.resize(4);
      if (find(LVPerm.begin(), LVPerm.end(), TLorentzVector()) != LVPerm.end()) continue;
      double PJes = a->RM->MinimizeP(LVPerm4);
      if (PJes < 0) continue;
      double PPerm = PJes * thispbtags;
      if (PPerm > BestP || IsTruePerm) {
        vector<vector<TLorentzVector> > LVSets;
        vector<vector<double> > PVector = a->RM->ReCalcPVector(LVSets);
        LVSets[0][4] = Jets.at(thisperm.at(4));
        LVSets[1][4] = Jets.at(thisperm.at(4));
        if (PPerm > BestP) {
          BestP = PPerm;
          Pre.SetLV(LVSets[0]);
          Post.SetLV(LVSets[1]);
          Pre.Calculate(SampleType);
          Post.Calculate(SampleType);
          *Scales = PVector[0];
          *PScales = PVector[1];
          *PPreMass = PVector[2];
          *PPostMass = PVector[3];
          *PBTags = a->JT->CalcPFlavorVector(thisperm, BTags);
          //Consistency check
          double pbtagstemp = 1;
          for (unsigned ij = 0; ij < PBTags->size(); ++ij) pbtagstemp *= PBTags->at(ij);
          if (pbtagstemp != thispbtags) cout << "\n Product of vector = " << pbtagstemp << " ; while overall = " << thispbtags <<endl;

          for (unsigned ij = 0; ij < 4; ++ij) {
            if (Pre.Observables().at(ij) != LVPerm.at(ij)) cout << endl << "Pre-scale " << ij<< " th jet is inconsistent" <<endl;
          }
          hadb = LVPerm4[2];
        }
        if (IsTruePerm) {
          Match.SetLV(LVSets[1]);
          Match.Calculate(SampleType);
          *Scales_Match = PVector[0];
          *PScales_Match = PVector[1];
          *PPreMass_Match = PVector[2];
          *PPostMass_Match = PVector[3];
          *PBTags_Match = a->JT->CalcPFlavorVector(thisperm, BTags);
        }
        else  {
          Match.Reset();
          Scales_Match->clear();
          PScales_Match->clear();
          PPreMass_Match->clear();
          PPostMass_Match->clear();
          PBTags_Match->clear();
        }
      }
    }
    // cout <<endl;
    // (hadb - Pre.HadB).Print();
    // cout <<endl;
    // (Pre.HadB * Scales->at(2) - Post.HadB).Print();
    // cout <<endl;
    a->Tree_Fill();
    if (BestP <= 0) continue;
    WPMass->Fill(Post.FL_WPMass,Post.LL_WPMass);
    if (AllFilled) WPMassMatched->Fill(Match.FL_WPMass,Match.LL_WPMass);

    if (Jets.size() > 5 && SampleType != 2) {
      vector<TLorentzVector> RemainJets;
      for (TLorentzVector jet : Jets) {
        if (jet==Pre.LF0 || jet==Pre.LF1 || jet==Pre.HadB || jet==Pre.LepB) continue;
        RemainJets.push_back(jet);
      }
      if (RemainJets.size() > Jets.size() - 4) {
        cout << "How can you find so many remaining jets?" <<endl;
        cout << "Jets: " << Jets.size() << " ; RemainJets: " << RemainJets.size() << endl;
        // cout << "Jets List: " <<endl;
        // for (TLorentzVector jet: Jets) {
        //   jet.Print();
        //   cout << endl;
        // }
        // cout << endl << "Remains:" <<endl;
        // for (TLorentzVector jet: RemainJets) {
        //   jet.Print();
        //   cout <<endl;
      }
      TLorentzVector WPTop;
      if (SampleType == 0) WPTop = Post.HadT;
      if (SampleType == 1) WPTop = Post.LepT;
      TLorentzVector RecoWPB = RemainJets[0];
      for (auto jet : RemainJets) {
        if (jet.DeltaR(a->Gen.WPB) < RecoWPB.DeltaR(a->Gen.WPB)) RecoWPB = jet;
      }
      if (RecoWPB.DeltaR(a->Gen.WPB) < 0.4) {
        WPBPt->Fill(RecoWPB.Pt());
        WPBdPhi->Fill(RecoWPB.DeltaPhi(WPTop));
        for (auto jet : RemainJets) {
          if (jet == RecoWPB) continue;
          OtherRemainPt->Fill(jet.Pt());
          OtherRemaindPhi->Fill(jet.DeltaPhi(WPTop));
        }
      }
    }
  }

  a->Tree_Save();
  a->SaveOutput();

}
