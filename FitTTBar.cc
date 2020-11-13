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
  TString savename = "FitTTBar";
  double dRToMatch = 0.4;
  cout << "Running " << savename << endl;
  a->SetupROOTMini();
  a->SetOutput(savepath, savename);
  a->DebugMode(debug);
  a->CDOut();
  TH2F* WPMass = new TH2F("WPMass","WPMass; FL-like ; LL-like",1000,0,1000,1000,0,1000);
  TH2F* WPMassMatched = new TH2F("WPMassMatched","WPMassMatched; FL-like ; LL-like",1000,0,1000,1000,0,1000);
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
    vector<int> TruePerm = a->Tree_Reco();
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

    a->CountEvent("PassPreSelection");
    vector<TLorentzVector> Jets = a->LVJets;
    vector<bool> BTags = a->BTags;
    a->RM->SetLep(a->LVLeptons[0],a->LVMET);
    vector<vector<int> > perms = a->JT->MakePermutations5(Jets.size());
    double BestP = -1;
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
      }
    }
    a->Tree_Fill();
    if (BestP > 0) {
      WPMass->Fill(Post.FL_WPMass,Post.LL_WPMass);
      if (AllFilled) WPMassMatched->Fill(Match.FL_WPMass,Match.LL_WPMass);
    }
    // a->t->Fill();
  }

  a->Tree_Save();
  a->SaveOutput();

}
