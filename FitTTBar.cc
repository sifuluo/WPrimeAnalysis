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
  // TTree* t = new TTree("t","EventTree");
  vector<double>* PScales = new vector<double>; // LF0 LF1 HadB LepB
  vector<double>* PPreMass = new vector<double>; // HadW HadT LepW LepT
  vector<double>* PPostMass = new vector<double>; // HadW HadT LepW LepT
  vector<double>* PBTags = new vector<double>;
  // vector<double>* POthers = new vector<double>;
  vector<double>* Scales = new vector<double>;
  Hypothesis Pre, Post;
  a->t->Branch("PScales",&PScales);
  a->t->Branch("PPreMass",&PPreMass);
  a->t->Branch("PPostMass",&PPostMass);
  a->t->Branch("PBTags",&PBTags);
  // a->t->Branch("POthers",&POthers);
  a->t->Branch("Scales",&Scales);
  Pre.BookBranches(t, "Pre", false);
  Post.BookBranches(t, "Post", false);

  a->Tree_Init();

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue;
    a->Tree_Reco();
    PScales->clear();
    PPreMass->clear();
    PPostMass->clear();
    PBTags->clear();
    Scales->clear();
    a->CountEvent("PassPreSelection");
    vector<TLorentzVector> Jets = a->LVJets;
    vector<bool> BTags = a->BTags;
    a->RM->SetLep(a->LVLeptons[0],a->LVMET);
    vector<vector<int> > perms = a->JT->MakePermutations5(Jets.size());
    double BestP = -1;
    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      vector<int> thisperm = perms.at(iperm);
      double thispbtags = a->JT->CalcPFlavor(thisperm, BTags);
      vector<TLorentzVector> LVPerm = a->JT->GetPermutationLV(thisperm, Jets);
      vector<TLorentzVector> LVPerm4 = LVPerm;
      LVPerm4.resize(4);
      if (find(LVPerm.begin(), LVPerm.end(), TLorentzVector()) != LVPerm.end()) continue;
      double PJes = a->RM->MinimizeP(LVPerm4);
      if (PJes < 0) continue;
      double PPerm = PJes * thispbtags;
      if (PPerm > BestP) {
        BestP = PPerm;
        vector<vector<TLorentzVector> > LVSets;
        vector<vector<double> > PVector = a->RM->ReCalcPVector(LVSets);
        LVSets[0][4] = Jets.at(thisperm.at(4));
        LVSets[1][4] = Jets.at(thisperm.at(4));
        Pre.SetLV(LVSets[0]);
        Post.SetLV(LVSets[1]);
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
    }


    a->Tree_Fill();
    // a->t->Fill();
  }

  a->Tree_Save();
  a->SaveOutput();

}
