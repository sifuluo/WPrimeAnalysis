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

void MatchAll(int SampleType = 0, int irun = 1, int OptionCode = 0, int debug = -1) {
  cout << "start" << endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MatchAll";
  double dRToMatch = 0.4;
  cout << "Running " << savename << endl;

  a->SetupROOTMini();
  a->SetOutput(savepath, savename);
  a->DebugMode(debug);
  a->CDOut();
  // Histograms

  vector<int>* Conflict = new vector<int>;
  vector<double>* DeltaR = new vector<double>;
  vector<double>* PtRatio = new vector<double>;
  a->Tree_Init(); // Init tree and initialize gen and reco branch
  a->t->Branch("Conflict",&Conflict);
  a->t->Branch("DeltaR",&DeltaR);
  a->t->Branch("PtRatio",&PtRatio);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue;
    a->Tree_Reco();
    Conflict->clear();
    DeltaR->clear();
    PtRatio->clear();
    a->CountEvent("PassPreSelection");
    vector<TLorentzVector> Jets = a->LVJets;
    vector<TLorentzVector> GenJets = a->LVGenOutSort; // LF0, LF1, HadB, LepB, W' B
    vector<int> tmpvec;
    for (unsigned igen = 0; igen < GenJets.size(); ++igen) {
      map<double,int> matchcontainer; // Sort by deltaR before put in vector
      TLorentzVector lvgen = GenJets.at(igen);
      if (lvgen.Pt() == 0) {
        tmpvec.push_back(-2);
        continue;
      }
      for (unsigned ij = 0; ij < Jets.size(); ++ij) {
        TLorentzVector lvreco = Jets.at(ij);
        double dr = lvreco.DeltaR(lvgen);
        if (dr < dRToMatch) {
          matchcontainer[dr] = ij;
        }
      }
      for (const auto &it: matchcontainer) {
        tmpvec.push_back(it.second);
      }
      tmpvec.push_back(-2);
    }

    int indexgen = 0;
    int indexconf = 1;
    Conflict->resize(tmpvec.size(),0);
    for (unsigned iv = 0; iv < tmpvec.size(); ++iv) {
      bool has_conf = false;
      int indexreco = tmpvec[iv];
      if (indexreco == -2) {
        Conflict->at(iv) = -2;
        DeltaR->push_back(-2);
        PtRatio->push_back(-2);
        indexgen++;
        continue;
      }
      for (unsigned iv2 = iv + 1; iv2 < tmpvec.size(); ++iv2) {
        int indexreco2 = tmpvec[iv2];
        if (indexreco2 == indexreco) {
          has_conf = true;
          Conflict->at(iv2) = indexconf;
        }
      }
      if (has_conf) {
        Conflict->at(iv) = indexconf;
        indexconf++;
      }
      TLorentzVector jet = Jets.at(indexreco);
      TLorentzVector gen = GenJets.at(indexgen);
      DeltaR->push_back(jet.DeltaR(gen));
      PtRatio->push_back(jet.Pt() / gen.Pt());
    }
    // cout << Form("vector sizes: tmpvec = %d, Conflict = %d, DeltaR = %d",tmpvec.size(),Conflict->size(), DeltaR->size())<<endl;
    a->Tree_Fill();
  }

  a->Tree_Save();
  a->SaveOutput();



}
