#include "Utilities/Analyzer.cc"
#include "Utilities/JESTools.cc"
#include "Utilities/ROOTMini.cc"

// Delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include <TROOT.h>
#include <TH1.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>

void GenMinimizeAllPermTest(int SampleType = 0, int irun = 1, int debug = 0) {
  cout << "start" <<endl;
  //SetUpUtilities
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MinimizerAllPermTest";


  JESTools *b = new JESTools();
  TFile* PFile = new TFile("PFile/PFile.root");
  b->ReadJESPlots(PFile);

  ROOTMini *m = new ROOTMini(b);
  // m->SetMinimizer();

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);
  m->SetDebug(0);

  TH1F* BestPDis     = new TH1F("BestPDis"    ,"Best P distribution", 100,0,1);


  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;

    // Input set up
    vector<TLorentzVector> AllJets;
    // vector<int> GenOutQuark = a->GenOutQuark;
    vector<int> GenOutQuark = a->GenOutSort;
    bool allfound = true;
    for (unsigned ig = 0; ig < GenOutQuark.size(); ++ig) {
      if (GenOutQuark[ig] == -1) {
        allfound = false;
        AllJets.push_back(TLorentzVector());
      }
      else {
        AllJets.push_back(a->GetGenParticleP4(GenOutQuark[ig]));
      }
    }
    if (!allfound) {
      cout << "Set of Gen particles incomplete, skipping" <<endl;
      continue;
    }

    TLorentzVector Lepton = a->LVGenLep;
    TLorentzVector LVMET = a->LVGenNeu;
    m->SetLep(Lepton, LVMET);

    vector<double> PermPs;
    vector< vector<int> > perms = b->MakePermutations(AllJets.size());
    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      vector<int> thisperm = perms.at(iperm);
      vector<TLorentzVector> Jets;
      for (unsigned ip = 0; ip < thisperm.size(); ++ip) {
        Jets.push_back(AllJets.at(thisperm.at(ip)));
      }
      double ThisP = m->MinimizeP(Jets);

    }

    BestPDis->Fill(BestP);
    // cout << endl << "BestP = " << BestP << Form(" BestPerm is %i,%i,%i,%i",BestPerm[0],BestPerm[1],BestPerm[2],BestPerm[3]) <<endl;
    // if (!RightPerm) {
    //   WrongPDis->Fill(BestP);
    //   RightPDis->Fill(RightPermP);
    //   cout << "Wrong but Best P = " << BestP << ", Right Perm P = " << RightPermP <<endl;
    // }

  }
  a->SaveOutput();
}
