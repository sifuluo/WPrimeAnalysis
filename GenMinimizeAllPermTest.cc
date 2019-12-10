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
  TString savename = "MinimizerAllPerm4Test";


  JESTools *b = new JESTools();
  TFile* PFile = new TFile("PFile/PFile.root");
  b->ReadJESPlots(PFile);

  ROOTMini *m = new ROOTMini(b);
  // m->SetMinimizer();

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);
  m->SetDebug(0);

  vector<double> PermPs;
  vector< vector<int> > perms = b->MakePermutations(4);
  int nperm = perms.size();
  cout << endl;
  for (int ip = 0; ip < nperm; ++ip) {
    vector<int> thisperm = perms[ip];
    cout << Form("Perm: %i: (%i,%i,%i,%i)", ip, thisperm[0],thisperm[1],thisperm[2],thisperm[3]) <<endl;
  }
  TH1F* RightPermP = new TH1F("RightPermP", "P distribution of the correct perm",100,0,1);
  TH1F* BestPermP     = new TH1F("BestPermP"    ,"Best Perm P distribution", 100,0,1);
  TH1F* PermChoice = new TH1F("PermChoice","Permutation with the best P",nperm+1, -0.5,nperm+0.5);
  TH2F* BestPVsPerm  = new TH2F("BestPVsPerm", "Best P of the Best Permutation;PermIndex;P",nperm+1,-0.5,nperm+0.5,100,0.,1.);
  TH2F* BestPDiffVsPerm  = new TH2F("BestPDiffVsPerm", "Diff Best P of the Best Permutation and the right one;PermIndex;P",nperm+1,-0.5,nperm+0.5,100,0.,1.);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;

    // Input set up
    vector<TLorentzVector> AllJets;
    // vector<int> GenOutQuark = a->GenOutQuark;
    vector<int> GenOutQuark = a->GenOutSort;
    bool allfound = true;
    for (unsigned ig = 0; ig < 4; ++ig) {
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

    double RightP = 0;
    double BestP = 0;
    int BestPerm = 0;
    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      vector<int> thisperm = perms.at(iperm);
      vector<TLorentzVector> Jets;
      for (unsigned ip = 0; ip < thisperm.size(); ++ip) {
        Jets.push_back(AllJets.at(thisperm.at(ip)));
      }
      double ThisP = m->MinimizeP(Jets);
      if (iperm == 0) RightP = ThisP;
      if (ThisP > BestP) {
        BestP = ThisP;
        BestPerm = iperm;
      }
    }
    RightPermP->Fill(RightP);
    BestPermP->Fill(BestP);
    PermChoice->Fill(BestPerm);
    BestPVsPerm->Fill(BestPerm,BestP);
    if (BestPerm != 0) {
      BestPDiffVsPerm->Fill(BestPerm,BestP-RightP);
    }


  }
  a->SaveOutput();
}
