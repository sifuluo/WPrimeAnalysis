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

void MinimizerTestOutput(int SampleType = 0, int irun = 1, int testmode = 0) {
  cout << "Test Starts" <<endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MinimizerTestOutput";

  JESTool *b = new JESTool();
  TFile* PFile = new TFile("PFile/PFile.root");
  b->ReadJESPlots(PFile);
  ROOTMini *m = new ROOTMini(b);

  a->SetOutput(savepath,savename);
  a->DebugMode(testmode);
  m->SetDebug(1);
  m->SetTempOutput(1);

  vector<double> PermPs;
  vector< vector<int> > perms = b->MakePermutations(4);
  int nperm = perms.size();
  cout << endl;
  for (int ip = 0; ip < nperm; ++ip) {
    vector<int> thisperm = perms[ip];
    cout << Form("Perm: %i: (%i,%i,%i,%i)", ip, thisperm[0],thisperm[1],thisperm[2],thisperm[3]) <<endl;
  }

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    vector<TLorentzVector> AllJets;
    vector<int> GenOutQuark = a->GenOutSort;


  }


}
