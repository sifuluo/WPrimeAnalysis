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

void RecoTree(int SampleType = 0, int irun = 0, int OptionCode = 0, int debug = -2) {
  cout << "start" << endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "RecoTree";
  cout << "Running " << savename << endl;

  a->SetupROOTMini();
  a->SetOutput(savepath, savename);
  a->DebugMode(debug);

  a->CDOut();

  // a->AddPlot(new TH1F("GenOutSortMissing","GenOutSortMissings",5,-0.5,4.5));
  a->AddPlot(new TH1F("GenOutdR","GenOutdR",1010,-0.01,1.));
  a->AddPlot(new TH1F("GenOutPTRatio","GenOutdR",10000,0.5,1.5));
  a->AddPlot(new TH1F("GenOutDiffMag","GenOutDiffMag",1000,0.,10.));

  // a->Tree_Init();
  cout <<"Start Loop" <<endl;
  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    // if (a->RecoPass == -1) continue;
    map<string, TH1F*> p1d = a->Plots1D;
    map<string, TH2F*> p2d = a->Plots2D;

    // Verify quark missing in Gen Level.
    for (unsigned i = 0; i < 5; ++i) {
      if (a->GenOutSort[i] == -1) cout << endl<<i<<"th GenOutSort is -1" <<endl;
    }

    // From Gen to GenOutPart
    for (unsigned i = 0; i <5; ++i) {
      TLorentzVector lvgen = a->GenParticles[a->GenOutSort[i]]->P4();
      TLorentzVector lvgenout = a->LVGenOutSort[i];
      double dr = lvgen.DeltaR(lvgenout);
      if (dr >1.) dr = -0.005;
      p1d["GenOutdR"]->Fill(dr);

      double ratio = lvgenout.Pt() / lvgen.Pt();
      if (ratio >= 1.5) ratio = 1.5;
      if (ratio < 0.5) ratio = 0.5;
      p1d["GenOutPTRatio"]->Fill(ratio);

      double diffmag = (lvgenout - lvgen).M();
      p1d["GenOutDiffMag"]->Fill(diffmag);
    }
    // a->Tree_Fill();
  }
  a->SaveOutput();
}
