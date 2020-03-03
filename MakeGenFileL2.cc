#include "Utilities/Analyzer.cc"
// #include "Utilities/JetMatch.cc"
// #include "Utilities/JESTools.cc"
#include "Utilities/GenTools.cc"

#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TEfficiency.h>
#include <TString.h>

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>

using namespace std;

void MakeGenFileL2(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  // JESTools *b = new JESTools();
  TString savepath = "PFile/";
  TString savename = "GenFileL2";

  GenTools *g = new GenTools("PFile","GenFileL1");

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  TH1F* WPdPhi = new TH1F("WPdPhi", "dPhi = W\' daughters", 40, 0., 4.);
  TH1F* WPBPt = new TH1F("WPBPt", "Pt of W\' b", 100, 0, 1000);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->LVJets.size() < 5) continue;
    double wpdphi = fabs(a->LVGenWPB.DeltaPhi(a->LVGenWPT));
    double wpdphil1 = g->CalcP("WPdPhi",wpdphi,SampleType);
    double wpbpt = a->LVGenWPB.Pt();
    double wpbptl1 = g->CalcP("WPBPt",wpbpt,SampleType);
    WPdPhi->Fill(wpdphi, wpbptl1);
    WPBPt->Fill(wpbpt, wpdphil1);

  }
  a->SaveOutput();

}
