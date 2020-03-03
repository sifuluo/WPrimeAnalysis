#include "Utilities/Analyzer.cc"
// #include "Utilities/JetMatch.cc"
// #include "Utilities/JESTools.cc"

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

void MakeGenFile(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  // JESTools *b = new JESTools();
  TString savepath = "PFile/";
  TString savename = "GenFile";
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  TH1F* WPdR = new TH1F("WPdR", "dR of W\' daughters", 60, 0., 6.);
  TH1F* WPdPhi = new TH1F("WPdPhi", "dPhi = W\' daughters", 40, 0., 4.);
  TH1F* WPBPt = new TH1F("WPBPt", "Pt of W\' b", 100, 0, 1000);
  // Th1F* OtLeadingPt = new TH1F("OtLeadingPt", "Pt of leading particle other than W\'b",100,0,1000);
  TH1F* TopdR = new TH1F("TopdR", "dR of two tops", 60, 0., 6.);
  TH1F* TopdPhi = new TH1F("TopdPhi", "dPhi of tops", 40, 0., 4.);
  TH1F* TopPt = new TH1F("TopPt", "Pt of tops", 100, 0, 1000);
  TH1F* TopHadPt = new TH1F("TopHadPt", "Pt of hadronic tops", 100, 0, 1000);
  TH1F* TopLepPt = new TH1F("TopLepPt", "Pt of leptonic tops", 100, 0, 1000);
  TH1F* TopPtDiff = new TH1F("TopPtDiff", "Pt diff of HadT - LepT", 100, -500., 500.);
  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    TopdR->Fill(a->LVGenHadT.DeltaR(a->LVGenLepT));
    TopdPhi->Fill(a->LVGenHadT.DeltaPhi(a->LVGenLepT));
    TopHadPt->Fill(a->LVGenHadT.Pt());
    TopLepPt->Fill(a->LVGenLepT.Pt());
    TopPt->Fill(a->LVGenHadT.Pt());
    TopPt->Fill(a->LVGenLepT.Pt());
    TopPtDiff->Fill(a->LVGenHadT.Pt() - a->LVGenLepT.Pt());
    if (SampleType != 2) {
      WPdR->Fill(a->LVGenWPB.DeltaR(a->LVGenWPT));
      WPdPhi->Fill(a->LVGenWPB.DeltaPhi(a->LVGenWPT));
      WPBPt->Fill(a->LVGenWPB.Pt());
      // OtLeadingPt->Fill(a->LV)
    }
  }

  a->SaveOutput();

}
