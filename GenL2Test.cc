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

void GenL2Test(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  // JESTools *b = new JESTools();
  TString savepath = "results/";
  TString savename = "GenL2Test";

  GenTools *g1 = new GenTools("PFile","GenFileL1");
  GenTools *g2 = new GenTools("PFile","GenFileL2");

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);
  g2->GetHist("WPBPt",0)->Clone();
  g2->GetHist("WPdPhi",0)->Clone();
  TH1F* L1Multi = new TH1F("L1Multi","L1 distribution multiplied", 150,0.,1.5);
  TH1F* SeqWPdPhi = new TH1F("SeqWPdPhi","L1 WPBPt * L2 WPdPhi", 150, 0., 1.5);
  TH1F* SeqWPBPt = new TH1F("SeqWPBPt","L1 WPdPhi * L2 WPBPt", 150, 0., 1.5);
  TH1F* SeqDiff = new TH1F("SeqDiff","L1 WPBPt * L2 WPdPhi - L1WPdPhi * L2 WPBPt", 300, -1.5, 1.5);
  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->LVJets.size() < 5) continue;

    double dphi = fabs(a->LVGenWPB.DeltaPhi(a->LVGenWPT));
    double pt = a->LVGenWPB.Pt();
    double ptl1 = g1->CalcP("WPBPt",pt,SampleType);
    double ptl2 = g2->CalcP("WPBPt",pt,SampleType);
    double dphil1 = g1->CalcP("WPdPhi",dphi,SampleType);
    double dphil2 = g2->CalcP("WPdPhi",dphi,SampleType);

    L1Multi->Fill(ptl1*dphil1);
    SeqWPdPhi->Fill(ptl1*dphil2);
    SeqWPBPt->Fill(ptl2*dphil1);
    SeqDiff->Fill(ptl1*dphil2 - ptl2*dphil1);

  }
  a->SaveOutput();

}
