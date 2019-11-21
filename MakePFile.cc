#include "Utilities/Analyzer.cc"
#include "Utilities/JetMatch.cc"
#include "Utilities/EtaPtBins.cc"

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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>

using namespace std;

void MakePFile(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  EtaPtBins *b = new EtaPtBins();
  TString savepath = "PFile/";
  TString savename = "PFile";
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  TProfile2D *pJetProfile = new TProfile2D("pJetProfile2D","Jet Response TProfile2D; recoPt; recoEta", 600, 0, 300, 60, 0, 6.0);
  b->MakeJES();
  b->GetPlot(1.0,40.)->Fill(2);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    // a->AssignGenParticles();
    a->MatchJets();

    for (auto mapit = a->JetMatchMap.begin(); mapit != a->JetMatchMap.end(); ++mapit) {
      TLorentzVector lvgen = a->LVOutPart.at((*mapit).first);
      TLorentzVector lvjet;
      int j1 = (*mapit).second.at(0);
      int j2 = (*mapit).second.at(1);
      if (j1 == j2) lvjet = a->LVJets.at(j1);
      else lvjet = a->LVJets.at(j1) + a->LVJets.at(j2);
      b->GetPlot(lvjet.Eta(), lvjet.Pt())->Fill(lvgen.Pt() / lvjet.Pt());
      pJetProfile->Fill(lvjet.Pt(),lvjet.Eta(),lvgen.Pt() / lvjet.Pt());
    }

  }

  a->SaveOutput();

}
