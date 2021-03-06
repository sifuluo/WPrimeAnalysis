#include "Utilities/Analyzer.cc"
#include "Utilities/JetMatch.cc"
#include "Utilities/JESTools.cc"

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

void MakeJESFile(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  JESTools *b = new JESTools();
  TString savepath = "PFile/";
  TString savename = "JESFile";
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  b->MakeJESPlots();

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
      double jetpt = lvjet.Pt();
      double jeteta = lvjet.Eta();
      if (jetpt < 30) continue;
      pair<int,int> bins = b->CalcBins(jeteta,jetpt);
      double rsp = lvgen.Pt() / jetpt;
      b->FillJESPlot(rsp, bins.first, bins.second);
    }
  }
  a->SaveOutput();

}
