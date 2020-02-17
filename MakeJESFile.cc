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

void MakePFile(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  JESTools *b = new JESTools();
  TString savepath = "NewPFile/";
  TString savename = "PFile";
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  b->MakeJESPlots();
  // TProfile2D *pJetProfile2D = new TProfile2D("pJetProfile2D","Jet Response TProfile2D; recoPt; recoEta", 600, 0, 300, 60, 0, 6.0);
  // vector<TProfile*> pJetProfiles;
  // for (unsigned i = 0; i < b->EtaBins().size() - 1; ++i) {
  //   vector<double> bins = b->PtBins(i);
  //   int nb = bins.size();
  //   double arr[nb];
  //   copy(bins.begin(),bins.end(),arr);
  //   pJetProfiles.push_back(new TProfile(Form("pJetProfile_%d",i),"Jet Response Profile; recoPt;genPt / recoPt",nb-1,arr,"S"));
  // }
  // b->MakeAdditionalPlots();

  TH1F* WPdR = new TH1F();
  TH1F* WPdPhi = new TH1F();
  TH1F* WPBPt = new TH1F();
  TH1F* TopdRBG = new TH1F();
  TH1F* TopdR = new TH1F();
  TH1F* WPTPtBG = new TH1F();
  TH1F* WPTPt = new TH1F();
  TH1F* WPTPtDiffBG = new TH1F();
  TH1F* WPTPtDiff = new TH1F();
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
      // pJetProfile2D->Fill(jetpt,jeteta,rsp);
      // pJetProfiles[bins.first]->Fill(jetpt,rsp);

    }
    if (SampleType == 2) {
      b->AddVector.at()
    };
    if (a->AssignGenParticles() == -1) continue;
    double dr = a->LVGenWPB.DeltaR(a->LVGenWPT);
    b->AddVector.at(0)->Fill(dr);
    double dphi = a->LVGenWPB.DeltaPhi(a->LVGenWPT);
    b->AddVector.at(1)->Fill(dphi);
    b->AddVector.at(2)->Fill(a->LVGenWPB.Pt());
    if (SampleType != 2) {

    }

  }

  a->SaveOutput();

}
