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

  b->MakeJES();
  TProfile2D *pJetProfile2D = new TProfile2D("pJetProfile2D","Jet Response TProfile2D; recoPt; recoEta", 600, 0, 300, 60, 0, 6.0);
  vector<TProfile*> pJetProfiles;
  for (unsigned i = 0; i < b->EtaBins().size() - 1; ++i) {
    vector<double> bins = b->PtBins(i);
    int nb = bins.size();
    double arr[nb];
    copy(bins.begin(),bins.end(),arr);
    pJetProfiles.push_back(new TProfile(Form("pJetProfile_%d",i),"Jet Response Profile; recoPt;genPt / recoPt",nb-1,arr,"S"));
  }
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
      b->CalcBins(jeteta,jetpt);
      double rsp = lvgen.Pt() / jetpt;
      b->FillPlot(rsp);
      pJetProfile2D->Fill(jetpt,jeteta,rsp);
      pJetProfiles[b->GetiEta()]->Fill(jetpt,rsp);

    }

  }

  a->SaveOutput();

}
