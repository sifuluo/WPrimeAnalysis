#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "analyzer.cc"
#include "tools.cc"

#include <TROOT.h>
#include <TString.h>
#include <TVector2.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <TProfile.h>
#include <TEfficiency.h>
#include <string>

using namespace std;

void wprime(int Sampletype, int irun = 0) {
  // gSystem->Load("libDelphes");
  vector<TString> samplebasepaths;
  samplebasepaths.push_back("/fdata/hepx/store/user/aoverton0342/madGraph/");

  vector<TString> folders;
  TString savepath = "/fdata/hepx/store/user/siluo/wprime/results/";
  TString savename;
  if (Sampletype == 0) {
    folders.push_back("TDual_FormerLeptonic/");
    savename = "TDual_FormerLeptonic_Optimized";
  }
  else if (Sampletype == 1) {
    folders.push_back("TDual_LatterLeptonic/");
    savename = "TDual_LatterLeptonic_Optimized";
  }
  else {
    folders.push_back("ttbar2/");
    savename = "Background";
  }

  Analyzer *a = new Analyzer(samplebasepaths,folders,irun, 30);
  a->SetOutput(savepath,savename);
  // bool CreateProbability = true;
  // a->GetProbability(CreateProbability);
  // if (CreateProbability) return;
  // ofstream &logfile = a->logfile;
  // TFile* PFile= a->PFile;
  bool debug = false;
  if (debug) {
    a->SetStartEntry(0);
    a->ProcessEntries(5000);
  }
  TH2F* JetMatches = new TH2F("JetMatches","JetMatches;NoMatch,PartGenMatch,GenGenJetMatch; Entries/100",3,-0.5,2.5,80,-0.5,79.5);
  // TH1F* PartonDeltaR = new TH1F("PartonDeltaR","Parton Delta R; Delta R; Percentage",1000,0.,100.);

  // Long64_t StartEntry = a->GetStartEntry();
  // Long64_t EndEntry = a->GetEndEntry();

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry){
    if (a->ReadEvent(entry,debug) == -2) continue;
    // a->AssignGenParticles();
    a->MakeMatchMaps();
    double fff = floor(entry/100);

    map<int,int>::iterator it1, it2;
    // int c1(0),c2(0),c3(0);
    for (unsigned i = 0; i < a->LVOutPart.size(); ++i) {
      it1 = a->OutPartGenJetMap.find(i);
      if (it1 == a->OutPartGenJetMap.end() ) {
        JetMatches->Fill(0.,fff);
      }
      else {
        it2 = a->GenJetJetMap.find(it1->second);
        if (it2 == a->GenJetJetMap.end() ) {
          JetMatches->Fill(1.,fff);
        }
        else {
          JetMatches->Fill(2.,fff);
        }
      }
    }

    // vector<TLorentzVector> out = a->LVOutPart;
    // for (unsigned iout = 0; iout < out.size(); ++iout) {
    //   for (unsigned iout2 = 0; iout2 < out.size(); ++iout2) {
    //     if (iout2 == iout) continue;
    //     PartonDeltaR->Fill(out[iout].DeltaR(out[iout2]));
    //   }
    // }

  }//Looped over all entries

  // JetMatches->Scale(1./JetMatches->GetEntries());
  a->SaveOutput();
}
