#include "Utilities/Analyzer.cc"
#include "Utilities/JESTools.cc"
#include "Utilities/GenTools.cc"
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

void RecoMinimize(int SampleType = 0, int irun = 1, int testmode = 0) {
  cout << "start" <<endl;
  //SetUpUtilities
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "RecoMinimize";


  JESTools *b = new JESTools();
  GenTools *g = new GenTools("PFile","GenFile");
  TFile* JESFile = new TFile("PFile/JESFile.root");
  b->ReadJESPlots(JESFile);
  ROOTMini *m = new ROOTMini(b);
  // m->SetMinimizer();

  a->SetOutput(savepath,savename);
  a->DebugMode(testmode);
  m->SetDebug(0);
  m->SetTempOutput(0);
  bool output_ = false;

  TH1F* hFLBestProb = new TH1F("FLBestProb","FL Best Probability", 100,0.,1.);
  TH1F* hFLRecoWprimeMass = new TH1F("FLRecoWprimeMass", "FL Reconstructed W\' mass",1000,0,1000);
  TH2F* hFLPickVsProb = new TH2F("FLPickVsProb","FL Correct Pick Vs Probability; Probability; Parton Index",100,0,1.0,5,-0.5,4.5);
  TH1F* hFLSampleFlags = new TH1F("FLSampleFlags","FL Sample Type determined by minimized results",3,-0.5,2.5);

  TH1F* hLLBestProb = new TH1F("LLBestProb","LL Best Probability", 100,0.,1.);
  TH1F* hLLRecoWprimeMass = new TH1F("LLRecoWprimeMass", "LL Reconstructed W\' mass",1000,0,1000);
  TH2F* hLLPickVsProb = new TH2F("LLPickVsProb","LL Correct Pick Vs Probability; Probability; Parton Index",100,0,1.0,5,-0.5,4.5);
  TH1F* hLLSampleFlags = new TH1F("LLSampleFlags","LL Sample Type determined by minimized results",3,-0.5,2.5);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    // if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue; //lepton != 1, jet < 5 : discard

    //Set up inputs
    vector<TLorentzVector> Jets = a->LVJets;
    vector<bool> BTags = a->BTags;
    TLorentzVector Lepton = a->LVLeptons.at(0);
    TLorentzVector LVMET = a->LVMET;
    m->SetLep(Lepton, LVMET);

    vector< vector<int> > perms = b->MakePermutations(Jets.size());
    double BestProb = -1;
    int BestPerm = 0;
    int BestWPB = -1;
    int SampleFlag = -1;
    vector<TLorentzVector> BestParticles;
    BestParticles.clear();
    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      vector<int> thisperm = perms.at(iperm);
      double PBTags = b->CalcPFlavor(thisperm, BTags);
      vector<TLorentzVector> PermJets = b->GetPermutationLV(thisperm, Jets);
      double PermP = m->MinimizeP(PermJets);
      if (PermP == -1) continue;
      PermP *= PBTags;

      //Start to deal with W'B
      vector<TLorentzVector> PermParticles;
      m->ReCalcP(PermParticles);
      TLorentzVector hadt = PermParticles[0]+PermParticles[1]+PermParticles[2];
      TLorentzVector lept = PermParticles[3]+PermParticles[4]+PermParticles[5];

      // Probability to be samples by Tops
      double pTopFL = 1.;
      double pTopLL = 1.;
      double pTopBG = 1.; // Background likelihood correction

      double topdr = hadt.DeltaR(lept);
      pTopFL *= g->CalcP("TopdR",topdr,0);
      pTopLL *= g->CalcP("TopdR",topdr,1);
      pTopBG *= g->CalcP("TopdR",topdr,2);

      double topptdiff = hadt.Pt() - lept.Pt();
      pTopFL *= g->CalcP("TopPtDiff",topptdiff,0);
      pTopLL *= g->CalcP("TopPtDiff",topptdiff,1);
      pTopBG *= g->CalcP("TopPtDiff",topptdiff,2);

      pTopFL *= pTopBG;
      pTopLL *= pTopBG;

      // W'b Searching
      vector<int> WPBCand = b->FindWPB(Jets.size(),thisperm);
      // Processing variables
      double PermPWPB = -1;
      int PermFlag = -1;
      int PermiWPB = -1;
      TLorentzVector LVPermWPB = TLorentzVector();

      for (unsigned iwpb = 0; iwpb < WPBCand.size(); ++iwpb) {
        double thispwpb;
        int thisflag;
        TLorentzVector ThisWPB = Jets.at(WPBCand.at(iwpb));
        double pWPBTag = b->CalcBTag(iwpb, BTags, true);
        double phadwpb = g->CalcP("WPdPhi",ThisWPB.DeltaPhi(hadt),0) * g->CalcP("WPBPt",ThisWPB.Pt(),0);
        phadwpb *= pTopFL * pWPBTag;
        double plepwpb = g->CalcP("WPdPhi",ThisWPB.DeltaPhi(lept),1) * g->CalcP("WPBPt",ThisWPB.Pt(),1);
        plepwpb *= pTopLL * pWPBTag;
        if (phadwpb > plepwpb) {
          thisflag = 0;
          thispwpb = phadwpb;
        }
        else {
          thisflag = 1;
          thispwpb = plepwpb;
        }
        if (thispwpb > PermPWPB) {
          PermFlag = thisflag;
          PermPWPB = thispwpb;
          PermiWPB = WPBCand.at(iwpb);
          LVPermWPB = ThisWPB;
        }
      }
      //WPB for this permutation is calculated
      PermParticles.push_back(LVPermWPB);
      PermP *= PermPWPB;
      // PermParticles orders as LFJet1, LFjet2, HadB, LepB, lepton, neutrino, w'b

      // if (iperm == 0 ) BestProb = PermP;
      if (PermP > BestProb || iperm == 0) {
        BestProb = PermP;
        BestPerm = iperm;
        BestWPB = PermiWPB;
        SampleFlag = PermFlag;
        BestParticles = PermParticles;
      }
    }

    if (SampleFlag == 0) {
      hFLBestProb->Fill(BestProb);
      double wpmass = (BestParticles[0]+BestParticles[1]+BestParticles[2]+BestParticles[6]).M();
      hFLRecoWprimeMass->Fill(wpmass);
    }
    else if (SampleFlag == 1) {
      hLLBestProb->Fill(BestProb);
      double wpmass = (BestParticles[0]+BestParticles[1]+BestParticles[2]+BestParticles[6]).M();
      hLLRecoWprimeMass->Fill(wpmass);
    }

  }
  a->SaveOutput();

}
