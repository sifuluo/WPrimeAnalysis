#include "Utilities/Analyzer.cc"
#include "Utilities/JetMatch.cc"

#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TEfficiency.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>

using namespace std;
void MatchHypothesis(int SampleType = 0, int irun = 1, int debug = -2) {
  cout << "start"<< endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MatchedHypothesis";
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  TH1F* EventCounter = new TH1F("EventCounter","Number of Events;Passcode",10,-0.5,9.5);
  //TopMassDis
  TH1F* TopdR = new TH1F("TopdR","TopdR;dR",80,0,4);
  TH1F* TopPtDiff = new TH1F("TopPtDiff","W' t minus Other t, hadronic t minus leptonic t (bg)",400,-1000,1000);
  TH1F* WPTPt = new TH1F("WPTPt","W' t Pt",200,0,1000);
  TH1F* OtTPt = new TH1F("OtTPt","Other t Pt",200,0,1000);

  //Ws
  TH1F* WdR = new TH1F("WdR","W dR;dR",200,0,10);
  TH1F* WPtDiff = new TH1F("WPtDiff","W Pt difference",200,0,1000);

  //observables
  TH1F* WPBPt = new TH1F("WPBPt","W' b Pt",200,0,1000);
  TH1F* WPTBPt = new TH1F("WPTBPt","W' t b Pt",200,0,1000);
  TH1F* OtTBPt = new TH1F("OtTBPt","Other t b Pt",200,0,1000);
  TH1F* LFPt = new TH1F("LFPt","Light flavor Pt",200,0,1000);
  TH1F* WPBLeading = new TH1F("WPBLeading","W' b leading?",2,-0.5,1.5);
  TH1F* JetNumber = new TH1F("JetNumber","Number of Jets",20,-0.5,19.5);
  TH1F* BJetNumber = new TH1F("BJetNumber","Number of b-Jets",20,-0.5,19.5);
  TH1F* NBJetNumber = new TH1F("NBJetNumber","Number of Non-b-Jets",20,-0.5,19.5);
  TH2F* NBJetsVsBJets = new TH2F("NBJetsVsBJets","Non-B-Jets Vs B-Jets; N B-jets; B-jets",20,-0.5,19.5,20,-0.5,19.5);
  TH1F* LeadingIsB = new TH1F("LeadingIsB","Leading Jet is a B-jet",2,-0.5,1.5);

  //Additional
  TH1F* WPBDiffLeading = new TH1F("WPBDiffLeading","W' b minus other leading parton",400,-1000,1000);
  TH1F* WPBDiffLeadingB = new TH1F("WPBDiffLeadingB","W' b minus other leading b",400,-1000,1000);
  TH1F* WPTBDiffOtB = new TH1F("WPTBDiffOtB","W'-t-b minus the other t-b",400,-1000,1000);
  TH2F* WPBDiffTB = new TH2F("WPBDiffTB","W' b minus t-b; W'-b - W'-t-b; W'-b - Other t-b",200,-1000,1000,200,-1000,1000);
  TH2F* TPtVsLepPt = new TH2F("TPtVsLepPt","Leptonic t Pt Vs lepton Pt;lepton Pt;t Pt",100,0,1000,100,0,1000);
  TH2F* TPtVsMETPt = new TH2F("TPtVsMETPt","Leptonic t Pt Vs MET Pt;MET Pt;t Pt",100,0,1000,100,0,1000);
  TH1F* TPtMinusLepPt = new TH1F("TPtMinusLepPt","Leptonic t Pt minus Lepton Pt",200,0,1000);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->Jets.size() < 5 || a->GenW.size() != 2 || a->LVLeptons.size() != 1) continue;
    EventCounter->Fill(0);

    a->AssignGenParticles();
    a->GetRecoHypothesis();

  }
}
