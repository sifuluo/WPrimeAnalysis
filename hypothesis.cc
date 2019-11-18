#include "Utilities/Analyzer.cc"

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

using namespace std;

void hypothesis(int SampleType = 0, int irun = 1, int debug = -2) {
  cout << "start"<< endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "hypothesis";

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);
  // a->ProcessEntries(100);

  //Numbers
  TH1F* EventCounter = new TH1F("EventCounter","Number of Events;Passcode",10,-0.5,9.5);
  TH1F* NumberW = new TH1F("NumberW","Number of W",10,-0.5,9.5);
  TH1F* Numbert = new TH1F("Numbert","Number of t",10,-0.5,9.5);
  TH1F* Numberb = new TH1F("Numberb","Number of b",10,-0.5,9.5);
  TH1F* NumberWP = new TH1F("NumberWP","Number of Wp",10,-0.5,9.5);

  //Tops
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
    //Numbers
    EventCounter->Fill(0);
    NumberW->Fill(a->GenW.size());
    Numbert->Fill(a->GenT.size());
    Numberb->Fill(a->GenB.size());
    NumberWP->Fill(a->GenWP.size());
    JetNumber->Fill(a->Jets.size());
    BJetNumber->Fill(a->BJets.size());
    NBJetNumber->Fill(a->NBJets.size());
    NBJetsVsBJets->Fill(a->NBJets.size(),a->BJets.size());

    // a->logfile << "Entry: " << entry << "  GenPass = " << a->GenPass<< endl;

    if (a->Jets.size() == 0) continue;
    // if (a->RecoPass < 0) continue;
    if (a->LVLeptons.size() != 1) {
      EventCounter->Fill(1);
      continue;
    }
    if (a->Jets.size() >= 5) {
      EventCounter->Fill(2);
      if (a->BJets.size() >=2 ) EventCounter->Fill(3);
    }
    // a->logfile << Form("For Entry: %i, RecoPass = %i, Jetsize = %i, Leptonsize = %i", (int) entry, (int)a->RecoPass, (int) a->Jets.size(), (int) (a->LVSoftLep.size() + a->LVLeptons.size()) ) << endl;
    // continue;
    if (a->GenPass < 0) continue;
    // EventCounter->Fill(2);
    // if (a->BJets.size() < 3) continue;
    // EventCounter->Fill(3);
    // a->PrintGenParticles("1");
    a->AssignGenParticles();

    //Tops
    TopdR->Fill(a->LVGenHadT.DeltaR(a->LVGenLepT));

    if (SampleType != 2) {
      WPTPt->Fill(a->LVGenWPT.Pt());
      OtTPt->Fill(a->LVGenOtT.Pt());
      TopPtDiff->Fill(a->LVGenWPT.Pt() - a->LVGenOtT.Pt());
    }
    else {
      WPTPt->Fill(a->LVGenHadT.Pt());
      OtTPt->Fill(a->LVGenLepT.Pt());
      TopPtDiff->Fill(a->LVGenHadT.Pt() - a->LVGenLepT.Pt());
    }

    //Ws
    WdR->Fill(a->LVGenHadW.DeltaR(a->LVGenLepW));
    WPtDiff->Fill(fabs(a->LVGenHadW.Pt() - a->LVGenLepW.Pt()));

    //observables
    LFPt->Fill(a->LVGenLFJet.at(0).Pt());
    LFPt->Fill(a->LVGenLFJet.at(1).Pt());
    if (SampleType != 2) {
      WPBPt->Fill(a->LVGenWPB.Pt());
      WPTBPt->Fill(a->LVGenWPTB.Pt());
      OtTBPt->Fill(a->LVGenOtB.Pt());
      bool leading = (a->LVGenWPB.Pt() > a->LVGenHadB.Pt() && a->LVGenWPB.Pt() > a->LVGenLepB.Pt() && a->LVGenWPB.Pt() > a->LVGenLFJet.at(0).Pt() && a->LVGenWPB.Pt() > a->LVGenLFJet.at(1).Pt() );
      if (leading) WPBLeading->Fill(1);
      else WPBLeading->Fill(0);
    }
    else {
      WPTBPt->Fill(a->LVGenHadB.Pt());
      OtTBPt->Fill(a->LVGenLepB.Pt());
    }

    map<double,int> JetOrder;
    for (unsigned it = 0; it < a->Jets.size(); ++it) {
      JetOrder.insert(pair<double,int>(a->Jets.at(it)->PT,it));
    }
    int nit = 0;
    for (auto it = JetOrder.rbegin(); it != JetOrder.rend(); ++it) {
      if (nit != (*it).second) a->logfile << "Jetorder was wrong at " << nit << " th jet of Event:" << entry <<endl;
      ++nit;
    }
    LeadingIsB->Fill(a->Jets.at((*(JetOrder.rbegin())).second)->BTag);

    //Additional
    vector<double> ptout;
    for (unsigned iout = 0; iout < a->LVGenOutSort.size(); ++iout) {
      ptout.push_back(a->LVGenOutSort.at(iout).Pt());
    }
    double wpbpt = a->LVGenWPB.Pt();
    double leadingpt = *max_element(ptout.begin(),ptout.begin()+4);
    double leadingbpt = *max_element(ptout.begin()+2,ptout.begin()+4);

    if (SampleType != 2) {
      WPBDiffLeading->Fill(wpbpt - leadingpt);
      WPBDiffLeadingB->Fill(wpbpt - leadingbpt);
      WPTBDiffOtB->Fill(a->LVGenWPTB.Pt() - a->LVGenOtB.Pt());
      WPBDiffTB->Fill(wpbpt - a->LVGenWPTB.Pt(), wpbpt - a->LVGenOtB.Pt());
    }
    else {
      WPTBDiffOtB->Fill(a->LVGenHadB.Pt() - a->LVGenLepB.Pt());
    }
    TPtVsLepPt->Fill(a->LVGenLep.Pt(), a->LVGenLepT.Pt());
    TPtVsMETPt->Fill(a->LVMET.Pt(), a->LVGenLepT.Pt());
    TPtMinusLepPt->Fill(a->LVGenLepT.Pt() - a->LVGenLep.Pt());
  }
  if (SampleType == 2)  WPTBDiffOtB->SetTitle("Hadronic b pt - leptonic b pt");
  a->SaveOutput();
};
