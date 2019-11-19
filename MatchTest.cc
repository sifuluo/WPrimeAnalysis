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

void MatchTest(int SampleType = 0, int irun = 1, int debug = -2) {
  cout << "start"<< endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MatchTest";
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  TH1F* EventCounter = new TH1F("EventCounter","Number of Events;Passcode",10,-0.5,9.5);
  TH1F* JetMatchRatio = new TH1F("JetMatchRatio","Jet Match Pt Ratio",600,0.,6.);
  TH1F* AdvJetMatchRatio = new TH1F("AdvJetMatchRatio","Advanced Jet Match Pt Ratio",600,0.,6.);
  TH1F* AdvJetMatch2Rate = new TH1F("AdvJetMatch2Rate","# of jet Advanced Jet Match matches to",2,0.5,2.5);
  TH1F* NoMatchPar = new TH1F("NoMatch","Particle has no Match in JetMatch; Other parton, LFjet1, LFjet2, Had b, Lep b, W' b", 6,-1.5,4.5);
  TH1F* NoAdvMatchPar = new TH1F("NoAdvMatch","Particle has no Match in Advanced Match; Other parton, LFjet1, LFjet2, Had b, Lep b, W' b",6,-1.5,4.5);
  TH2F* NoMatchParPtVsMinDr = new TH2F("NoMatchParPtVsMinDr","Particle Attr. has no Match in Match;Pt;Min#Delta R with other Parton",100,0.,500., 60,0.,6.);
  TH2F* NoAdvMatchParPtVsMinDr = new TH2F("NoAdvMatchParPtVsMinDr","Particle Attr. has no Match in Advanced Match;Pt;Min#Delta R with other Parton",100,0.,500., 60,0.,6.);
  TH1F* LeptonResolution = new TH1F("LeptonResolution", "Lepton Resoluiton; Lepton Pt/ GenLepton Pt", 600,0.,6.);

  TH2F* LowerMergedJets = new TH2F("LowerMergedJets", "Merged Jets of lower Pt Attr. ; Pt ; DeltaR",600,0.,600.,60,0.,6.);
  TH2F* HigherMergedJets = new TH2F("HigherMergedJets", "Merged Jets of higher Pt Attr. ; Pt ; DeltaR",600,0.,600.,60,0.,6.);
  // Missing plots are:
  //Lost match in Merged Jets.
  //Additional Jets in Merged Jets
  int leptonb = 0;
  int otherb = 0;
  int advleptonb = 0;
  int advotherb = 0;
  int wnot2 = 0;
  int nojet = 0;
  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    EventCounter->Fill(0);
    if (a->Jets.size() == 0) {
      nojet++;
      continue;
    }
    if (a->GenW.size() != 2) {
      wnot2++;
      continue;
    }
    //Preselections
    if (a->LVLeptons.size() != 1 ) {
      continue;
    }
    EventCounter->Fill(1);
    if (a->LVJets.size() < 5) {
      continue;
    }
    EventCounter->Fill(2);

    a->AssignGenParticles();
    vector<TLorentzVector> allout = a->LVOutPart;
    vector<TLorentzVector> gen = a->LVGenOutSort;
    // GenOutSort ordered as LFjet1, LFjet2, HadB, LepB, W'B
    vector<TLorentzVector> reco = a->LVJets;
    double maxdr;
    vector<int> LVOutPartType; // vector of allout particle types(which particle in hypothesis)
    LVOutPartType.clear();

    //find if a particle in LVOutPart is one of desired sorted particles
    for (unsigned  i1 = 0; i1 < allout.size(); ++i1) {
      int igen = -1;
      TLorentzVector v1 = allout.at(i1);
      for (unsigned i2 = 0; i2 < gen.size(); ++i2) {
        TLorentzVector v2 = gen.at(i2);
        if (v1 == v2) {
          igen = i2;
          LVOutPartType.push_back(i2);
          break;
        }
      }
      if (igen == -1) LVOutPartType.push_back(-1);
    }
    // Testing part for LVOutPartType
    // a->logfile << "The match list is: ";
    // for (unsigned i1 = 0; i1 < allout.size(); ++i1) {
    //   a->logfile << LVOutPartType.at(i1)<<" , ";
    // }
    // a->logfile <<endl;

    map<int,int> jetmatch = JetMatch(allout, reco, maxdr, true, false, true);
    map<int, vector<int> > advjetmatch = AdvJetMatch(allout, reco, maxdr, true, false, true);


    for (unsigned it = 0; it <allout.size(); ++it) {
      auto jetmatchit = jetmatch.find(it);
      if (jetmatchit != jetmatch.end()){
        int match = (*jetmatchit).second;
        JetMatchRatio->Fill(reco.at(match).Pt() / allout.at(it).Pt());
      }
      else {
        NoMatchPar->Fill(LVOutPartType.at(it));
        if (LVOutPartType.at(it) != -1){
          if (gen.at(LVOutPartType.at(it)) != allout.at(it)) a->logfile << "Wrong Match!!!" <<endl;
        }
        if (LVOutPartType.at(it) == 3) leptonb++;
        else otherb++;
        double minDr = 100;
        int closegen = -1;
        bool init = true;
        for (unsigned itgen = 0; itgen < allout.size(); ++itgen) {
          if (itgen == it) continue;
          double dr = allout[it].DeltaR(allout[itgen]);
          if (init) {
            init = false;
            minDr = dr;
            closegen = itgen;
            continue;
          }
          if (dr < minDr) {
            minDr = dr;
            closegen = itgen;
          }
        }
        NoMatchParPtVsMinDr->Fill(allout[it].Pt(), minDr);
      }

      auto advjetmatchit = advjetmatch.find(it);
      if (advjetmatchit != advjetmatch.end()){
        int match1 = (*advjetmatchit).second.at(0);
        int match2 = (*advjetmatchit).second.at(1);
        if (match1 == match2) {
          AdvJetMatchRatio->Fill( (reco.at(match1)).Pt() / allout.at(it).Pt());
          AdvJetMatch2Rate->Fill(1);
        }
        else {
          AdvJetMatchRatio->Fill( (reco.at(match1)+reco.at(match2)).Pt() / allout.at(it).Pt());
          AdvJetMatch2Rate->Fill(2);
          double pt1 = reco.at(match1).Pt();
          double pt2 = reco.at(match2).Pt();
          double dr = reco.at(match1).DeltaR(reco.at(match2));
          if (pt1 >= pt2){
            LowerMergedJets->Fill(pt2,dr);
            HigherMergedJets->Fill(pt1,dr);
          }
          else {
            LowerMergedJets->Fill(pt1,dr);
            HigherMergedJets->Fill(pt2,dr);
          }
        }

      }
      else {
        NoAdvMatchPar->Fill(LVOutPartType.at(it));
        double minDr = 100;
        int closegen = -1;
        bool init = true;
        for (unsigned itgen = 0; itgen < allout.size(); ++itgen) {
          if (itgen == it) continue;
          double dr = allout[it].DeltaR(allout[itgen]);
          if (init) {
            init = false;
            minDr = dr;
            closegen = itgen;
            continue;
          }
          if (dr < minDr) {
            minDr = dr;
            closegen = itgen;
          }
        }
        NoAdvMatchParPtVsMinDr->Fill(allout[it].Pt(), minDr);
      }
    }

    // if (a->LVLeptons.size() != (a->GenE.size() + a->GenMu.size())) {
    //   a->logfile << "LVLepton size = " << a->LVLeptons.size() << " GenLep Size = " << a->GenE.size() + a->GenMu.size() << " lvoutpart size = " << a->LVOutPart.size() << endl;
    // }
    // a->logfile << endl;
    if (a->LVLeptons.size() != 1) continue;
    if (a->GenE.size() == 1) {
      LeptonResolution->Fill(a->LVLeptons.at(0).Pt() / a->GenParticles.at(a->GenE.at(0))->PT);
    }
    else if (a->GenMu.size() ==1) {
      LeptonResolution->Fill(a->LVLeptons.at(0).Pt() / a->GenParticles.at(a->GenMu.at(0))->PT);
    }
  }
  a->logfile << "leptonb : " <<leptonb << " others : " <<otherb <<endl;
  a->logfile << "nojet: " <<nojet << " wnot2 : " <<wnot2 <<endl;
  a->SaveOutput();
}
