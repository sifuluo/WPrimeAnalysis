#include "Utilities/Analyzer.cc"
// #include "Utilities/JetMatch.cc"
// #include "Utilities/JESTools.cc"
#include "Utilities/GenTools.cc"
#include "Utilities/ROOTMini.cc"

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

void GenTemplateTest(int SampleType = 0, int irun = 1, int debug = 0) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "GenTemplateTest";

  GenTools *g1 = new GenTools("PFile","GenFileL1");
  GenTools *g2 = new GenTools("PFile","GenFileL2");
  JESTools *j  = new JESTools();
  TFile* PFile = new TFile("PFile/JESFile.root");
  j->ReadJESPlots(PFile);

  ROOTMini *m = new ROOTMini(j);

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);
  m->SetDebug(0);

  TH1F* EventCounter = new TH1F("EventCounter","Event Counter",10, -0.5,9.5);
  TH1F* JESCorRate = new TH1F("JESCorRate","Correct Rate with only JES",6,-0.5,5.5);
  TH1F* JESCorParticle = new TH1F("JESCorParticle","Correct Particle with only JES",6,-0.5,5.5);
  TH1F* FullCorRate = new TH1F("FullCorRate","Correct Rate with all discriminants",6,-0.5,5.5);
  TH1F* FullCorParticle = new TH1F("FullCorParticle","Correct Particle with all discriminants",6,-0.5,5.5);
  TH1F* hFLBestProb = new TH1F("FLBestProb","FL Best Probability;Log(P)", 201,-20,0.1);
  TH1F* hLLBestProb = new TH1F("LLBestProb","LL Best Probability;Log(P)", 201,-20,0.1);
  TH1F* hFLRecoNeudR = new TH1F("hFLRecoNeudR","FL Reconstructed Neutrino dR",40,0.,4.);
  TH1F* hLLRecoNeudR = new TH1F("hLLRecoNeudR","LL Reconstructed Neutrino dR",40,0.,4.);
  TH1F* hRecoHadTdR = new TH1F("hRecoHadTdR","Reconstructed hadronic t dR with Gen",40,0.,4.);
  TH1F* hRecoHadTMass = new TH1F("hRecoHadTMass","Reconstructed hadronic t mass",400,0,400);
  TH1F* hRecoLepTdR = new TH1F("hRecoLepTdR","Reconstructed leptonic t dR with Gen",40,0.,4.);
  TH1F* hRecoLepTMass = new TH1F("hRecoLepTMass","Reconstructed leptonic t mass",400,0,400);
  TH1F* hFLRecoLepWMass = new TH1F("hFLRecoLepWMass","FL Reconstructed leptonic W mass", 150,0.,150.);
  TH1F* hLLRecoLepWMass = new TH1F("hLLRecoLepWMass","LL Reconstructed leptonic W mass", 150,0.,150.);
  TH1F* hFLRecoLepWdR = new TH1F("hFLRecoLepWdR","FL Reconstructed leptonic W dR", 40,0.,4.);
  TH1F* hLLRecoLepWdR = new TH1F("hLLRecoLepWdR","LL Reconstructed leptonic W dR", 40,0.,4.);
  TH1F* hFLRecoWprimeMass = new TH1F("FLRecoWprimeMass", "FL Reconstructed W\' mass",1000,0,1000);
  TH1F* hLLRecoWprimeMass = new TH1F("LLRecoWprimeMass", "LL Reconstructed W\' mass",1000,0,1000);
  TH1F* hCorRecoWprimeMass = new TH1F("CorRecoWprimeMass", "ReCorrect flaging Reconstructed W\' mass",1000,0,1000);
  // TH1F* LeadingAssignCorrect = new TH1F("LeadingAssignCor","Leading Particle Correct Assignment",5,-0.5,4.5);
  // TH1F* LeadingAssignWrong = new TH1F("LeadingAssignWro","Leading Particle Wrong Assignment",5,-0.5,4.5);


  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    EventCounter->Fill(0);
    if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue; //lepton != 1, jet < 5 : discard

    vector<TLorentzVector> Jets = a->LVOutPart;
    vector<bool> BTags = a->GenOutBTags;
    m->SetLep(a->LVGenLep, a->LVGenNeu);
    if (BTags.size() != Jets.size()) {
      cout <<endl;
      cout << "Jets size = " << Jets.size() <<endl;
      cout << "BTags size = " << BTags.size() <<endl;
    }
    if (Jets.size() < 5) continue;

    EventCounter->Fill(1);

    // Find the correct permutation using Gen info
    vector<int> GenOutQuark = a->GenOutQuark;
    vector<int> GenOutSort = a->GenOutSort;
    vector<int> hypoindex;
    for (unsigned i = 0; i < GenOutSort.size(); ++i) {
      int hypo = GenOutSort.at(i);
      if (hypo == -1) {
        hypoindex.push_back(-1);
        continue;
      }
      for (unsigned j = 0; j < GenOutQuark.size(); ++j) {
        if (GenOutQuark[j] == hypo) {
          hypoindex.push_back(j);
          break;
        }
      }
    }
    // hypoindex points to the correct permutation in Jets, Meaning hypoindex is the correct perm.

    vector< vector<int> > perms = j->MakePermutations(Jets.size());
    double BestP = -1;
    vector<int> BestPerm;
    int SampleFlag = -1;
    vector<TLorentzVector> BestParticles;
    double BestPJes = -1;
    vector<int> BestPermJES;

    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      vector<int> thisperm = perms.at(iperm);
      double PBTags = j->CalcPFlavor(thisperm, BTags);
      vector<TLorentzVector> PermJets = j->GetPermutationLV(thisperm, Jets);
      double PJes = m->MinimizeP(PermJets);
      if (PJes == -1) continue;
      // cout << endl<<"Worked out PJes = " << PJes <<endl;
      double PPerm = PJes * PBTags;

      if (PPerm > BestPJes) {
        BestPJes = PPerm;
        BestPermJES = thisperm;
      }

      vector<TLorentzVector> PermParticles;
      m->ReCalcP(PermParticles);
      TLorentzVector hadt = PermParticles[0]+PermParticles[1]+PermParticles[2];
      TLorentzVector lept = PermParticles[3]+PermParticles[4]+PermParticles[5];

      double pTopFL(1.), pTopLL(1.), pTopBG(1.);
      double topdphi = fabs(hadt.DeltaPhi(lept));
      double topdpt = hadt.Pt() - lept.Pt();

      double pTopFLdPhi = g1->CalcP("TopdPhi",topdphi,0);
      double pTopLLdPhi = g1->CalcP("TopdPhi",topdphi,1);
      double pTopBGdPhi = g1->CalcP("TopdPhi",topdphi,2);
      pTopFLdPhi -= pTopBGdPhi;
      pTopLLdPhi -= pTopBGdPhi;
      if (pTopFLdPhi < 0) pTopFLdPhi = 0;
      if (pTopLLdPhi < 0) pTopLLdPhi = 0;
      pTopFL *= pTopFLdPhi;
      pTopLL *= pTopLLdPhi;

      double pTopFLPtDiff = g1->CalcP("TopPtDiff",topdpt,0);
      double pTopLLPtDiff = g1->CalcP("TopPtDiff",topdpt,1);
      double pTopBGPtDiff = g1->CalcP("TopPtDiff",topdpt,2);
      pTopFLPtDiff -= pTopBGPtDiff;
      pTopLLPtDiff -= pTopBGPtDiff;
      if (pTopFLPtDiff < 0) pTopFLPtDiff = 0;
      if (pTopLLPtDiff < 0) pTopLLPtDiff = 0;
      pTopFL *= pTopFLPtDiff;
      pTopLL *= pTopLLPtDiff;

      vector<int> WPBCand = j->FindWPB(Jets.size(),thisperm);
      double BestpWPB = 0;
      int BestWPBflag = -1;
      int BestWPB = -1;
      TLorentzVector BestLVWPB;

      for (unsigned iwpb = 0; iwpb < WPBCand.size(); ++iwpb) {
        double pWPb = 1;
        int wpbflag = -1;
        TLorentzVector LVWPb = Jets.at(WPBCand.at(iwpb));
        double pWPBTag = j->CalcBTag(iwpb, BTags, true);
        double pWPBFL = g1->CalcP("WPBPt",LVWPb.Pt(),0) * g2->CalcP("WPdPhi", fabs(LVWPb.DeltaPhi(hadt)),0) * pTopFL;
        double pWPBLL = g1->CalcP("WPBPt",LVWPb.Pt(),1) * g2->CalcP("WPdPhi", fabs(LVWPb.DeltaPhi(hadt)),1) * pTopLL;
        if (pWPBFL > pWPBLL) {
          wpbflag = 0;
          pWPb = pWPBFL * pWPBTag;
        }
        else {
          wpbflag = 1;
          pWPb = pWPBLL * pWPBTag;
        }
        if (pWPb > BestpWPB) {
          BestpWPB = pWPb;
          BestWPBflag = wpbflag;
          BestWPB = iwpb;
          BestLVWPB = LVWPb;
        }
      }

      PPerm *= BestpWPB;
      if (PPerm > BestP) {
        thisperm.push_back(BestWPB);
        PermParticles.push_back(BestLVWPB);
        BestP = PPerm;
        BestPerm = thisperm;
        SampleFlag = BestWPBflag;
        BestParticles = PermParticles;
      }
    }

    if (BestPJes > 0) {
      int JesCorCount = 0;
      if (hypoindex[0] == -1 || hypoindex[0] == BestPermJES[0] || hypoindex[0] == BestPermJES[1]) {
        JesCorCount++;
        JESCorParticle->Fill(0);
      }
      if (hypoindex[1] == -1 || hypoindex[1] == BestPermJES[0] || hypoindex[1] == BestPermJES[1]) {
        JesCorCount++;
        JESCorParticle->Fill(1);
      }
      if (hypoindex[2] == -1 || hypoindex[2] == BestPermJES[2]) {
        JesCorCount++;
        JESCorParticle->Fill(2);
      }
      if (hypoindex[3] == -1 || hypoindex[3] == BestPermJES[3]) {
        JesCorCount++;
        JESCorParticle->Fill(3);
      }
      JESCorRate->Fill(JesCorCount);
    }

    if (BestP == -1) continue;
    int FullCorCount = 0;
    if (hypoindex[0] == -1 || hypoindex[0] == BestPerm[0] || hypoindex[0] == BestPerm[1]) {
      FullCorCount++;
      FullCorParticle->Fill(0);
    }
    if (hypoindex[1] == -1 || hypoindex[1] == BestPerm[0] || hypoindex[1] == BestPerm[1]) {
      FullCorCount++;
      FullCorParticle->Fill(1);
    }
    if (hypoindex[2] == -1 || hypoindex[2] == BestPerm[2]) {
      FullCorCount++;
      FullCorParticle->Fill(2);
    }
    if (hypoindex[3] == -1 || hypoindex[3] == BestPerm[3]) {
      FullCorCount++;
      FullCorParticle->Fill(3);
    }
    FullCorRate->Fill(FullCorCount);

    TLorentzVector hadt = BestParticles[0]+BestParticles[1]+BestParticles[2];
    TLorentzVector lept = BestParticles[3]+BestParticles[4]+BestParticles[5];
    hRecoHadTdR->Fill(hadt.DeltaR(a->LVGenHadT));
    hRecoHadTMass->Fill(hadt.M());
    hRecoLepTdR->Fill(lept.DeltaR(a->LVGenLepT));
    hRecoLepTMass->Fill(lept.M());

    if (SampleFlag == 0) {
      hFLBestProb->Fill(log10(BestP));
      double wpmass = (BestParticles[0]+BestParticles[1]+BestParticles[2]+BestParticles[6]).M();
      hFLRecoWprimeMass->Fill(wpmass);
      hFLRecoNeudR->Fill(BestParticles[5].DeltaR(a->LVGenNeu));
      hFLRecoLepWMass->Fill((BestParticles[4]+BestParticles[5]).M());
      hFLRecoLepWdR->Fill((BestParticles[4]+BestParticles[5]).DeltaR(a->LVGenLepW));
    }
    else if (SampleFlag == 1) {
      hLLBestProb->Fill(log10(BestP));
      double wpmass = (BestParticles[3]+BestParticles[4]+BestParticles[5]+BestParticles[6]).M();
      hLLRecoWprimeMass->Fill(wpmass);
      hLLRecoNeudR->Fill(BestParticles[5].DeltaR(a->LVGenNeu));
      hLLRecoLepWMass->Fill((BestParticles[4]+BestParticles[5]).M());
      hLLRecoLepWdR->Fill((BestParticles[4]+BestParticles[5]).DeltaR(a->LVGenLepW));
    }

  }
  a->SaveOutput();

}
