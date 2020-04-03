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

void GenTemplate2d(int SampleType = 0, int irun = 1, int debug = 0) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "GenTemplate2d";

  GenTools *g1 = new GenTools("PFile","GenFileL1");
  GenTools *g2 = new GenTools("PFile","GenFileL2");
  JESTools *j  = new JESTools();
  TFile* PFile = new TFile("PFile/JESFile.root");
  j->ReadJESPlots(PFile);

  ROOTMini *m = new ROOTMini(j);

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);
  m->SetDebug(0);

  TH1F* EventCounter = new TH1F("EventCounter","Event Counter",10,-0.5,9.5);
  TH2F* WPrimeMasses2D = new TH2F("WPrimeMasses2D","W\' Masses; Mass if FL; Mass if LL",1000,0,1000,1000,0,1000);
  TH1F* BestPJes = new TH1F("BestPJes","Best Probability of ttbar hypothesis; Log(P)", 201,-20,0.1);
  TH2F* Probabilities = new TH2F("Probabilities","Best Probabilities; Log(Prob) if FL; Log(Prob) if LL",201,-20,0.1,200,-20,0.1);
  TH1F* WPrimeMassFL = new TH1F("WPrimeMassFL","W\' mass Reconstructed As FL",1000,0,1000);
  TH1F* WPrimeMassLL = new TH1F("WPrimeMassLL","W\' mass Reconstructed As LL",1000,0,1000);
  TH1F* WPrimeMassFlagFL = new TH1F("WPrimeMassFlagFL","W\' mass Flagged As FL",1000,0,1000);
  TH1F* WPrimeMassFlagLL = new TH1F("WPrimeMassFlagLL","W\' mass Flagged As LL",1000,0,1000);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    EventCounter->Fill(0); // Total Events
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
    if (a->RecoPass == -1) continue; //lepton != 1, jet < 5 : discard

    vector<TLorentzVector> Jets = a->LVOutPart; //Inputs
    vector<bool> BTags = a->GenOutBTags;
    m->SetLep(a->LVGenLep, a->LVGenNeu);
    if (BTags.size() != Jets.size()) {
      cout <<endl;
      cout << "Jets size = " << Jets.size() <<endl;
      cout << "BTags size = " << BTags.size() <<endl;
    }
    if (Jets.size() < 5) continue;
    EventCounter->Fill(1); // Event Passed selection

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
    double BestP = -1.;
    vector<int> BestPerm;
    vector<TLorentzVector> BestParticles;

    for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
      vector<int> thisperm = perms.at(iperm);
      double PBTags = j->CalcPFlavor(thisperm, BTags);
      vector<TLorentzVector> PermJets = j->GetPermutationLV(thisperm, Jets);
      double PJes = m->MinimizeP(PermJets);
      if (PJes == -1) continue;
      double PPerm = PJes * PBTags;
      if (PPerm > BestP) {
        BestP = PPerm;
        BestPerm = thisperm;
        m->ReCalcP(BestParticles);
      }
    } // Best Permutation ttbar hypothesis consists of is determined
    if (BestP == -1.) continue;

    EventCounter->Fill(2); //Event has a viable ttbar hypothesis
    BestPJes->Fill(log10(BestP));

    TLorentzVector hadt = BestParticles[0] + BestParticles[1] + BestParticles[2];
    TLorentzVector lept = BestParticles[3] + BestParticles[4] + BestParticles[5];

    // Deterimine the sample by together Top dPhi, Top dPt and W' top Pt.
    double pTopFL(1.), pTopLL(1.);
    double topdphi = fabs(hadt.DeltaPhi(lept));
    double topdpt = hadt.Pt() - lept.Pt();
    pTopFL *= g1->CalcP("TopdPhi",topdphi,0) * g1->CalcP("TopPtDiff",topdpt,0) * g1->CalcP("TopHadPt",hadt.Pt(),0);
    pTopLL *= g1->CalcP("TopdPhi",topdphi,1) * g1->CalcP("TopPtDiff",topdpt,1) * g1->CalcP("TopLepPt",lept.Pt(),1);

    double BestPFL(-1.), BestPLL(-1.);
    TLorentzVector BestFLWPB, BestLLWPB;
    vector<int> WPBCand = j->FindWPB(Jets.size(),BestPerm);
    for (unsigned iwpb = 0; iwpb < WPBCand.size(); ++iwpb) {
      double pWPB = 1.;
      TLorentzVector LVWPb = Jets.at(WPBCand.at(iwpb));
      double pWPBTag = j->CalcBTag(WPBCand.at(iwpb), BTags, true);
      double pWPBFL = g1->CalcP("WPBPt",LVWPb.Pt(),0) * g2->CalcP("WPdPhi", fabs(LVWPb.DeltaPhi(hadt)),0) * pTopFL;
      double pWPBLL = g1->CalcP("WPBPt",LVWPb.Pt(),1) * g2->CalcP("WPdPhi", fabs(LVWPb.DeltaPhi(hadt)),1) * pTopLL;

      if (pWPBFL > BestPFL) {
        BestPFL = pWPBFL;
        BestFLWPB = LVWPb;
      }
      if (pWPBLL > BestPLL) {
        BestPLL = pWPBLL;
        BestLLWPB = LVWPb;
      }
    }
    TLorentzVector FLWp = BestFLWPB + hadt;
    TLorentzVector LLWp = BestLLWPB + lept;
    double FLWpMass = FLWp.M();
    double LLWpMass = LLWp.M();
    WPrimeMasses2D->Fill(FLWpMass, LLWpMass);
    Probabilities->Fill(log10(BestP*BestPFL),log10(BestP*BestPLL));
    WPrimeMassFL->Fill(FLWpMass);
    WPrimeMassLL->Fill(LLWpMass);
    int SampleFlag = 0;
    if (BestPLL > BestPFL) SampleFlag = 1;
    if (!SampleFlag) {
      WPrimeMassFlagFL->Fill(FLWpMass);
    }
    else {
      WPrimeMassFlagLL->Fill(LLWpMass);
    }


  }
  a->SaveOutput();

}
