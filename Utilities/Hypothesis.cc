#ifndef HYPOTHESIS_CC
#define HYPOTHESIS_CC
#include <TLorentzVector.h>
#include <TString.h>
#include <TTree.h>
#include <vector>

class Hypothesis {
public:
  Hypothesis() { // TLorentzVectors should be stored inside of this class, pointers are only for branch use
    MakePointers();
  }
  TLorentzVector WPB, HadB, LF0, LF1, LepB, Lep, Neu, HadW, HadT, LepW, LepT, WP;
  TLorentzVector *pWPB, *pHadB, *pLF0, *pLF1, *pLepB, *pLep, *pNeu, *pHadW, *pHadT, *pLepW, *pLepT, *pWP;
  void MakePointers() {
    pWPB= &WPB;
    pHadB= &HadB;
    pLF0= &LF0;
    pLF1= &LF1;
    pLepB= &LepB;
    pLep= &Lep;
    pNeu= &Neu;
    pHadW= &HadW;
    pHadT= &HadT;
    pLepW= &LepW;
    pLepT= &LepT;
    pWP= &WP;
  }
  void Reset() {
    LF0 = LF1 = HadB = LepB = WPB = Lep = Neu = TLorentzVector();
    HadW = HadT = LepW = LepT = TLorentzVector();
    WP = TLorentzVector();
  }
  void Calculate(int SampleType) {
    HadW = LF0 + LF1;
    HadT = HadW + HadB;
    LepW = Lep + Neu;
    LepT = LepW + LepB;
    if (SampleType == 0) WP = WPB + HadT;
    else if (SampleType == 1) WP = WPB + LepT;
    else if (SampleType == 2) {
      if (WPB.DeltaR(HadT) > WPB.DeltaR(LepT)) WP = WPB + HadT;
      else WP = WPB + LepT;
    }
  }
  vector<TLorentzVector> Observables() {
    vector<TLorentzVector> out{LF0, LF1, HadB, LepB, WPB};
    return out;
  }
  vector<TLorentzVector*> pObservables() {
    vector<TLorentzVector*> out{pLF0, pLF1, pHadB, pLepB, pWPB};
    return out;
  }
  void SetLV(vector<TLorentzVector> vec) {
    const double &size = vec.size();
    if (size == 5 || size == 7) {
      LF0 = vec[0];
      LF1 = vec[1];
      HadB = vec[2];
      LepB = vec[3];
      WPB = vec[4];
    }
    if (size == 7) {
      Lep = vec[5];
      Neu = vec[6];
    }
    if (size < 5) {
      cout << endl << "Failed Setting Particles / Jets, Vector size less than 5" <<endl;
    }
    if (size == 6 || size > 7) {
      cout << endl << "Vector size = " << size <<" ? Is this intentional?" <<endl;
    }
  }
  void BookBranches(TTree* t, TString name, int complete) {
    t->Branch(name + "LF0",&(pLF0));
    t->Branch(name + "LF1",&(pLF1));
    t->Branch(name + "HadB",&(pHadB));
    t->Branch(name + "LepB",&(pLepB));
    t->Branch(name + "WPB",&(pWPB));
    t->Branch(name + "Lep",&(pLep));
    t->Branch(name + "Neu",&(pNeu));

    if (!complete) return;
    t->Branch(name + "WP",&(pWP));
    t->Branch(name + "HadT",&(pHadT));
    t->Branch(name + "HadW",&(pHadW));
    t->Branch(name + "LepT",&(pLepT));
    t->Branch(name + "LepW",&(pLepW));
  }


};

#endif
