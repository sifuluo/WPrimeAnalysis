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
  vector<TLorentzVector*> pHypoVector;
  vector<TString> HypoNames;
  double FL_WPMass, LL_WPMass;
  vector<bool> BTags;
  int SampleType;

  void SetSampleType(int st_) {
    SampleType = st_;
  }
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
    pHypoVector = vector<TLorentzVector*>{pLF0,pLF1,pHadB,pLepB,pWPB,pLep,pNeu,pWP,pHadT,pHadW,pLepT,pLepW};
    HypoNames = vector<TString>{"LF0","LF1","HadB","LepB","WPB","Lep","Neu","WP","HadT","HadW","LepT","LepW"};
    //                            0     1      2      3     4     5     6     7    8      9      10     11
  }
  void Reset(int st_ = -1) {
    if (st_ != -1) SampleType = st_;
    LF0 = LF1 = HadB = LepB = WPB = Lep = Neu = TLorentzVector();
    HadW = HadT = LepW = LepT = TLorentzVector();
    WP = TLorentzVector();
    FL_WPMass = LL_WPMass = 0;
    BTags.clear();
  }
  void Calculate(int st_ = -1) {
    if (st_ != -1) SampleType = st_;
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
    FL_WPMass = (WPB + HadT).M();
    LL_WPMass = (WPB + LepT).M();
  }
  vector<TLorentzVector> Observables() {
    vector<TLorentzVector> out{*pLF0, *pLF1, *pHadB, *pLepB, *pWPB};
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
  void BookBranches(TTree* t, TString name, int complete = true) {
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
  bool AllFilled() {
    TLorentzVector lv0 = TLorentzVector();
    if (LF0 == lv0 || LF1 == lv0 || HadB == lv0 || LepB == lv0 || (WPB == lv0 && SampleType != 2)) return false;
    else return true;
  }


};

#endif
