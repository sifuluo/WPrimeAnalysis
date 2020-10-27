#ifndef HYPOTHESIS_CC
#define HYPOTHESIS_CC
#include <TLorentzVector.h>
#include <vector>

class Hypothesis {
public:
  Hypothesis() { // TLorentzVectors should be stored inside of this class, pointers are only for branch use
    // VToP();
  }
  TLorentzVector WPB, HadB, LF0, LF1, LepB, Lep, Neu, HadW, HadT, LepW, LepT, WP;
  TLorentzVector *pWPB, *pHadB, *pLF0, *pLF1, *pLepB, *pLep, *pNeu, *pHadW, *pHadT, *pLepW, *pLepT, *pWP;
  // void PtoV() {
  //   WPB = *pWPB;
  //   HadB = *pHadB;
  //   LF0 = *pLF0;
  //   LF1 = *pLF1;
  //   LepB = *pLepB;
  //   Lep = *pLep;
  //   Neu = *pNeu;
  //   HadW = *pHadW;
  //   HadT = *pHadT;
  //   LepW = *pLepW;
  //   LepT = *pLepT;
  //   WP = *pWP;
  // }
  void VToP() {
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


};

#endif
