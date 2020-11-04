#ifndef GENTOOLS_CC
#define GENTOOLS_CC

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TLorentzVector.h>

#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

class GenTools{
public:
  GenTools(TString filepath, TString filename, int SampleType_ = 2) {
    SampleType = SampleType_;
    filename = filepath+"/"+filename+"_";
    fFL = new TFile(filename+"FL.root");
    fLL = new TFile(filename+"LL.root");
    if (SampleType > 1) fBG = new TFile(filename+"BG.root");
    Normalization();
    doFlags.clear();
    doFlags.resize(PNames.size(),true);
  };

  void Normalization() {
    for (int isample = 0; isample <= SampleType; ++isample){
      TFile *ftemp;
      if (isample == 0) ftemp = fFL;
      else if (isample == 1) ftemp = fLL;
      else ftemp = fBG;
      TList *ls = ftemp->GetListOfKeys();
      for (int iplot = 0; iplot < ls->GetSize(); ++iplot) {
        // TH1F* hist = (TH1F*) ftemp->Get(PNames.at(iplot));
        TH1F* hist = (TH1F*) ftemp->Get(ls->At(iplot)->GetName());
        hist->Scale(1./hist->GetMaximum());
      }
    }
  }

  double CalcP(TString pname_, double value_, int isample) {
    TFile *ftemp;
    if (isample == 0) ftemp = fFL;
    else if (isample == 1) ftemp = fLL;
    else ftemp = fBG;
    TH1F* hist = (TH1F*) ftemp->Get(pname_);
    double p = hist->GetBinContent(hist->FindBin(value_));
    if (p>1) cout << endl<<"Histogram: "<< pname_ << " at "<< value_<<" is " << p <<endl;
    return p;
  }

  TH1F* GetHist(TString pname_, int isample) {
    TFile* ftemp;
    if (isample == 0) ftemp = fFL;
    else if (isample == 1) ftemp = fLL;
    else ftemp = fBG;
    TH1F* hist = (TH1F*) ftemp->Get(pname_);
    return hist;
  }

  // double doValue(vector<>)

  vector<bool> doFlags;

private:
  int SampleType;
  TFile *fFL, *fLL, *fBG;
  // vector<TH1F*> WPdR, WPdPhi, WPBPt, TopdR, TopdPhi, TopPt, TopHadPt, TopLepPt, TopPtDiff;
  vector<TString> PNames = {"WPdR", "WPdPhi", "WPBPt", "TopdR", "TopdPhi", "TopPt", "TopHadPt", "TopLepPt", "TopPtDiff"};

  // vector<vector<TH1F*> > GenPlots;

};

#endif
