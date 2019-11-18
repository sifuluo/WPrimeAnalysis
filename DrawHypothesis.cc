#include "Utilities/DrawTools.cc"

#include <TROOT.h>
#include <TString.h>
#include <TVector2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TFile.h>

#include<iostream>

using namespace std;

void DrawHypothesis(){
  TFile *fFL = new TFile("results/hypothesis_FL_Add.root");
  TFile *fLL = new TFile("results/hypothesis_LL_Add.root");
  TFile *fBG = new TFile("results/hypothesis_BG_Add.root");

  DrawTH1(fFL, fLL, fBG, "TopdR", "TopdR", true, true);
}
