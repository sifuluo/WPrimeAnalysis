#ifndef CREATEPROBABILITY_CC
#define CREATEPROBABILITY_CC

#include "Analyzer.cc"
#include "ProbabilityBins.cc"

// ROOT
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

// std
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

void CreateProbability(int SampleType = 0, int irun = 0) {
  TString pfilename;
  double pt = 30;
  if (SampleType == 0) {
    pfilename = "TDual_FormerLeptonic";
  }
  else if (SampleType == 1) {
    pfilename = "TDual_LatterLeptonic";
  }
  else {
    pfilename = "Backgroud";
  }
  pfilename = pfilename + "_Probability.root";
  cout << "Creating Probability Histrograms" <<endl;
  Analyzer *a = new Analyzer(SampleType, irun, pt);
  TFile *PFile = new TFile(pfilename,"RECREATE");
  vector< vector<double> > ptbins;
  vector<int> npt;
  vector<double> etabins;

  ProbabilityBins(ptbins, npt, etabins);

}

#endif
