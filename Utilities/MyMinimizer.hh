#ifndef MYMINIMIZER_HH
#define MYMINIMIZER_HH

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TVector2.h>
#include <TString.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>

//For Minimizer
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class MyMinimizer{
public:
  MyMinimizer(vector<TLorentzVector> LVJets_, vector<bool> BTags_);
  void SetLepton(TLorentzVector lepton_);
  void SetMET(TLorentzVector met_);
  void SetPFile(TFile* PFile_);
  void SetEtaBins(vector<double> etabins_);
  void GetPermutations();
  void Optimize();
  double GetBestP();
  vector<int> GetBestPerm();

private:
  void Initialize();

  double GetpType(vector<int> perm_);
  static double SolveNeutrino(TLorentzVector scaledMET);
  static double GetpJES(const double  *scales);
  static double GetpMass(vector<TLorentzVector> lvjets);

}
