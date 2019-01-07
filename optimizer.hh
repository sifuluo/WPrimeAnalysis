#ifndef OPTIMIZER_HH
#define OPTIMIZER_HH

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

//Root Related
#include <TROOT.h>
#include <TString.h>
#include <TVector2.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TError.h>
#include <TFitResult.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>

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

class Optimizer{
public:
  Optimizer(vector<TLorentzVector> LVJets_, vector<bool> BTags_);
  void SetLepton(TLorentzVector lepton_);
  void SetMET(TLorentzVector MET_);
  void SetPFile(TFile* PFile_);
  void SetEtaBins(vector<double> etabins_);
  void GetPermutations();
  void Optimize();
  double GetBestP();
  vector<int> GetBestPerm();
  vector<double> GetBestScales();
  vector<TLorentzVector> GetBestLVJets();
  TLorentzVector GetBestNeutrino();
  TLorentzVector GetWPrime();

private:
  void Initialize();

  double GetpType(vector<int> perm_);
  static double SolveNeutrino(TLorentzVector ScaledMET);
  static double GetpJES(const double *scales);
  static double GetpMass(vector<TLorentzVector> lvjets);
  double OptimizeThisPerm(vector<double> &ThisScale);
  void GradientOptimize();
  void FindBestPerm();
  static vector<TLorentzVector> LVJets, OrderedJets, scaledjets;
  static vector<bool> BTags;
  static vector<double> etabins;
  static TFile* PFile;
  static vector<TH1F*> pJet;
  static vector<vector <int> > Permutations;
  static TLorentzVector LVLep, LVMET, LVNeu1, LVNeu2, LVNeu, BestNeutrino;
  static double LepWMass, BestP;
  static vector<int> BestPerm;
  static vector<double> BestScales;
  static ROOT::Math::Minimizer *mini;
  static TF1* stdfgaus, *TopMassDis, *WMassDis;
  static bool debug;
};


#endif
