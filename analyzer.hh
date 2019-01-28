#ifndef ANALYZER_HH
#define ANALYZER_HH

#include "tools.cc"
#include "optimizer.cc"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

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
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>

//For Minimizer
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

//For finding files
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>

class Analyzer{
public:
  Analyzer(vector<TString> basepaths, vector<TString> inputfolders,int irun, double pt = 30);
  void SetStartEntry(int starentry);
  void ProcessEntries(int nprocess);
  void SetEndEntry(int endentry);
  void SetOutput(TString outputfolder, TString outputname);
  void SetJetPtThreshold(double pt);
  Long64_t GetEntries();
  Long64_t GetStartEntry();
  Long64_t GetEndEntry();
  int ReadEvent(Int_t ievt, bool debug = false);
  void GetInfos();
  void SortInitialize();
  void SortGenParticle();
  void AssignGenParticles();
  void SaveOutput();

  TString GenParticleInfo(int iGen);
  void PrintGenParticles(string type = "");
  int BeRealMother(int igen);
  int BeRealDaughter(int igen);
  vector<int> WeHadAMother(vector<int> &list, int &mother);
  vector<int> WeHadADaughter(vector<int> &list, int &daughter);
  vector<int> MyDaughters(int &mother);
  vector<int> MyMothers(int &daughter);
  vector<int> WayToBeHadron(int igen);
  int ToBeHadron(int igen);

  GenParticle* GetGenParticle(int igen);
  Jet* GetJet(int ijet);
  Muon* GetMuon(int imuon);
  Electron* GetElectron(int ie);

  TLorentzVector RecoParticle(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, double target, int &p1, int &p2);
  map<int,int> JetMatch(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, double &maxDeltaR, bool cut = true);
  void GetProbability(bool ForceRecreate = false);
  void CreateProbability(TString pfilename);
  double Optimize();
  double Optimize(vector<int> &BestPerm, vector<TLorentzVector> &scaledjets, TLorentzVector &neutrino, TLorentzVector &WPrime, double &pjes, double &pmass, double &ptype);
  int OptiJetMatch(double &dR);

  // double CalculateProbability(vector<int> jets);
  // double CalculateTypeProbability(vector<int> jets);
  // double CalculateJESProbability(const double *scales);
  // double CalculateMassProbability(vector<TLorentzVector> lvjets);
  // double GetMaxProbability(vector<double> &scales);

  //Sorting Particles
  vector<GenParticle*> GenParticles;
  vector<Jet*> AllJets, Jets, BJets, NBJets;
  vector<Electron*> Electrons;
  vector<Muon*> Muons;
  MissingET* MET;
  vector<int> GenD, GenU, GenS, GenC, GenB, GenT, GenE, GenNuE, GenMu, GenNuMu, GenTau, GenG, GenGamma, GenW, GenWP, GenOut, GenOutSort;
  vector<TLorentzVector> LVAllJets, LVJets, LVBJets, LVNBJets, LVElectrons, LVMuons, LVLeptons, LVSoftLep, LVGenJets; //These are all reconstructed TLorentzVectors.
  TLorentzVector LVMET;

  //Assigning Particles
  int WP, GenWPB, GenWPT, GenWPTB, GenWPTW, GenOtW, GenOtT, GenOtB, GenLep, GenNeu;
  TLorentzVector LVWP, LVGenWPB, LVGenWPT, LVGenWPTB, LVGenWPTW, LVGenOtW, LVGenOtT, LVGenOtB, LVGenLep, LVGenNeu, LVGenWPTWJ1, LVGenWPTWJ2;
  vector<TLorentzVector> LVGenWPTWJ, LVGenOutSort, LVJetSort;
  vector<int> GenWPTWJ;
  map<int,int> AllJetMatchMap, OutJetMatchMap, OptiJetMatchMap;
  double AllJetMatchMaxDeltaR, OutJetMatchDeltaR;

  //The probability calculations
  TFile* PFile;
  vector<double> etabins, scales;
  vector< vector<double> > ptbins;
  vector<int> npt;
  TH1F *pGenTMass, *pGenWMass, *pGenWPdPhi;
  vector<TH1F*> pJet;
  vector<bool> BTags;

  //Optimizer Results
  vector<int> OptiBestPerm;
  vector<TLorentzVector> OptiScaledjets;
  TLorentzVector OptiNeutrino, OptiWPrime;
  double OptiP;

  ofstream logfile;



private:
  ExRootTreeReader *treeReader;
  TClonesArray *branchGen, *branchJet, *branchElectron, *branchMuon, *branchMPT;
  TFile *ofile;
  TChain* chain_;
  TString OutputName, OutputLogName;
  Long64_t nEntries, iEntry, StartEntry, EndEntry;
  double JetPtThreshold;
  bool verbose;
  int JetPtCutOff, JetNumberCutoff, irun;
};



#endif
