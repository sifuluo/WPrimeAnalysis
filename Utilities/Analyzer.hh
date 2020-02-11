#ifndef ANALYZER_HH
#define ANALYZER_HH
// Customized modules
#include "SearchFiles.cc"
#include "ProgressBar.cc"
#include "JetMatch.cc"
#include "GetPID.cc"
#include "JetMatch.cc"
// #include "Optimizer.cc"

// Delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

// ROOT
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

// std
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>

class Analyzer{
public:
  //Initialize Analysis
  Analyzer(int SampleType_, int irun_, double pt);
  void SetStartEntry(int startentry);
  void ProcessEntries(int nprocess);
  void SetEndEntry(int endentry);
  void SetOutput(TString outputfolder, TString outputname);
  void CDOut();
  void SetJetPtThreshold(double pt);
  Long64_t GetEntries();
  Long64_t GetStartEntry();
  Long64_t GetEndEntry();
  void DebugMode(int debug);
  void SaveOutput();
  int SampleType;
  // 0 for TDual_FormerLeptonic; 1 for TDual_LatterLeptonic; 2 for TTBar
  ofstream logfile;

  // Initialize Event
  void ReadEvent(Int_t ievt);
  void GetInfos();
  void GenParticleTypes();
  TString GenParticleInfo(int iGen);
  void PrintGenParticles(string type, bool verb);
  int BeRealMother(int igen);
  int BeRealDaughter(int igen);
  vector<int> WeHadAMother(vector<int> &list, int &mother);
  vector<int> WeHadADaughter(vector<int> &list, int &daughter);
  vector<int> MyDaughters(int &mother);
  vector<int> MyMothers(int &daughter);
  vector<int> WayToBeHadron(int igen);
  int ToBeHadron(int igen);
  int GenPass, RecoPass;
  vector<GenParticle*> GenParticles;
  vector<Jet*> AllJets, Jets, BJets, NBJets, GenJets;
  vector<Electron*> Electrons;
  vector<Muon*> Muons;
  MissingET* MET;
  vector<TLorentzVector> LVAllJets, LVJets, LVBJets, LVNBJets, LVElectrons, LVMuons, LVLeptons, LVSoftLep, LVOutPart, LVGenJets;
  vector<bool> BTags;
  vector<int> GenD, GenU, GenS, GenC, GenB, GenT, GenE, GenNuE, GenMu, GenNuMu, GenTau, GenG, GenGamma, GenW, GenWP, GenOutQuark, GenOutGluon;
  TLorentzVector LVMET;

  // Assigning Gen Particles
  int AssignGenParticles();
  int WP, GenWPB, GenHadT, GenHadB, GenHadW, GenLepT, GenLepB, GenLepW, GenLep, GenNeu;
  TLorentzVector LVWP, LVGenWPB, LVGenHadT, LVGenHadB, LVGenHadW, LVGenLepT, LVGenLepB, LVGenLepW, LVGenLep, LVGenNeu;
  vector<TLorentzVector> LVGenLFJet, LVGenOutSort;
  vector<int> GenLFJet, GenOutSort;
  int GenWPT, GenWPTB, GenWPTW, GenOtT, GenOtB, GenOtW;
  TLorentzVector LVGenWPT, LVGenWPTB, LVGenWPTW, LVGenOtT, LVGenOtB, LVGenOtW;

  TLorentzVector GetGenParticleP4(int igen);

  void MatchJets();
  map<int, vector<int> > JetMatchMap;
  double JetMatchMaxDeltaR;
  void GetRecoHypothesis();
  vector<TLorentzVector> RecoHypothesis;
  int MatchedHypo = 0;





private:
  ExRootTreeReader *treeReader;
  TClonesArray *branchGen, *branchJet, *branchElectron, *branchMuon, *branchMPT, *branchGenJet;
  TFile *ofile;
  TChain* chain_;
  TString OutputName, OutputLogName;
  Long64_t nEntries, iEntry, StartEntry, EndEntry;
  double JetPtThreshold, AlgodR;
  bool verbose;
  int debug;
  int JetPtCutOff, JetNumberCutoff, irun;
};
#endif
