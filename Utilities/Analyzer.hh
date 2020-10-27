#ifndef ANALYZER_HH
#define ANALYZER_HH
// Customized modules
#include "SearchFiles.cc"
#include "ProgressBar.cc"
#include "JetMatch.cc"
#include "GetPID.cc"
#include "ROOTMini.cc"
#include "JESTools.cc"
#include "GenTools.cc"
#include "Hypothesis.cc"

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
  Long64_t nEntries, iEntry, StartEntry, EndEntry;
  // 0 for TDual_FormerLeptonic; 1 for TDual_LatterLeptonic; 2 for TTBar
  ofstream logfile;

  // Initialize Event
  void ReadEvent(Int_t ievt);
  void CountEvent(TString lbl);
  set<TString> CounterLabels;
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
  vector<bool> GenOutBTags;
  TLorentzVector LVMET;

  // Prepare plots
  void AddPlot(TH1F* hh_);
  void AddPlot(TH2F* hh_);
  map<string, TH1F*> Plots1D;
  map<string, TH2F*> Plots2D;

  // Assigning Gen Particles
  int AssignGenParticles();
  int WP, GenWPB, GenHadT, GenHadB, GenHadW, GenLepT, GenLepB, GenLepW, GenLep, GenNeu;
  TLorentzVector LVGenWP, LVGenWPB, LVGenHadT, LVGenHadB, LVGenHadW, LVGenLepT, LVGenLepB, LVGenLepW, LVGenLep, LVGenNeu, LVGenLF0, LVGenLF1;
  vector<TLorentzVector> LVGenOutSort;
  vector<int> GenLFJet, GenOutSort;
  int GenWPT, GenWPTB, GenWPTW, GenOtT, GenOtB, GenOtW;
  TLorentzVector LVGenWPT, LVGenWPTB, LVGenWPTW, LVGenOtT, LVGenOtB, LVGenOtW;

  TLorentzVector GetGenParticleP4(int igen);

  GenTools *GT1, *GT2;
  void SetupGenTools();

  JESTools *JT;
  TFile* JESFile;
  void SetupJESTools();

  ROOTMini *RM;
  void SetupROOTMini();

  double JetMatchMaxDeltaR;
  void MatchJets();
  map<int, vector<int> > AdvJetMatchMap; // map of indices of LVOutPart and LVJets
  map<int, int> JetMatchMap;
  void GetGenCorrectPerm();
  vector<int> GenCorrectPerm; // correct set of indices in LVOutPart/GenOutQuark in the order of hypothesis
  void GetRecoCorrectPerm();
  vector<int> RecoCorrectPerm;
  void GetRecoHypothesis();
  vector<TLorentzVector> RecoHypothesis;
  int MatchedHypo = 0;

  pair<double, vector<TLorentzVector> > SolveTTbar(vector<TLorentzVector> jets_, vector<bool> btags_, vector<int> &bestperm_);

  TFile *TreeFile;
  TTree* t;
  // Branches to store in the tree and Reco Branch containers
  TLorentzVector *m_GenWP, *m_GenWPB, *m_GenHadT, *m_GenHadB, *m_GenHadW, *m_GenLF0, *m_GenLF1, *m_GenLepT, *m_GenLepB, *m_GenLepW, *m_GenLep, *m_GenNeu;
  TLorentzVector LVRecoWPB, LVRecoHadB, LVRecoLF0, LVRecoLF1, LVRecoLepB, LVRecoLep, LVRecoNeu, LVRecoHadW, LVRecoHadT, LVRecoLepW, LVRecoLepT, LVRecoWP;
  TLorentzVector *m_RecoWPB, *m_RecoHadB, *m_RecoLF0, *m_RecoLF1, *m_RecoLepB, *m_RecoLep, *m_RecoNeu, m_RecoHadW, m_RecoHadT, m_RecoLepW, m_RecoLepT, m_RecoWP;
  Hypothesis Gen, Reco;
  // End of Branch declaration
  void Tree_Init(int SaveTreeLevel);
  void Tree_Reco();
  void Tree_Fill();
  void Tree_Save();

private:
  ExRootTreeReader *treeReader;
  TClonesArray *branchGen, *branchJet, *branchElectron, *branchMuon, *branchMPT, *branchGenJet;
  TFile *ofile;
  TChain* chain_;
  TString OutputName, OutputLogName, outputname, outputfolder;
  double JetPtThreshold, LepIso;
  bool verbose;
  int debug;
  int JetPtCutOff, JetNumberCutoff, irun;
};
#endif
