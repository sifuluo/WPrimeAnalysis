#ifndef ROOTMINI_CC
#define ROOTMINI_CC

#include "JESTools.cc"

#include <TROOT.h>
#include <TH1.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TString.h>

#include <vector>
#include <utility>
#include <iostream>
#include <string>

//For Minimizer
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

class ROOTMini{
public:
  ROOTMini(JESTools *b_) {
    // cout << endl <<"Invoked Minimizer" <<endl;
    b = b_;
  };

  void SetLep(TLorentzVector Lepton_, TLorentzVector MET_) {
    Lepton = Lepton_;
    LVMET = MET_;
    GenNeu = TLorentzVector();
  }

  void SetGenNeu(TLorentzVector GenNeu_) {
    GenNeu = GenNeu_;
  }

  void SetDebug(int debug_) {
    debug = debug_;
  }

  void SetTempOutput(int TempOutput_) {
    TempOutput = TempOutput_;
  }

  static double MinimizePFunc(const double *scales) {
    FunctionCalls ++;
    // // Test funciton
    // const double x = scales[0];
    // const double y = scales[1];
    // const double tmp1 = y-x*x;
    // const double tmp2 = 1-x;
    // return 100*tmp1*tmp1+tmp2*tmp2 + scales[2]* scales[2] + scales[3] * scales[3];

    //Calculation of PScale
    // double PScale = b->CalcPScalesHist(Jets, scales);
    double PScale = b->CalcPScalesFunc(Jets, scales);
    //Scaling Jets and MET
    TLorentzVector ScaledMET;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, scales, LVMET, ScaledMET);
    //Solving Neutrino
    // vector<TLorentzVector> Neutrinos;
    // double PNeutrino = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);
    // InterPNeutrino = PNeutrino;
    // if (PNeutrino < 0) {
    //   return ((-1.0) * PNeutrino + 1);
    // }
    // //Calculation of P on leptonic side
    // TLorentzVector Neutrino = TLorentzVector();
    // double PLep = b->CalcPLep(ScaledJets.at(3), Lepton, Neutrinos, Neutrino);

    TLorentzVector Neutrino = TLorentzVector();
    double PLep = b->CalcPLep(ScaledJets.at(3), Lepton, ScaledMET, Neutrino);
    if (PLep < 0) {
      return ((-1.0) * PLep + 1);
    }

    //Calculation of P on hadronic side
    double PHad = b->CalcPHad(ScaledJets);


    // for (unsigned ineu = 0; ineu < Neutrinos.size(); ++ineu) {
    //   TLorentzVector NeutrinoTemp =  Neutrinos.at(ineu);
    //   double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
    //   if (PLepTMassTemp > PLep) {
    //     PLep = PLepTMassTemp;
    //     Neutrino = NeutrinoTemp;
    //   }
    // }
    // So far the P is the higher the better


    //Summing all and make it negative and plus 1
    double Prob = (PScale * PHad * PLep * (-1.0) + 1);
    //Intermediate outputs
    if (TempOutput) {
      InterScalesVector.clear();
      for (unsigned is = 0; is < Jets.size(); ++is) {
        InterScalesVector.push_back(scales[is]);
      }
      InterScaledJets = ScaledJets;
      InterScaledMET = ScaledMET;
      InterNeutrino = Neutrino;
      vector<double> tempp{PScale,PHad,PLep,Prob};
      InterProbs = tempp;
    }

    if (RedoOutput) {
      RedoParticles.clear();
      RedoParticles = ScaledJets;
      RedoParticles.push_back(Lepton);
      RedoParticles.push_back(Neutrino);
    }

    //Outputs for debugging
    if (debug) {
      cout << endl <<"------In Minimizer------" <<endl;
      // cout << "------Before Scale------" <<endl;
      double lepw = (Lepton + LVMET).M();
      double lept = (Lepton + LVMET + Jets[3]).M();
      double hadw = (Jets[0]+Jets[1]).M();
      double hadt = (Jets[0]+Jets[1]+Jets[2]).M();
      double plepw = b->CalcPWMass(lepw);
      double plept = b->CalcPTMass(lept);
      double phadw = b->CalcPWMass(hadw);
      double phadt = b->CalcPTMass(hadt);
      // cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
      // cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PHad = %f, PLep = %f",plepw,plept,phadw,phadt, phadw*phadt, plepw*plept) <<endl;
      cout << "------After Scale------" <<endl;
      lepw = (Lepton + Neutrino).M();
      lept = (Lepton + Neutrino + ScaledJets[3]).M();
      hadw = (ScaledJets[0]+ScaledJets[1]).M();
      hadt = (ScaledJets[0]+ScaledJets[1]+ScaledJets[2]).M();
      plepw = b->CalcPWMass(lepw);
      plept = b->CalcPTMass(lept);
      phadw = b->CalcPWMass(hadw);
      phadt = b->CalcPTMass(hadt);
      cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
      cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PHad = %f, PLep = %f",plepw,plept,phadw,phadt, phadw*phadt, plepw*plept) <<endl;
      cout << "------In Summary------" <<endl;
      cout << Form("Scales are: %f, %f, %f, %f", scales[0], scales[1], scales[2], scales[3]) << endl;
      cout << Form("PScale = %f, PHad = %f, PLep = %f, P = %f, PCriteria = %f",PScale,PHad, PLep, PScale * PHad * PLep, Prob)<<endl;
    }

    return Prob;
  }

  double MinimizeP(vector<TLorentzVector> Jets_) {
    FunctionCalls = 0;
    //Set up minimizer
    func = ROOT::Math::Functor(&MinimizePFunc,4);
    if (Jets_.size() != 4) return -1;
    mini->SetPrintLevel(0);
    if (debug) {
      mini->SetPrintLevel(3);
    }
    mini->SetStrategy(3);
    mini->SetMaxFunctionCalls(100000);
    mini->SetMaxIterations(10000);
    mini->SetTolerance(0.01);
    mini->SetErrorDef(0.5);
    mini->SetFunction(func);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    // SetDebug(0);
    // Setting Jets
    Jets = Jets_;
    vector< pair<double, double> > ScaleLimits;
    for (unsigned ij = 0; ij < Jets.size(); ++ij) {
      TLorentzVector Jet = Jets.at(ij);
      if (Jet == TLorentzVector()) return -2;
      pair<double, double> ThisLimits = b->CalcLimitsFunc(Jet.Eta(),Jet.Pt());
      mini->SetLimitedVariable(ij,Form("Scale_%i",ij),1.0,0.01,ThisLimits.first,ThisLimits.second);
    }
    if (debug) cout << endl<< "Before Minimizing" <<endl;
    // cout << endl<< "Minimizing started" <<endl;
    mini->Minimize();
    if (debug) cout << endl<< "Minimizer Successfully Ran" <<endl;
    MinimizedScales.clear();
    double Prob = 0;
    if (!(mini->Status())) {
      if (debug) {
        cout << endl<< "Successfully minimized" <<endl;
      }
      for (unsigned is = 0; is < Jets.size(); ++is) {
        MinimizedScalesArray[is] = mini->X()[is];
        MinimizedScales.push_back(mini->X()[is]);
      }
      Prob = (1. - (mini->MinValue()));
    }
    else {
      if (debug) {
        cout << endl<< "Minimizing failed with code: " << mini->Status() <<endl;
      }
      Prob = -1;
    }
    if (debug) {
      cout << "FunctionCalls = " << FunctionCalls <<endl;
    }
    return Prob;
  }

  static double MinimizePLepFunc(const double *scales) {
    FunctionCalls++;
    double PScale = b->CalcPScalesFunc(LepJets, scales);
    // double PScale = b->CalcPScaleFunc(JetLepB.Eta(),JetLepB.Pt(),scales[0],0);
    TLorentzVector ScaledMET;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(LepJets,scales,LVMET,ScaledMET);
    TLorentzVector Neutrino = TLorentzVector();
    double PLep;
    if (GenNeu == TLorentzVector()) PLep = b->CalcPLep(ScaledJets.at(0), Lepton, ScaledMET, Neutrino);
    else PLep = b->CalcPLep(ScaledJets.at(0), Lepton, ScaledMET, GenNeu, Neutrino);
    if (PLep < 0) {
      return ((-1.0) * PLep +1);
    }
    double Prob = (PScale * PLep * (-1.0) + 1);
    return Prob;

  }

  double MinimizePLep(TLorentzVector LVLepB, double &LepBScale) {
    FunctionCalls = 0;
    func = ROOT::Math::Functor(&MinimizePLepFunc,1);
    mini->SetPrintLevel(0);
    if (debug) {
      mini->SetPrintLevel(3);
    }
    mini->SetStrategy(3);
    mini->SetMaxFunctionCalls(100000);
    mini->SetMaxIterations(10000);
    mini->SetTolerance(0.01);
    mini->SetErrorDef(0.5);
    mini->SetFunction(func);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

    LepJets.clear();
    LepJets.push_back(LVLepB);

    pair<double, double> ThisLimits = b->CalcLimitsFunc(LVLepB.Eta(),LVLepB.Pt());
    mini->SetLimitedVariable(0,"Scale_LepTB",1.0,0.01,ThisLimits.first,ThisLimits.second);

    mini->Minimize();

    double Prob = 0;
    if (!(mini->Status())) {
      Prob = (1. - (mini->MinValue()));
      LepBScale = mini->X()[0];
    }
    else Prob = -1;

    return Prob;

  }
  // Unfinished outputs
  vector<double> GetInterScalesVector() {
    return InterScalesVector;
  }

  // Outputs
  vector<double> GetScales() {
    return MinimizedScales;
  }

  double* GetScalesArray() {
    return MinimizedScalesArray;
  }

  vector<TLorentzVector> GetJets(TLorentzVector& met_) {
    return b->ScaleJets(Jets, MinimizedScalesArray, LVMET, met_ );
  }

  // vector<TLorentzVector> MinimizedParticles() {
  //   RedoOutput = 1;
  //   MinimizePFunc(MinimizedScalesArray);
  //   RedoOutput = 0;
  //   return RedoParticles;
  // }

  vector<double> ReCalcP(vector<TLorentzVector> &Particles) {
    double PScale = b->CalcPScalesFunc(Jets, MinimizedScalesArray);
    TLorentzVector ScaledMET, Neutrino;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, MinimizedScalesArray, LVMET, ScaledMET);

    vector<TLorentzVector> Neutrinos;
    double radical = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);

    double PHad = b->CalcPHad(ScaledJets);
    double PLep = 0;
    for (unsigned ineu = 0; ineu < Neutrinos.size(); ++ineu) {
      TLorentzVector NeutrinoTemp =  Neutrinos.at(ineu);
      double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
      if (PLepTMassTemp > PLep) {
        PLep = PLepTMassTemp;
        Neutrino = NeutrinoTemp;
      }
    }

    vector<double> Probs{PScale, PHad, PLep, radical};
    vector<TLorentzVector> parts = ScaledJets;
    parts.push_back(Lepton);
    parts.push_back(Neutrino);

    Particles = parts;
    //parts are LFJet1, LFJet2, HadB, LepB, Lepton, Neutrino
    return Probs;
  }

  // Unfinished outputs
  static vector<double> InterScalesVector;
  static vector<TLorentzVector> InterScaledJets;
  static TLorentzVector InterScaledMET;
  static TLorentzVector InterNeutrino;
  static vector<double> InterProbs;
  static double InterPNeutrino;

  //Outputs
  double * MinimizedScalesArray = new double[4];
  vector<double> MinimizedScales;
  static vector<TLorentzVector> RedoParticles;

private:
  static JESTools *b; // Base tool

  //Minimizer components
  static ROOT::Math::Minimizer* mini;
  static ROOT::Math::Functor func;

  //Inputs
  static TLorentzVector Lepton, LVMET, GenNeu;
  static vector<TLorentzVector> Jets;
  static vector<TLorentzVector> LepJets;

  //Intermediate
  // static TLorentzVector Neutrino;
  static int debug;
  static int TempOutput;
  static int RedoOutput;
  static int FunctionCalls;
};

//Initialization of static variables;
JESTools * ROOTMini::b;
ROOT::Math::Minimizer* ROOTMini::mini = ROOT::Math::Factory::CreateMinimizer("TMinuit");
ROOT::Math::Functor ROOTMini::func;
TLorentzVector ROOTMini::Lepton;
TLorentzVector ROOTMini::LVMET;
TLorentzVector ROOTMini::GenNeu;
vector<TLorentzVector> ROOTMini::Jets;
vector<TLorentzVector> ROOTMini::LepJets;
int ROOTMini::FunctionCalls;

//Intermediate outputs
vector<double>  ROOTMini:: InterScalesVector;
vector<TLorentzVector> ROOTMini::InterScaledJets;
TLorentzVector ROOTMini::InterScaledMET;
TLorentzVector ROOTMini::InterNeutrino;
vector<double> ROOTMini::InterProbs;
double ROOTMini::InterPNeutrino;

vector<TLorentzVector> ROOTMini::RedoParticles;

// TLorentzVector ROOTMini::Neutrino;
int ROOTMini::debug = 0;
int ROOTMini::TempOutput = 0;
int ROOTMini::RedoOutput = 0;

#endif
