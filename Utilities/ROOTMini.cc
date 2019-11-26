#ifndef ROOTMINI_CC
#define ROOTMINI_CC

#include "JESTools.cc"

#include <TROOT.h>
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
  }

  // void SetDebug(int debug_) {
  //   debug = debug_;
  // }

  static double CalcP(const double *scales) {
    FunctionCalls ++;
    //Calculation of PScale
    double PScale = b->CalcPScales(Jets, scales);

    //Scaling Jets and MET
    TLorentzVector ScaledMET;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, scales, LVMET, ScaledMET);

    //Solving Neutrino
    vector<TLorentzVector> Neutrinos;
    double PNeutrino = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);

    if (PNeutrino < 0) return (-1.0) * PNeutrino;

    //Calculation of P on hadronic side
    double PHad = b->CalcPHad(ScaledJets);

    //Calculation of P on leptonic side
    double PLep = 0;
    TLorentzVector Neutrino = TLorentzVector();
    for (unsigned in = 0; in < Neutrinos.size(); ++in) {
      TLorentzVector NeutrinoTemp =  Neutrinos.at(in);
      double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
      if (PLepTMassTemp > PLep) {
        PLep = PLepTMassTemp;
        Neutrino = NeutrinoTemp;
      }
    }
    // So far the P is the higher the better


    //Summing all and make it negative
    double Prob = PScale * PHad * PLep * (-1.0);

    //Outputs for debugging
    if (true) {
      cout << endl <<"------In Minimizer------" <<endl;
      cout << "------Before Scale------" <<endl;
      double lepw = (Lepton + LVMET).M();
      double lept = (Lepton + LVMET + Jets[3]).M();
      double hadw = (Jets[0]+Jets[1]).M();
      double hadt = (Jets[0]+Jets[1]+Jets[2]).M();
      double plepw = b->CalcPWMass(lepw);
      double plept = b->CalcPTMass(lept);
      double phadw = b->CalcPWMass(hadw);
      double phadt = b->CalcPTMass(hadt);
      cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
      cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PTotal = %f",plepw,plept,phadw,phadt, plept*phadw*phadt) <<endl;
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
      cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PTotal = %f",plepw,plept,phadw,phadt, plept*phadw*phadt) <<endl;
      cout << "------In Summary------" <<endl;
      cout << Form("PScale = %f, PHad = %f, PLep = %f, P = %f",PScale,PHad, PLep, Prob)<<endl;
      cout << Form("Scales are: %f, %f, %f, %f", scales[0], scales[1], scales[2], scales[3]) << endl;
    }

    return Prob;
  }

  double MinimizeP(vector<TLorentzVector> Jets_) {
    FunctionCalls = 0;
    //Set up minimizer
    mini->SetPrintLevel(0);
    mini->SetStrategy(1);
    mini->SetMaxFunctionCalls(1000);
    mini->SetMaxIterations(1000);
    mini->SetTolerance(0.01);
    mini->SetErrorDef(0.5);
    mini->SetFunction(func);

    //Setting Jets
    Jets = Jets_;
    vector< pair<double, double> > ScaleLimits;
    for (unsigned ij = 0; ij < Jets.size(); ++ij) {
      TLorentzVector Jet = Jets.at(ij);
      pair<double, double> ThisLimits = b->CalcLimits(Jet.Eta(),Jet.Pt());
      for (unsigned is = 0; is < Jets.size(); ++is) {
        mini->SetLimitedVariable(is,Form("Scale_%i",is),1.0,0.01,ThisLimits.first,ThisLimits.second);
      }
    }

    // cout << endl<< "Minimizing started" <<endl;
    mini->Minimize();
    MinimizedScales.clear();
    double Prob = 0;
    if (!(mini->Status())) {
      cout << endl<< "Successfully minimized" <<endl;
      // MinimizedScalesArray = mini->X();
      for (unsigned is = 0; is < Jets.size(); ++is) {
        MinimizedScalesArray[is] = mini->X()[is];
        MinimizedScales.push_back(mini->X()[is]);
      }
      Prob = (-1 * (mini->MinValue()));
    }
    else {
      cout << endl<< "Minimizing failed with code: " << mini->Status() <<endl;
      Prob = -1.5;
    }
    cout << "FunctionCalls = " << FunctionCalls <<endl;
    return Prob;
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

  vector<double> GetPs(TLorentzVector &Neutrino) {
    double PScale = b->CalcPScales(Jets, MinimizedScalesArray);

    TLorentzVector ScaledMET;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, MinimizedScalesArray, LVMET, ScaledMET);

    vector<TLorentzVector> Neutrinos;
    double radical = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);

    double PHad = b->CalcPHad(ScaledJets);
    double PLep = 0;
    for (unsigned in = 0; in < Neutrinos.size(); ++in) {
      TLorentzVector NeutrinoTemp =  Neutrinos.at(in);
      double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
      if (PLepTMassTemp > PLep) {
        PLep = PLepTMassTemp;
        Neutrino = NeutrinoTemp;
      }
    }

    vector<double> Probs{PScale, PHad, PLep, radical};

    return Probs;
  }

private:
  static JESTools *b; // Base tool

  //Minimizer components
  static ROOT::Math::Minimizer* mini;
  static ROOT::Math::Functor func;

  //Inputs
  static TLorentzVector Lepton, LVMET;
  static vector<TLorentzVector> Jets;

  //Intermediate
  // static TLorentzVector Neutrino;
  // static int debug;
  static int FunctionCalls;

  //Outputs
  double * MinimizedScalesArray;
  vector<double> MinimizedScales;
  vector<TLorentzVector> MinizedJets;
};

//Initialization of static variables;
JESTools * ROOTMini::b;
ROOT::Math::Minimizer* ROOTMini::mini = ROOT::Math::Factory::CreateMinimizer("TMinuit");
ROOT::Math::Functor ROOTMini::func = ROOT::Math::Functor(&CalcP,4);
TLorentzVector ROOTMini::Lepton;
TLorentzVector ROOTMini::LVMET;
vector<TLorentzVector> ROOTMini::Jets;
int ROOTMini::FunctionCalls;
// TLorentzVector ROOTMini::Neutrino;
// int ROOTMini::debug;

#endif
