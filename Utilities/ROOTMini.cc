#ifndef ROOTMINI_CC
#define ROOTMINI_CC

#include "JESTools.cc"

#include <TLorentzVector.h>

#include <vector>
#include <utility>

//For Minimizer
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class ROOTMini{
public:
  ROOTMini(JESTools *b_) {
    b = b_;
  };

  void SetMinimizer() {
    mini = ROOT::Math::Factory::CreateMinimizer("TMinuit");
    func = ROOT::Math::Functor(&CalcP,4);
    mini->SetPrintLevel(0);
    mini->SetStrategy(2);
    mini->SetMaxFunctionCalls(400);
    mini->SetMaxIterations(400);
    mini->SetTolerance(0.01);
    mini->SetErrorDef(0.5);
    mini->SetFunction(func);
  }

  void SetLep(TLorentzVector Lepton_, TLorentzVector MET_) {
    Lepton = Lepton_;
    LVMET = MET_;
  }

  void SetJets(vector<TLorentzVector> Jets_) {
    Jets = Jets_;
    vector< pair<double, double> > ScaleLimits;
    for (unsigned ij = 0; ij < Jets.size(); ++ij) {
      TLorentzVector Jet = Jets.at(ij);
      pair<double, double> ThisLimits = b->CalcLimits(Jet.Eta(),Jet.Pt());
      for (unsigned is = 0; is < Jets.size(); ++is) {
        mini->SetLimitedVariable(is,Form("Scale_%i",is),1.0,0.01,ThisLimits.first,ThisLimits.second);
      }
    }
  }

  static double CalcP(const double * scales) {
    //Calculation of PScale
    double PScale = b->CalcPScales(Jets, scales);

    //Scaling Jets and MET
    TLorentzVector ScaledMET;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, scales, LVMET, ScaledMET);

    //Solving Neutrino
    vector<TLorentzVector> Neutrinos;
    double PNeutrino = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);

    if (PNeutrino < 0) return PNeutrino;

    //Calculation of P on hadronic side
    double PHad = b->CalcPHad(ScaledJets);

    //Calculation of P on leptonic side
    double PLep = 0;
    // Neutrino = TLorentzVector();
    for (unsigned in = 0; in < Neutrinos.size(); ++in) {
      TLorentzVector NeutrinoTemp =  Neutrinos.at(in);
      double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
      if (PLepTMassTemp > PLep) {
        PLep = PLepTMassTemp;
        // Neutrino = NeutrinoTemp;
      }
    }
    // So far the P is the higher the better

    //Summing all and make it negative
    double Prob = PScale * PHad * PLep * (-1.0);

    return Prob;
  }

  double MinimizeP() {
    mini->Minimize();
    MinimizedScales.clear();
    double Prob = 0;
    if (!(mini->Status())) {
      // MinimizedScalesArray = mini->X();
      for (unsigned is = 0; is < Jets.size(); ++is) {
        MinimizedScalesArray[is] = mini->X()[is];
        MinimizedScales.push_back(mini->X()[is]);
      }
      Prob = (-1 * (mini->MinValue()));
    }
    return Prob;
  }

  vector<double> GetScales() {
    return MinimizedScales;
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
  ROOT::Math::Minimizer* mini;
  ROOT::Math::Functor func;

  //Inputs
  static TLorentzVector Lepton, LVMET;
  static vector<TLorentzVector> Jets;

  //Intermediate
  // static TLorentzVector Neutrino;

  //Outputs
  double * MinimizedScalesArray;
  vector<double> MinimizedScales;
  vector<TLorentzVector> MinizedJets;
};

//Initialization of static variables;
JESTools * ROOTMini::b;
// ROOT::Math::Minimizer* ROOTMini::mini;
// ROOT::Math::Functor ROOTMini::func;
TLorentzVector ROOTMini::Lepton;
TLorentzVector ROOTMini::LVMET;
vector<TLorentzVector> ROOTMini::Jets;
// TLorentzVector ROOTMini::Neutrino;

#endif
