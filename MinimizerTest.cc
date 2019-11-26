#ifndef MINIMIZERTEST_CC
#define MINIMIZERTEST_CC

#include "Utilities/Analyzer.cc"
#include "Utilities/JESTools.cc"
#include "Utilities/ROOTMini.cc"

// Delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include <TROOT.h>
#include <TH1.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>

void MinimizerTest(int SampleType = 0, int irun = 1, int debug = 0) {
  cout << "start" <<endl;
  //SetUpUtilities
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MinimizerTest";


  JESTools *b = new JESTools();
  TFile* PFile = new TFile("PFile/PFile.root");
  b->ReadJESPlots(PFile);

  ROOTMini *m = new ROOTMini(b);
  m->SetMinimizer();

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    // if (a->Jets.size() < 5 || a->GenW.size() != 2 || a->LVLeptons.size() != 1) {
    //   cout <<"Did not pass Preselection" <<endl;
    //   continue;
    // }
    a->AssignGenParticles();
    m->SetLep(a->LVGenLep,a->LVGenNeu);

    // Gen as Jets
    vector<int> GenOutSort = a->GenOutSort;
    vector<TLorentzVector> ScaledJets;
    bool allfound = true;
    for (unsigned ig = 0; ig < GenOutSort.size(); ++ig) {
      if (GenOutSort[ig] == -1) {
        allfound = false;
        ScaledJets.push_back(TLorentzVector());
      }
      else {
        ScaledJets.push_back(a->GetGenParticleP4(GenOutSort[ig]));
      }
    }
    if (!allfound) {
      cout << "Set of Gen particles incomplete, skipping" <<endl;
      continue;
    }

    m->SetJets(ScaledJets);
    double scales[4] = {1,1,1,1};
    vector<bool> BTags{false,false,true,true,true};
    double lepw = (a->LVGenLep + a->LVGenNeu).M();
    double lept = (a->LVGenLep + a->LVGenNeu + ScaledJets[3]).M();
    double hadw = (ScaledJets[0]+ScaledJets[1]).M();
    double hadt = (ScaledJets[0]+ScaledJets[1]+ScaledJets[2]).M();
    double plepw = b->CalcPWMass(lepw);
    double plept = b->CalcPTMass(lept);
    double phadw = b->CalcPWMass(hadw);
    double phadt = b->CalcPTMass(hadt);

    cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
    cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PTotal = %f",plepw,plept,phadw,phadt, plept*phadw*phadt) <<endl;

    vector< vector<int> > perms = b->MakePermutations(ScaledJets.size());
    cout << "GenPar size = " <<ScaledJets.size()  << " perm size = " << perms.size() <<endl;

    vector<int> perm = perms.at(0);
    double PFlavor = b->CalcPFlavor(perm, BTags);

    vector<TLorentzVector> Neutrinos;
    double PNeutrino = b->SolveNeutrinos(a->LVGenLep, a->LVGenNeu, Neutrinos);
    cout << "This Perm is " << Form("%i,%i,%i,%i",perm[0],perm[1],perm[2],perm[3]) << " FLavor P = " << PFlavor;
    cout << "PNeutrino = " << PNeutrino <<endl;
    double lepw1 = (a->LVGenLep + Neutrinos[0]).M();
    double lept1 = (a->LVGenLep + Neutrinos[0] + ScaledJets[3]).M();
    double lepw2 = (a->LVGenLep + Neutrinos[1]).M();
    double lept2 = (a->LVGenLep + Neutrinos[1] + ScaledJets[3]).M();
    cout << Form("lepw1 mass = %f, plepw1 = %f, lept1 mass = %f, plept1 = %f, neuDR = %f", lepw1, b->CalcPWMass(lepw1), lept1, b->CalcPTMass(lept1), Neutrinos[0].DeltaR(a->LVGenNeu))<<endl;
    cout << Form("lepw2 mass = %f, plepw2 = %f, lept2 mass = %f, plept2 = %f, neuDR = %f", lepw2, b->CalcPWMass(lepw2), lept2, b->CalcPTMass(lept2), Neutrinos[1].DeltaR(a->LVGenNeu))<<endl;
    double phad = b->CalcPHad(ScaledJets);
    TLorentzVector Neutrino;
    if (b->CalcPTMass(lept1) > b->CalcPTMass(lept2)) Neutrino = Neutrinos[0];
    else Neutrino = Neutrinos[1];
    double plep = b->CalcPLep(ScaledJets[3], a->LVGenLep, Neutrino);
    cout << Form("The Prob after solve neutrino is: PHad = %f, PLep = %f, P = %f", phad, plep, phad * plep) << endl;

    double ToolP = m->CalcP(scales);
    cout << Form("JESTool give the probability as: %f", ToolP) <<endl;

    // vector<TLorentzVector> LVPerm = b->GetPermutationLV(perm,ScaledJets);
    // m->SetJets(LVPerm);
    // double Prob = m->MinimizeP();
    // vector<double> Scales = m->GetScales();
    // TLorentzVector Neutrino;
    // vector<double> Ps = m->GetPs(Neutrino);
    //
    // cout << "This Perm P = " << Prob <<endl;
    // cout << "PScale = " << Ps[0] << " PHad = " << Ps[1] << " PLep = " << Ps[2] << " radical = " << Ps[3] <<endl;


    // for (unsigned ip = 0; ip < perms.size(); ++ip) {
    //   vector<int> perm = perms.at(ip);
    //   double PFlavor = b->CalcPFlavor(perm, BTags);
    //   vector<TLorentzVector> LVPerm = b->GetPermutationLV(perm,ScaledJets);
    //   m->SetJets(LVPerm);
    //   double Prob = m->MinimizeP();
    //   vector<double> Scales = m->GetScales();
    //   TLorentzVector Neutrino;
    //   vector<double> Ps = m->GetPs(Neutrino);
    //   cout << "This Perm is " << Form("%i,%i,%i,%i",perm[0],perm[1],perm[2],perm[3]) <<endl;
    //   cout << "This Perm P = " << Prob <<endl;
    //   cout << "PScale = " << Ps[0] << " PHad = " << Ps[1] << " PLep = " << Ps[2] << " radical = " << Ps[3] <<endl;
    //   break;

      // ProbDis->Fill(Prob);
      // if (Prob > BestP) BestP = Prob;
    // }
    // ProbDis->Fill(BestP);
    // if (BestP > 0) EventCounter->Fill(2);
    // else EventCounter->Fill(3);
  }

  // a->SaveOutput();
}



#endif
