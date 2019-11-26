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

void GenMinimizeTest(int SampleType = 0, int irun = 1, int debug = 0) {
  cout << "start" <<endl;
  //SetUpUtilities
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename = "MinimizerTest";


  JESTools *b = new JESTools();
  TFile* PFile = new TFile("PFile/PFile.root");
  b->ReadJESPlots(PFile);

  ROOTMini *m = new ROOTMini(b);
  // m->SetMinimizer();

  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    a->AssignGenParticles();
    // Gen as Jets
    vector<int> GenOutSort = a->GenOutSort;
    vector<TLorentzVector> Jets;
    bool allfound = true;
    for (unsigned ig = 0; ig < GenOutSort.size(); ++ig) {
      if (GenOutSort[ig] == -1) {
        allfound = false;
        Jets.push_back(TLorentzVector());
      }
      else {
        Jets.push_back(a->GetGenParticleP4(GenOutSort[ig]));
      }
    }
    if (!allfound) {
      cout << "Set of Gen particles incomplete, skipping" <<endl;
      continue;
    }
    double scalestemp[4] = {1,1,1,1};
    double * scales = scalestemp;
    TLorentzVector Lepton = a->LVGenLep;
    TLorentzVector LVMET = a->LVGenNeu;


    vector< vector<int> > perms = b->MakePermutations(Jets.size());
    vector<int> perm = perms.at(0);
    double PScale = b->CalcPScales(Jets, scales);
    TLorentzVector ScaledMET;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, scales, a->LVGenNeu, ScaledMET);
    vector<TLorentzVector> Neutrinos;
    double PNeutrino = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);

    if (PNeutrino < 0) continue;
    double PHad = b->CalcPHad(ScaledJets);
    double PLep = 0;
    TLorentzVector Neutrino = TLorentzVector();
    for (unsigned in = 0; in < Neutrinos.size(); ++in) {
      TLorentzVector NeutrinoTemp =  Neutrinos.at(in);
      // double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
      double PLepTMassTemp = b->CalcPTMass((ScaledJets.at(3)+Lepton+NeutrinoTemp).M());
      if (PLepTMassTemp > PLep) {
        PLep = PLepTMassTemp;
        Neutrino = NeutrinoTemp;
      }
    }
    double Prob = PScale * PHad * PLep * (-1.0);

    //Outputs;
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

    // cout << Form("PLep1 = %f, PLep2 = %f", b->CalcPTMass((ScaledJets[3]+Lepton+Neutrinos[0]).M()), b->CalcPTMass((ScaledJets[3]+Lepton+Neutrinos[1]).M()) )<<endl;

    // Using the minimizer
    m->SetLep(Lepton, LVMET);
    // m->SetJets(Jets);
    // m->SetMinimizer();
    m->MinimizeP(Jets);
    // double* minscales = m->GetScalesArray();
    // m->SetDebug(1);
    // m->CalcP(scales);



  }



}
