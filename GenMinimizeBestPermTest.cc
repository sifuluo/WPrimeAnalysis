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
  m->SetDebug(1);
  bool output_ = true;

  TH1F* GenHadWMass = new TH1F("GenHadWMass","GenHadWMass", 300,50., 110.);
  TH1F* GenLepWMass = new TH1F("GenLepWMass","GenLepWMass", 300,50., 110.);
  TH1F* GenHadTMass = new TH1F("GenHadTMass","GenHadTMass", 600,110., 230.);
  TH1F* GenLepTMass = new TH1F("GenLepTMass","GenLepTMass", 600,110., 230.);
  TH1F* GenPDis     = new TH1F("GenPDis"    ,"GenP distribution", 100,0,1);

  TH1F* SolLepTMass = new TH1F("SolLepTMass","SolLepTMass", 600,110., 230.);
  TH1F* SolPDis     = new TH1F("SolPDis"    ,"SolP distribution", 100,0,1);

  TH1F* MinHadWMass = new TH1F("MinHadWMass","MinHadWMass", 300,50., 110.);
  TH1F* MinLepWMass = new TH1F("MinLepWMass","MinLepWMass", 300,50., 110.);
  TH1F* MinHadTMass = new TH1F("MinHadTMass","MinHadTMass", 600,110., 230.);
  TH1F* MinLepTMass = new TH1F("MinLepTMass","MinLepTMass", 600,110., 230.);
  TH1F* MinPDis     = new TH1F("MinPDis"    ,"MinP distribution", 100,0,1);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    if (a->AssignGenParticles() == -1) continue;
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
    TLorentzVector Lepton = a->LVGenLep;
    TLorentzVector LVMET = a->LVGenNeu;

    vector< vector<int> > perms = b->MakePermutations(Jets.size());
    vector<int> perm = perms.at(0);

    // Showing results before minimization
    double scalestemp[4] = {1,1,1,1};
    double * scales = scalestemp;
    TLorentzVector ScaledMET,Neutrino;
    vector<TLorentzVector> ScaledJets = b->ScaleJets(Jets, scales, a->LVGenNeu, ScaledMET);
    vector<TLorentzVector> Neutrinos;
    double PNeutrino = b->SolveNeutrinos(Lepton, ScaledMET, Neutrinos);
    double lepw,lept,hadw,hadt,plepw,plept,phadw,phadt,PScale,PHad,PLep,Prob;

    if (PNeutrino < 0) {
      cout << "Gen PNeutrino = " << PNeutrino <<endl;
    }
    else {
      PScale = b->CalcPScalesFunc(Jets, scales);
      PHad = b->CalcPHad(ScaledJets);
      PLep = 0;
      Neutrino = TLorentzVector();
      for (unsigned in = 0; in < Neutrinos.size(); ++in) {
        TLorentzVector NeutrinoTemp =  Neutrinos.at(in);
        // double PLepTMassTemp = b->CalcPLep(ScaledJets.at(3), Lepton, NeutrinoTemp);
        double PLepTMassTemp = b->CalcPTMass((ScaledJets.at(3)+Lepton+NeutrinoTemp).M());
        if (PLepTMassTemp > PLep) {
          PLep = PLepTMassTemp;
          Neutrino = NeutrinoTemp;
        }
      }
      Prob = PScale * PHad * PLep * (1.0);

      //Outputs;
      lepw = (Lepton + LVMET).M();
      lept = (Lepton + LVMET + Jets[3]).M();
      hadw = (Jets[0]+Jets[1]).M();
      hadt = (Jets[0]+Jets[1]+Jets[2]).M();
      plepw = b->CalcPWMass(lepw);
      plept = b->CalcPTMass(lept);
      phadw = b->CalcPWMass(hadw);
      phadt = b->CalcPTMass(hadt);
      GenHadWMass->Fill(hadw);
      GenLepWMass->Fill(lepw);
      GenHadTMass->Fill(hadt);
      GenLepTMass->Fill(lept);
      GenPDis->Fill(plepw*plept*phadw*phadt);


      if (output_) {
        cout << endl<<endl;
        cout << "----Before Minimizer-----" <<endl;
        cout << "------Before Scale------" <<endl;
        cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
        cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PHad = %f, PLep = %f, PTotal = %f",plepw,plept,phadw,phadt, phadw*phadt, plepw*plept,plepw*plept*phadw*phadt) <<endl;
      }
      lepw = (Lepton + Neutrino).M();
      lept = (Lepton + Neutrino + ScaledJets[3]).M();
      hadw = (ScaledJets[0]+ScaledJets[1]).M();
      hadt = (ScaledJets[0]+ScaledJets[1]+ScaledJets[2]).M();
      plepw = b->CalcPWMass(lepw);
      plept = b->CalcPTMass(lept);
      phadw = b->CalcPWMass(hadw);
      phadt = b->CalcPTMass(hadt);
      SolLepTMass->Fill(lept);
      SolPDis->Fill(plepw*plept*phadw*phadt);


      if (output_) {
        cout << "------After Scale------" <<endl;
        cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
        cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PHad = %f, PLep = %f, PTotal = %f",plepw,plept,phadw,phadt, phadw*phadt, plepw*plept,plepw*plept*phadw*phadt) <<endl;
        cout << "------In Summary------" <<endl;
        cout << Form("PScale = %f, PHad = %f, PLep = %f, P = %f",PScale,PHad, PLep, Prob)<<endl;
        cout << Form("Scales are: %f, %f, %f, %f", scales[0], scales[1], scales[2], scales[3]) << endl;
      }
    }



    // cout << Form("PLep1 = %f, PLep2 = %f", b->CalcPTMass((ScaledJets[3]+Lepton+Neutrinos[0]).M()), b->CalcPTMass((ScaledJets[3]+Lepton+Neutrinos[1]).M()) )<<endl;

    // Using the minimizer
    m->SetLep(Lepton, LVMET);
    double miniP = m->MinimizeP(Jets);
    double tempminscales[4] = {m->InterScalesArray[0],m->InterScalesArray[1],m->InterScalesArray[2],m->InterScalesArray[3]};
    scales = tempminscales;
    vector<double> miniScales = m->MinimizedScales;
    TLorentzVector miniNeutrino;
    vector<double> miniPs = m->GetPs(miniNeutrino);
    // cout << Form("Minscales: %f, %f, %f,%f",scales[0],scales[1],scales[2],scales[3] )<<endl;
    //After Minimizer


    PNeutrino = m->InterPNeutrino;

    if (PNeutrino < 0) {
      cout << "Min PNeutrino = " << PNeutrino <<endl;
    }
    else {
      PScale = m->InterProbs[0];
      double InterPHad = m->InterProbs[1];
      double InterPLep = m->InterProbs[2];
      double InterProb = m->InterProbs[3];
      ScaledJets = m->InterScaledJets;
      Neutrino = m->InterNeutrino;
      PHad = b->CalcPHad(ScaledJets);
      PLep = b->CalcPTMass((ScaledJets.at(3)+Lepton+Neutrino).M());

      Prob = PScale * PHad * PLep * (1.0);
      lepw = (Lepton + LVMET).M();
      lept = (Lepton + LVMET + Jets[3]).M();
      hadw = (Jets[0]+Jets[1]).M();
      hadt = (Jets[0]+Jets[1]+Jets[2]).M();
      plepw = b->CalcPWMass(lepw);
      plept = b->CalcPTMass(lept);
      phadw = b->CalcPWMass(hadw);
      phadt = b->CalcPTMass(hadt);
      if (output_) {
        cout << endl<<endl;
        cout << "----After Minimizer-----" <<endl;
        cout << "------Before Scale------" <<endl;
        cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
        cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PHad = %f, PLep = %f, PTotal = %f",plepw,plept,phadw,phadt, phadw*phadt, plepw*plept,plepw*plept*phadw*phadt) <<endl;
      }
      lepw = (Lepton + Neutrino).M();
      lept = (Lepton + Neutrino + ScaledJets[3]).M();
      hadw = (ScaledJets[0]+ScaledJets[1]).M();
      hadt = (ScaledJets[0]+ScaledJets[1]+ScaledJets[2]).M();
      plepw = b->CalcPWMass(lepw);
      plept = b->CalcPTMass(lept);
      phadw = b->CalcPWMass(hadw);
      phadt = b->CalcPTMass(hadt);
      MinHadWMass->Fill(hadw);
      MinLepWMass->Fill(lepw);
      MinHadTMass->Fill(hadt);
      MinLepTMass->Fill(lept);
      MinPDis->Fill(plepw*plept*phadw*phadt);

      if (output_) {
        cout << "------After Scale------" <<endl;
        cout << Form("LepWMass = %f, LepTMass = %f, HadWMass = %f, HadTMass = %f",lepw,lept,hadw,hadt) <<endl;
        cout << Form("PLepWMass = %f, PLepTMass = %f, PHadWMass = %f, PHadTMass = %f, PHad = %f, PLep = %f, PTotal = %f",plepw,plept,phadw,phadt, phadw*phadt, plepw*plept,plepw*plept*phadw*phadt) <<endl;
        cout << "------In Summary------" <<endl;
        cout << Form("Scales are: %f, %f, %f, %f", scales[0], scales[1], scales[2], scales[3]) << endl;
        cout << Form("PScale = %f, PHad = %f, PLep = %f, P = %f",PScale,PHad, PLep, Prob)<<endl;
        cout << "Minimizer gives:" <<endl;
        cout << Form("PScale = %f, PHad = %f, PLep = %f, P = %f",PScale, InterPHad, InterPLep, InterProb)<<endl;
        cout << "------Direct Mini-----" <<endl;
        cout << Form("MinValue = %f",miniP) <<endl;
        cout << Form("Scales are: %f, %f, %f, %f", miniScales[0], miniScales[1], miniScales[2], miniScales[3]) <<endl;
        cout << Form("PScale = %f, PHad = %f, PLep = %f, P = %f",miniPs[0],miniPs[1], miniPs[2], miniPs[0]*miniPs[1]*miniPs[2])<<endl;
      }
    }


  }

  a->SaveOutput();

}
