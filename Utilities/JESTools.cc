#ifndef JESTOOLS_CC
#define JESTOOLS_CC

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TString.h>
#include <TLorentzVector.h>

#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

class JESTools{
public:
  JESTools(){
    // cout << endl <<"Invoked JESTools" <<endl;
    // JESVector.clear();
    // TempiEta = 0;
    // TempiPt = 0;
  };

  // const double RNBMTag(0.01 / 0.99), RNBTag(1.0), RBMTag(0.3/0.7), RBTag(1.0);
  const double RNBMTag = 0.01 / 0.99;
  const double RNBTag = 1.0;
  const double RBMTag = 0.3/0.7;
  const double RBTag = 1.0;

  const vector<double> etabins{0., 1.3, 2.5, 3.0, 5.2};

  const vector<vector<double> > ptbins{
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130., 150.,180.,220., 260., 300.,350.,400.,500.,1000.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150., 180.,220.,260., 300.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150., 180.,220.,260.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.}
  };

  TF1* TopMassDis = new TF1("TBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,300.0);
  TF1* WMassDis = new TF1("WBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,200.0);

  TFile* JESFile;

  vector< vector<TH1F*> > JESVector;

  vector<TH1F*> AddVector;

  vector< vector<TF1*> > JESFuncVector;


  pair<int,int> TempiBin;
  int TempiEta, TempiPt;
  TH1F* TempHist;

  vector<double> EtaBins(){
    return etabins;
  };

  vector<double> PtBins(int i){
    return ptbins.at(i);
  };

  double EtaBinLow(int i){
    return etabins.at(i);
  };

  double EtaBinHigh(int i){
    return etabins.at(i+1);
  };

  double PtBinLow(int i, int j){
    return ptbins.at(i).at(j);
  };

  double PtBinHigh(int i, int j) {
    return ptbins.at(i).at(j+1);
  };

  int CalcEtaBin(double eta_) {
    int iEta = 0;
    for (unsigned ieta = 0; ieta < (etabins.size() -1); ++ ieta ) {
      if (fabs(eta_) < etabins.at(ieta+1) ) {
        iEta = ieta;
        break;
      }
      if (ieta == (etabins.size() -2)) iEta = ieta;
    }
    return iEta;
  }

  int CalcPtBin(double eta_, double pt_) {
    int iEta = CalcEtaBin(eta_);
    int iPt = 0;
    for (unsigned ipt = 0; ipt < (ptbins.at(iEta).size() -1); ++ ipt ) {
      if (pt_ < ptbins.at(iEta).at(ipt+1) ) {
        iPt = ipt;
        break;
      }
      if (ipt == (ptbins.at(iEta).size() -2)) iPt = ipt;
    }
    return iPt;
  }

  pair<int,int> CalcBins(double eta_, double pt_) {
    TempiEta = 0;
    for (unsigned ieta = 0; ieta < (etabins.size() -1); ++ ieta ) {
      if (fabs(eta_) < etabins.at(ieta+1) ) {
        TempiEta = ieta;
        break;
      }
      if (ieta == (etabins.size() -2)) TempiEta = ieta;
    }
    // cout <<"eta = " << eta_ << " iEta = " << TempiEta<<endl;
    TempiPt = 0;
    for (unsigned ipt = 0; ipt < (ptbins.at(TempiEta).size() -1); ++ ipt ) {
      if (pt_ < ptbins.at(TempiEta).at(ipt+1) ) {
        TempiPt = ipt;
        break;
      }
      if (ipt == (ptbins.at(TempiEta).size() -2)) TempiPt = ipt;
    }
    // cout <<"pt = " << pt_ << " iPt = " << TempiPt<<endl;
    TempiBin = pair<int,int>(TempiEta, TempiPt);
    return TempiBin;
  }

  int GetiEta(){
    return TempiEta;
  }

  int GetiPt(){
    return TempiPt;
  }

  TH1F* GetPlot() {
    TempHist = JESVector.at(TempiEta).at(TempiPt);
    return TempHist;
  }

  vector< vector<TH1F*> > MakeJESPlots() {
    vector<vector <TH1F*> > jes;
    jes.clear();
    for (unsigned ieta = 0; ieta < EtaBins().size()-1; ++ieta) {
      vector<TH1F*> jeseta;
      jeseta.clear();
      for (unsigned ipt = 0; ipt < PtBins(ieta).size() -1; ++ipt){
        TString sn = Form("eta%d_pt%d", ieta, ipt);
        TString st = Form("eta%.1fto%.1f_pt%dto%d;Pt_{Gen}/Pt_{Reco}",EtaBinLow(ieta),EtaBinHigh(ieta), int(PtBinLow(ieta,ipt)), int(PtBinHigh(ieta, ipt)) );
        jeseta.push_back(new TH1F(sn,st,600,0,6));
      }
      jes.push_back(jeseta);
    }
    JESVector = jes;
    return jes;
  }

  void FillJESPlot(double fill, int ieta_, int ipt_) {
    JESVector.at(ieta_).at(ipt_)->Fill(fill);
  }

  void FillJESPlot(double fill, double eta_, double pt_) {
    CalcBins(eta_, pt_);
    JESVector.at(TempiEta).at(TempiPt)->Fill(fill);
  }

  //Below is for Minimizer Interface
  vector< vector<TH1F*> > ReadJESPlots(TFile* f_) {
    JESFile = f_;
    SetUpMassFunctions();
    vector< vector<TH1F*> > jes;
    vector< vector<TF1*> > fjes;
    jes.clear();
    fjes.clear();
    for (unsigned ieta = 0; ieta < EtaBins().size()-1; ++ieta) {
      vector<TH1F*> jeseta;
      vector<TF1*> fjeseta;
      jeseta.clear();
      fjeseta.clear();
      for (unsigned ipt = 0; ipt < PtBins(ieta).size() -1; ++ipt){
        TString sn = Form("eta%d_pt%d", ieta, ipt);
        TH1F* h1 = (TH1F*)(f_->Get(sn));
        jeseta.push_back(h1); //Histogram might not be accessible after TFile being closed
        //Fitting histogram
        TString fsn = Form("f_eta%d_pt%d", ieta, ipt);
        TF1* f1 = new TF1(fsn,"[0]*TMath::BreitWigner(x,[1],[2])",0,2);
        f1->SetParameters(h1->GetMaximum()/4., 1.0, h1->GetStdDev());
        h1->Fit(fsn,"RMQ0","",0.,2.);
        fjeseta.push_back(f1);
      }
      jes.push_back(jeseta);
      fjes.push_back(fjeseta);
    }
    JESVector = jes;
    JESFuncVector = fjes;
    cout << "Finished Fitting" <<endl;
    return jes;
  }

  void SetUpMassFunctions() {
    TopMassDis = new TF1("TBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,300.0);
    WMassDis = new TF1("WBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,200.0);
    // TopMassDis->SetParameters(100.,172.7,1.32); // Gen-level
    TopMassDis->SetParameters(100.,171.186,26.76); // Reco-Level Leptonic Top
    // WMassDis->SetParameters(100,80.385,2.085);
    WMassDis->SetParameters(100,80.385,2.738);
    TopMassDis->SetParameter(0,100./TopMassDis->Eval(171.186)); // normalized it to peak at y = 1;
    WMassDis->SetParameter(0,100./WMassDis->Eval(80.385)); // normalized it to peak at y = 1;
  }

  vector< vector<int> > MakePermutations5(int jetsize) {
    vector< vector<int> > Permutations;
    Permutations.clear();
    for (int ihad1 = 0; ihad1 < jetsize - 1; ++ihad1) {
      for (int ihad2 = ihad1 + 1; ihad2 < jetsize; ++ihad2) {
        for (int ihadb = 0; ihadb < jetsize; ++ihadb) {
          if (ihadb == ihad1 || ihadb == ihad2) continue;
          for (int ilepb = 0; ilepb < jetsize; ++ ilepb) {
            if (ilepb == ihad1 || ilepb == ihad2 || ilepb == ihadb) continue;
            for (int iwpb = 0; iwpb < jetsize; ++ iwpb) {
              if (iwpb == ihad1 || iwpb == ihad2 || iwpb == ihadb || iwpb == ilepb) continue;
              vector<int> perm{ihad1,ihad2,ihadb,ilepb,iwpb};
              Permutations.push_back(perm);
            }
          }
        }
      }
    } // Permutations are all stored
    return Permutations;
  }

  vector< vector<int> > MakePermutations(int jetsize) {
    vector< vector<int> > Permutations;
    Permutations.clear();
    for (int ihad1 = 0; ihad1 < jetsize - 1; ++ihad1) {
      for (int ihad2 = ihad1 + 1; ihad2 < jetsize; ++ihad2) {
        for (int ihadb = 0; ihadb < jetsize; ++ihadb) {
          if (ihadb == ihad1 || ihadb == ihad2) continue;
          for (int ilepb = 0; ilepb < jetsize; ++ ilepb) {
            if (ilepb == ihad1 || ilepb == ihad2 || ilepb == ihadb) continue;
            vector<int> perm{ihad1,ihad2,ihadb,ilepb};
            Permutations.push_back(perm);
          }
        }
      }
    } // Permutations are all stored
    return Permutations;
  }

  vector<int> FindWPB(int jetsize, vector<int> perm) {
    vector<int> bperm;
    for (int ij = 0; ij < jetsize; ++ij) {
      if (find(perm.begin(),perm.end(),ij) == perm.end()) bperm.push_back(ij);
    }
    return bperm;
  }

  double CalcPWPB(TLorentzVector wpb, TLorentzVector wpt) {
    double dr = wpb.DeltaR(wpt);
    TH1F* hist = AddVector.at(0);
    double n = hist->GetBinContent(hist->FindBin(dr));
    double norm = hist->GetBinContent(hist->GetMaximumBin());
    double p = n / norm;
    return p;
  }

  double CalcPFlavor(vector<int> perm_, vector<bool> BTags_) {
    vector<bool> btags;
    for (unsigned it = 0; it < perm_.size(); ++it) {
      btags.push_back(BTags_.at(perm_.at(it)));
    }
    return CalcPFlavor(btags);
  }

  double CalcPFlavor(vector<bool> BTags_) {
    double pf = 1;
    if (BTags_.at(0)) pf *= RNBMTag; // Non-b-jet is tagged to be a b;
    else pf *= RNBTag; // Non-b-jet tagged non-b;
    if (BTags_.at(1)) pf *= RNBMTag; // Non-b-jet is tagged to be a b;
    else pf *= RNBTag; // Non-b-jet tagged non-b;
    if (BTags_.at(2)) pf *= RBTag; // b-jet tagged as a b;
    else pf *= RBMTag; // b-jet tagged to be a non-b;
    if (BTags_.at(3)) pf *= RBTag; // b-jet tagged as a b;
    else pf *= RBMTag; // b-jet tagged to be a non-b;
    if(BTags_.size() > 4){
      if (BTags_.at(4)) pf *= RBTag; // b-jet tagged as a b;
      else pf *= RBMTag; // b-jet tagged to be a non-b;
    }
    return pf;
  }

  vector<double> CalcPFlavorVector(vector<int> perm_, vector<bool> BTags_) {
    vector<bool> btags;
    for (unsigned it = 0; it < perm_.size(); ++it) {
      btags.push_back(BTags_.at(perm_.at(it)));
    }
    return CalcPFlavorVector(btags);
  }

  vector<double> CalcPFlavorVector(vector<bool> BTags_) {
    vector<double> out;
    if (BTags_.at(0)) out.push_back(RNBMTag); // Non-b-jet is tagged to be a b;
    else out.push_back(RNBTag); // Non-b-jet tagged non-b;
    if (BTags_.at(1)) out.push_back(RNBMTag); // Non-b-jet is tagged to be a b;
    else out.push_back(RNBTag); // Non-b-jet tagged non-b;
    if (BTags_.at(2)) out.push_back(RBTag); // b-jet tagged as a b;
    else out.push_back(RBMTag); // b-jet tagged to be a non-b;
    if (BTags_.at(3)) out.push_back(RBTag); // b-jet tagged as a b;
    else out.push_back(RBMTag); // b-jet tagged to be a non-b;
    if(BTags_.size() > 4){
      if (BTags_.at(4)) out.push_back(RBTag); // b-jet tagged as a b;
      else out.push_back(RBMTag); // b-jet tagged to be a non-b;
    }
    return out;
  }

  double CalcBTag(int it, vector<bool> BTags_, bool isB = true) {
    double pf = 1;
    if (isB) {
      if (BTags_.at(it)) pf *= RBTag;
      else pf *= RBMTag;
    }
    else {
      if (BTags_.at(it)) pf *= RNBMTag;
      else pf *= RNBTag;
    }
    return pf;
  }

  vector<TLorentzVector> GetPermutationLV(vector<int> perm_, vector<TLorentzVector> LVJets_) {
    vector<TLorentzVector> permlv;
    for (unsigned it = 0; it < perm_.size(); ++it) {
      permlv.push_back(LVJets_.at(perm_.at(it)));
    }
    // permlv.push_back(Lepton);
    // permlv.push_back(LVMET);
    return permlv;
  }
  //For minimizer

  double SolveNeutrinos(TLorentzVector LVLep, TLorentzVector ScaledMET, vector<TLorentzVector>& LVNeu_, bool debug_ = false) {
    bool debug = debug_;
    LVNeu_.clear();
    const double LepWMass = 80.4;
    double xn = ScaledMET.X();
    double yn = ScaledMET.Y();
    double xe = LVLep.X();
    double ye = LVLep.Y();
    double ze = LVLep.Z();
    double te = LVLep.T();
    double tar = LepWMass;
    double zn1, tn1, zn2, tn2;

    double tar2(tar*tar),tar4(tar2*tar2), xexn(xe*xn), yeyn(ye*yn), xe2(xe*xe), ye2(ye*ye), te2(te*te), xnye(xn*ye), xeyn(xe*yn);
    if (te2 - xe2 - ye2 - ze*ze != 0 && debug) cout << "Lepton mass is not zero!! : " << te2 - xe2 - ye2 - ze*ze <<endl;
    double v1 = 4.*pow((xnye - xeyn),2);
    double v2 = 4.*tar2*(xexn + yeyn);
    double radical = (tar4 - v1 + v2)*te2;

    if(debug) {
      cout << Form("Lepton X= %f, Y= %f, Z= %f, T= %f, MET X= %f Y= %f", xe, ye, ze, te, xn, yn)<<endl;
      cout << Form("tar2= %f, tar4= %f, xexn= %f, yeyn= %f, xe2= %f, ye2= %f, te2= %f, xnye= %f, xeyn= %f",tar2,tar4,xexn,yeyn,xe2,ye2,te2,xnye,xeyn)<<endl;
      cout << Form("radical = (tar4 - v1 + v2)*te2, here tar4 = %f, v1 = %f, v2 = %f, te2 = %f", tar4, v1, v2, te2)<<endl;
    }

    if (radical < 0) {
      LVNeu_.push_back(TLorentzVector());
      LVNeu_.push_back(TLorentzVector());
      return radical;
    }

    double coe = 1. / (2.*(xe2 + ye2));
    double a = (tar2 + 2.*xexn + 2.*yeyn) * ze;
    double b = sqrt(radical);
    zn1 = coe*(a + b);
    zn2 = coe*(a - b);
    tn1 = sqrt(zn1*zn1 + xn*xn + yn*yn);
    tn2 = sqrt(zn2*zn2 + xn*xn + yn*yn);
    TLorentzVector LVNeu1, LVNeu2;
    LVNeu1.SetXYZT(xn,yn,zn1,tn1);
    LVNeu2.SetXYZT(xn,yn,zn2,tn2);
    LVNeu_.push_back(LVNeu1);
    LVNeu_.push_back(LVNeu2);
    return 1;
  }

  vector<TLorentzVector> ScaleJets(vector<TLorentzVector> lvjets, const double *scales, TLorentzVector LVMET, TLorentzVector& ScaledMET) {
    vector<TLorentzVector> ScaledJets;
    ScaledMET = LVMET;
    ScaledJets.clear();
    for (unsigned ij = 0; ij < lvjets.size(); ++ij) {
      TLorentzVector newjet = lvjets.at(ij) * scales[ij];
      ScaledJets.push_back(newjet);
      ScaledMET = ScaledMET + lvjets.at(ij) - newjet;
    }
    return ScaledJets;
  }

  pair<double,double> CalcLimitsHist(double eta_, double pt_) {
    pair<int,int> bins = CalcBins(eta_, pt_);
    TH1F* hist = JESVector.at(bins.first).at(bins.second);
    double limits[2];
    double quantiles[2] = {0.023,0.977}; //range of 2 sigma : 95% events are within.
    hist->GetQuantiles(2,limits,quantiles);
    return pair<double, double>(limits[0],limits[1]);
  }

  pair<double,double> CalcLimitsFunc(double eta_, double pt_) {
    pair<int,int> bins = CalcBins(eta_, pt_);
    TF1* f1 = JESFuncVector.at(bins.first).at(bins.second);
    double mean = f1->GetParameter(1);
    double sigma = f1->GetParameter(2);
    return pair<double, double>(mean-2*sigma,mean+2*sigma);
  }

  double CalcPScaleHist(double eta_, double pt_, double scale_, int debug_ = 0) {
    pair<int,int> bins = CalcBins(eta_, pt_);
    TH1F* hist = JESVector.at(bins.first).at(bins.second);
    double n = hist->GetBinContent(hist->FindBin(scale_));
    double norm = hist->GetBinContent(hist->GetMaximumBin());
    double p = n/norm;
    if (debug_) cout <<endl << Form("Jet Pt = %f (%d), Eta = %f (%d), Scale = %f, n = %f, norm = %f, p = %f", pt_,bins.second, eta_, bins.first, scale_, n, norm, p) <<endl;
    return p;
  }

  double CalcPScaleFunc(double eta_, double pt_, double scale_, int debug_ = 0) {
    pair<int,int> bins = CalcBins(eta_, pt_);
    TF1* func = JESFuncVector.at(bins.first).at(bins.second);
    double n = func->Eval(scale_);
    double norm = func->GetMaximum();
    double p = n / norm;
    if (debug_) cout <<endl << Form("Jet Pt = %f (%d), Eta = %f (%d), Scale = %f, n = %f, norm = %f, p = %f", pt_,bins.second, eta_, bins.first, scale_, n, norm, p) <<endl;
    return p;
  }

  double CalcPScalesHist(vector<TLorentzVector> LVJets_, const double * scales, int debug_ = 0) {
    double PScale = 1;
    for (unsigned ij = 0; ij < LVJets_.size(); ++ij) {
      PScale *= CalcPScaleHist(LVJets_.at(ij).Eta(),LVJets_.at(ij).Pt(),scales[ij], debug_);
    }
    return PScale;
  }

  double CalcPScalesFunc(vector<TLorentzVector> LVJets_, const double * scales, int debug_ = 0) {
    double PScale = 1;
    for (unsigned ij = 0; ij < LVJets_.size(); ++ij) {
      PScale *= CalcPScaleFunc(LVJets_.at(ij).Eta(),LVJets_.at(ij).Pt(),scales[ij], debug_);
    }
    return PScale;
  }

  vector<double> CalcPScalesVectorFunc(vector<TLorentzVector> LVJets_, const double * scales, int debug_ = 0) {
    vector<double> out;
    for (unsigned ij = 0; ij < LVJets_.size(); ++ij) {
      double p = CalcPScaleFunc(LVJets_.at(ij).Eta(),LVJets_.at(ij).Pt(),scales[ij], debug_);
      out.push_back(p);
    }
    return out;
  }

  vector<double> CalcPMassVector(vector<TLorentzVector> jets) {
    double lepw = (jets[6] + jets[5]).M();
    double lept = (jets[6] + jets[5] + jets[3]).M();
    double hadw = (jets[0]+jets[1]).M();
    double hadt = (jets[0]+jets[1]+jets[2]).M();
    double plepw = CalcPWMass(lepw);
    double plept = CalcPTMass(lept);
    double phadw = CalcPWMass(hadw);
    double phadt = CalcPTMass(hadt);
    vector<double> out{phadw,phadt,plepw,plept};
    return out;
  }

  vector<double> CalcPMassVector(vector<TLorentzVector> jets_, TLorentzVector lep, TLorentzVector neu) {
    vector<TLorentzVector> jets = jets_;
    jets.push_back(lep);
    jets.push_back(neu);
    return CalcPMassVector(jets);
  }

  double CalcPW(TLorentzVector vv) {
    return WMassDis->Eval(vv.M());
  }

  double CalcPWMass(double mm) {
    return WMassDis->Eval(mm);
  }

  double CalcPT(TLorentzVector vv) {
    return TopMassDis->Eval(vv.M());
  }

  double CalcPTMass(double mm) {
    return TopMassDis->Eval(mm);
  }

  double CalcPLep(TLorentzVector LepB_, TLorentzVector Lep_, TLorentzVector Neu_) {
    return TopMassDis->Eval((Lep_ + Neu_ + LepB_).M()) * WMassDis->Eval((Lep_ + Neu_).M());
  }

  double CalcPLep(TLorentzVector LepB_, TLorentzVector Lep_, TLorentzVector LVMET, TLorentzVector &NeuOut) {
    vector<TLorentzVector> Neutrinos;
    double PNeutrino = SolveNeutrinos(Lep_, LVMET, Neutrinos);
    if (PNeutrino < 0) return PNeutrino;
    double PLep = 0;
    for (unsigned ineu = 0; ineu < Neutrinos.size(); ++ineu) {
      TLorentzVector Neu_ = Neutrinos.at(ineu);
      double PLepTMassTemp = TopMassDis->Eval((Lep_ + Neu_ + LepB_).M()) * WMassDis->Eval((Lep_ + Neu_).M());
      if (PLepTMassTemp > PLep) {
        PLep = PLepTMassTemp;
        NeuOut = Neu_;
      }
    }
    // PLep *= CalcPdPhi(LepB_.DeltaPhi(Lep_ + NeuOut));
    return PLep;
  }

  double CalcPLep(TLorentzVector LepB_, TLorentzVector Lep_, TLorentzVector LVMET, TLorentzVector GenNeu, TLorentzVector &NeuOut) {
    vector<TLorentzVector> Neutrinos;
    double PNeutrino = SolveNeutrinos(Lep_, LVMET, Neutrinos);
    if (PNeutrino < 0) return PNeutrino;
    double PLep = 0;
    TLorentzVector Neu_ = Neutrinos.at(0);
    if (Neutrinos.at(1).DeltaR(GenNeu) < Neu_.DeltaR(GenNeu)) Neu_ = Neutrinos.at(1);
    PLep = TopMassDis->Eval((Lep_ + Neu_ + LepB_).M()) * WMassDis->Eval((Lep_ + Neu_).M());
    NeuOut = Neu_;
    // PLep *= CalcPdPhi(LepB_.DeltaPhi(Lep_ + Neu_));
    return PLep;
  }

  double CalcPHad(vector<TLorentzVector> ScaledJets) {
    double PHadWMass = WMassDis->Eval((ScaledJets.at(0) + ScaledJets.at(1)).M());
    double PHadTMass = TopMassDis->Eval((ScaledJets.at(0) + ScaledJets.at(1) + ScaledJets.at(2)).M());
    double PHad = PHadWMass * PHadTMass;
    return PHad;
  }

  vector<int> Opt;
  void SetExtra(int opt_) {
    Opt.clear();
    while (opt_) {
      Opt.push_back(opt_%10);
      opt_ = opt_ / 10;
    }
  }

  // vector<double> CalcPExtra(vector<TLorentzVector> jets) {
  //   for (unsigned iopt = 0; iopt < Opt.size(); ++iopt) {
  //
  //   }
  // }

  // Snippets for evaluating dPhi of Leptonic W and b
  TFile* fdPhi;
  TH1D* hdPhi;
  void SetfdPhi() {
    fdPhi = new TFile("results/LepTopTreeRecoT_Add.root","READ");
    TH2F* hdPhi2D = (TH2F*) fdPhi->Get("WbdPhiVsMatching");
    hdPhi = hdPhi2D->ProjectionY("hdPhi",4,4);
  }
  double CalcPdPhi(double dphi) {
    dphi = fabs(dphi);
    int binx = hdPhi->GetXaxis()->FindBin(dphi);
    double p = hdPhi->GetBinContent(binx) / hdPhi->GetMaximum();
    return p;
  }

};

#endif
