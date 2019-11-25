#ifndef JESTOOLS_CC
#define JESTOOLS_CC

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>

#include <utility>
#include <vector>
#include <cmath>

using namespace std;

class JESTools{
public:
  JESTools(){
    // JESVector.clear();
    // TempiEta = 0;
    // TempiPt = 0;
  };

  const vector<double> etabins{0., 1.3, 2.5, 3.0, 5.2};

  const vector<vector<double> > ptbins{
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130., 150.,180.,220., 260., 300.,350.,400.,500.,1000.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150., 180.,220.,260., 300.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150., 180.,220.,260.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.}
  };

  TF1* TopMassDis = new TF1("TBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,300.0);
  TF1* WMassDis = new TF1("WBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,200.0);

  vector< vector<TH1F*> > JESVector;


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
      if (eta_ < etabins.at(ieta+1) ) {
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
      if (eta_ < etabins.at(ieta+1) ) {
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

  vector< vector<TH1F*> > MakeJESPlots(){
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

  TH1F* GetPlot() {
    TempHist = JESVector.at(TempiEta).at(TempiPt);
    return TempHist;
  }

  void FillPlot(double fill, int ieta_ = TempiEta, int ipt_ = TempiPt) {
    JESVector.at(ieta_).at(ipt)->Fill(fill);
  }

  //Below is for Minimizers
  vector< vector<TH1F*> > ReadJESPlots(TFile* f) {
    SetUpMassFunctions();
    vector< vector<TH1F*> > jes;
    jes.clear();
    for (unsigned ieta = 0; ieta < EtaBins().size()-1; ++ieta) {
      vector<TH1F*> jeseta;
      jeseta.clear();
      for (unsigned ipt = 0; ipt < PtBins(ieta).size() -1; ++ipt){
        TString sn = Form("eta%d_pt%d", ieta, ipt);
        jeseta.push_back( (TH1F*)(f->Get(sn))); //Histogram might not be accessible after TFile being closed
      }
      jes.push_back(jeseta);
    }
    JESVector = jes;
    return jes;
  }

  void SetUpMassFunctions() {
    TopMassDis = new TF1("TBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,300.0);
    WMassDis = new TF1("WBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,200.0);
    TopMassDis->SetParameters(100.,172.7,6.7);
    WMassDis->SetParameters(100,80.385,2.085);
    TopMassDis->SetParameter(0,100./TopMassDis->Eval(172.7)); // normalized it to peak at y = 1;
    WMassDis->SetParameter(0,100./WMassDis->Eval(80.385)); // normalized it to peak at y = 1;
  }

  double SolveNeutrinos(TLorentzVector LVLep, TLorentzVector ScaledMET, vector<TLorentzVector>& LVNeu_) {
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
    // if (te2 - xe2 - ye2 - ze2 != 0 && verbose) cout << "Lepton mass is not zero!! : " << te2 - xe2 - ye2 - ze2 <<endl;
    double radical = (tar4 - 4.*pow((xnye - xeyn),2) + 4.*tar2*(xexn + yeyn))*te2;
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

  vector<TLorentzVector> ScaleJets(vector<TLorentzVector> lvjets, double *scales, TLorentzVector LVMET, TLorentzVector& ScaledMET) {
    vector<TLorentzVector> scaledjets;
    ScaledMET = LVMET;
    scaledjets.clear();
    for (unsigned ij = 0; ij < lvjets.size(); ++ij) {
      TLorentzVector newjet = lvjets.at(ij) * scales[ij];
      scaledjets.push_back(newjet);
      ScaledMET = ScaledMET + lvjets.at(ij) - newjet;
    }
    return scaledjets;
  }

  vector< vector<int> > MakePermutations5(int jetsize) {
    vector< vector<int> > Permutations;
    Permutations.clear();
    for (unsigned ihad1 = 0; ihad1 < jetsize - 1; ++ihad1) {
      for (unsigned ihad2 = ihad1 + 1; ihad2 < jetsize; ++ihad2) {
        for (unsigned ihadb = 0; ihadb < jetsize; ++ihadb) {
          if (ihadb == ihad1 || ihadb == ihad2) continue;
          for (unsigned ilepb = 0; ilepb < jetsize; ++ ilepb) {
            if (ilepb == ihad1 || ilepb == ihad2 || ilepb == ihadb) continue;
            for (unsigned iwpb = 0; iwpb < jetsize; ++ iwpb) {
              if (iwpb == ihad1 || iwpb == ihad2 || iwpb == ihadb || iwpb == ilepb) continue;
              vector<int> perm{ihad1, ihad2, ihadb, ilepb, iwpb};
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
    for (unsigned ihad1 = 0; ihad1 < jetsize - 1; ++ihad1) {
      for (unsigned ihad2 = ihad1 + 1; ihad2 < jetsize; ++ihad2) {
        for (unsigned ihadb = 0; ihadb < jetsize; ++ihadb) {
          if (ihadb == ihad1 || ihadb == ihad2) continue;
          for (unsigned ilepb = 0; ilepb < jetsize; ++ ilepb) {
            if (ilepb == ihad1 || ilepb == ihad2 || ilepb == ihadb) continue;
            vector<int> perm{ihad1, ihad2, ihadb, ilepb};
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

  double CalcPFlavor(vector<int> perm, vector<bool> BTags) {
    const double RNBMTag(0.01), RNBTag(0.99), RBMTag(0.3), RBTag(0.7);
    double pf = 1;
    if (Btags.at(perm_.at(0))) pf *= RNBMTag; // Non-b-jet is tagged to be a b;
    else pf *= RNBTag; // Non-b-jet tagged non-b;
    if (Btags.at(perm_.at(1))) pf *= RNBMTag; // Non-b-jet is tagged to be a b;
    else pf *= RNBTag; // Non-b-jet tagged non-b;
    if (BTags.at(perm_.at(2))) pf *= RBTag; // b-jet tagged as a b;
    else pf *= RBMTag; // b-jet tagged to be a non-b;
    if (BTags.at(perm_.at(3))) pf *= RBTag; // b-jet tagged as a b;
    else pf *= RBMTag; // b-jet tagged to be a non-b;
    if(perm_.size() > 4){
      if (BTags.at(perm_.at(4))) pf *= RBTag; // b-jet tagged as a b;
      else pf *= RBMTag; // b-jet tagged to be a non-b;
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

  pair<double,double> CalcLimits(double eta_, double pt_) {
    pair<int,int> bins = CalcBins(eta_, pt_);
    TH1F* hist = JESVector.at(bins.first).at(bins.second);
    double limits[2];
    double quantiles[2] = {0.023,0.977}; //range of 2 sigma : 95% events are within.
    hist->GetQuantiles(2,limits,quantiles);
    return pair<double, double>(limits[0],limits[1]);
  }

  double CalcPScale(double eta_, double pt_, double scale_) {
    pair<int,int> bins = CalcBins(eta_, pt_);
    TH1F* hist = JESVector.at(bins.first).at(bins.second);
    double n = hist->GetBinContent(hist->FindBin(scale_));
    double norm = hist->GetBinContent(hist->GetMaximumBin());
    double p = n/norm;
    return p;
  }

  double CalcPScales(vector<TLorentzVector> LVJets_, double * scales) {
    double PScale = 1;
    for (unsigned ij = 0; ij < LVJets_.size(); ++ij) {
      PScale *= CalcPScale(LVJets_.at(ij).Eta(),LVJets_.at(ij).Pt(),scales[ij]);
    }
    return PScale;
  }

  double CalcPLep(TLorentzVector LepB_, TLorentzVector Lep_, TLorentzVector Neu_) {
    return TopMassDis->Eval((Lep_ + Neu_ + LepB_).M());
  }

  double CalcPHad(vector<TLorentzVector> ScaledJets) {
    double PHadWMass = WMassDis->Eval((ScaledJets.at(0) + ScaledJets.at(1)).M());
    double PHadTMass = TopMassDis->Eval((ScaledJets.at(0) + ScaledJets.at(1) + ScaledJets.at(2)).M());
    double PHad = PHadWMass * PHadTMass;
    return PHad;
  }


};

#endif
