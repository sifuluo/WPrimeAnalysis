#ifndef ETAPTBINS_CC
#define ETAPTBINS_CC

#include<vector>
#include<cmath>

using namespace std;

class EtaPtBins{
public:
  EtaPtBins(){};

  const vector<double> etabins{0., 1.3, 2.5, 3.0, 5.2};

  const vector<vector<double> > ptbins{
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130., 150.,180.,220., 260., 300.,350.,400.,500.,1000.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150., 180.,220.,260., 300.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150., 180.,220.,260.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.}
  };

  vector< vector<TH1F*> > JESVector;

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

  int iEta(double eta_) {
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

  pair<int,int> iBin(double eta_, double pt_) {
    int iEta(0), iPt(0);

    for (unsigned ieta = 0; ieta < (etabins.size() -1); ++ ieta ) {
      if (eta_ < etabins.at(ieta+1) ) {
        iEta = ieta;
        break;
      }
      if (ieta == (etabins.size() -2)) iEta = ieta;
    }

    for (unsigned ipt = 0; ipt < (ptbins.at(iEta).size() -1); ++ ipt ) {
      if (pt_ < ptbins.at(iEta).at(ipt+1) ) {
        iPt = ipt;
        break;
      }
      if (ipt == (ptbins.at(iEta).size() -2)) iPt = ipt;
    }

    return pair<int,int>(iEta, iPt);
  }

  vector< vector<TH1F*> > MakeJES(){
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

  vector< vector<TH1F*> > ReadJES(TFile* f) {
    vector< vector<TH1F*> > jes;
    jes.clear();
    for (unsigned ieta = 0; ieta < EtaBins().size()-1; ++ieta) {
      vector<TH1F*> jeseta;
      jeseta.clear();
      for (unsigned ipt = 0; ipt < PtBins(ieta).size() -1; ++ipt){
        TString sn = Form("eta%d_pt%d", ieta, ipt);
        jeseta.push_back( (TH1F*)(f->Get(sn))); //Histogram might not be accessible after TFile being closed
        // TH1F hjes = *((TH1F*)(f->Get(sn)));
        // jeseta.push_back(&hjes);
      }
      jes.push_back(jeseta);
    }
    JESVector = jes;
    return jes;
  }

  TH1F* GetPlot(double eta_, double pt_) {
    return JESVector.at(iBin(eta_, pt_).first).at(iBin(eta_,pt_).second);
  }

};



#endif
