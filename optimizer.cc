#include "optimizer.hh"

//Passed Variables
vector<TLorentzVector> Optimizer::LVJets;
vector<bool> Optimizer::BTags;
TLorentzVector Optimizer::LVLep;
TLorentzVector Optimizer::LVMET;
vector<double> Optimizer::etabins;
TFile* Optimizer::PFile;
vector<TH1F*> Optimizer::pJet;

//Constant Variables
double Optimizer::LepWMass;

//Intermediate Variables
vector<vector<int> > Optimizer::Permutations;
TLorentzVector Optimizer::LVNeu;
TLorentzVector Optimizer::LVNeu1;
TLorentzVector Optimizer::LVNeu2;
vector<TLorentzVector> Optimizer::OrderedJets;
vector<TLorentzVector> Optimizer::scaledjets;

//Optimize Result
double Optimizer::BestP;
double Optimizer::BestPJES;
double Optimizer::BestPMass;
double Optimizer::BestPType;
vector<int> Optimizer::BestPerm;
vector<double> Optimizer::BestScales;
TLorentzVector Optimizer::BestNeutrino;
ROOT::Math::Minimizer* Optimizer::mini = ROOT::Math::Factory::CreateMinimizer("TMinuit");
TF1* Optimizer::stdfgaus = new TF1("fgaus","gaus");
TF1* Optimizer::TopMassDis = new TF1("TBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,300.0);
TF1* Optimizer::WMassDis = new TF1("WBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,200.0);
bool Optimizer::resulting = false;
bool Optimizer::debug = false;




Optimizer::Optimizer(vector<TLorentzVector> LVJets_, vector<bool> BTags_){
  Initialize();
  LVJets = LVJets_;
  BTags = BTags_;
  GetPermutations();
  if (!mini) cout << "Fail to create minimizer" <<endl;
}

void Optimizer::Initialize(){
  TLorentzVector LVZero = TLorentzVector(0,0,0,0);
  // Passed Variables
  LVJets.clear();
  BTags.clear();
  LVLep = LVZero;
  LVMET = LVZero;
  // Constants
  LepWMass = 80.385;

  //Intermediate Variables
  Permutations.clear();
  LVNeu = LVZero;
  LVNeu1 = LVZero;
  LVNeu2 = LVZero;
  OrderedJets.clear();
  scaledjets.clear();

  //Optimize Results.
  BestP = 0;
  BestPerm.clear();
  BestScales.clear();
  BestNeutrino = LVZero;
  stdfgaus->SetParameters(1,0,1);
  TopMassDis->SetParameters(40.4,172.7,6.7);
  WMassDis->SetParameters(40.4,80.385,2.085);
}

void Optimizer::SetLepton(TLorentzVector lepton_) {
  LVLep = lepton_;
}

void Optimizer::SetMET(TLorentzVector MET_){
  LVMET = MET_;
}

void Optimizer::SetPFile(TFile* PFile_) {
  PFile = PFile_;
}

void Optimizer::SetEtaBins(vector<double> etabins_){
  etabins = etabins_;
}

void Optimizer::GetPermutations(){
  int jetsize = LVJets.size();
  for (int had1 = 0; had1 < jetsize - 1; ++had1) {
    for (int had2 = had1 + 1; had2 < jetsize; ++had2) {
      for (int hadtb = 0; hadtb < jetsize; ++hadtb) {
        if (hadtb == had1 || hadtb == had2) continue;
        for (int leptb = 0; leptb < jetsize; ++leptb) {
          if (leptb == had1 || leptb == had2 || leptb == hadtb) continue;
          for (int wpb = 0; wpb < jetsize; ++wpb) {
            if (wpb == leptb || wpb == hadtb || wpb == had1 || wpb == had2) continue;
            vector<int> perm;
            perm.push_back(had1);
            perm.push_back(had2);
            perm.push_back(hadtb);
            perm.push_back(leptb);
            perm.push_back(wpb);
            Permutations.push_back(perm);
          }
        }
      }
    }
  }
}

double Optimizer::GetpType(vector<int> perm_) {
  double pnbtag = 0.99;
  double pbtag = 0.7;
  double pnbmistag = 1.0 - pnbtag;
  double pbmistag = 1.0 - pbtag;
  double Ptype = (BTags[perm_[0]] ? pnbmistag : pnbtag) * (BTags[perm_[1]]? pnbmistag : pnbtag) * (BTags[perm_[2]] ? pbtag : pbmistag) * (BTags[perm_[3]] ? pbtag : pbmistag) * (BTags[perm_[4]] ? pbtag : pbmistag);
  return Ptype;
}

double Optimizer::SolveNeutrino(TLorentzVector ScaledMET) {
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

  double coe = 1. / (2.*(xe2 + ye2));
  double a = (tar2 + 2.*xexn + 2.*yeyn) * ze;
  double radical = (tar4 - 4.*pow((xnye - xeyn),2) + 4.*tar2*(xexn + yeyn))*te2;
  if (radical < 0) {
    // radical *= -1.0;
    return radical;
  }
  double b = sqrt(radical);
  zn1 = coe*(a + b);
  zn2 = coe*(a - b);
  tn1 = sqrt(zn1*zn1 + xn*xn + yn*yn);
  tn2 = sqrt(zn2*zn2 + xn*xn + yn*yn);
  LVNeu1.SetXYZT(xn,yn,zn1,tn1);
  LVNeu2.SetXYZT(xn,yn,zn2,tn2);
  return 1;
}

double Optimizer::GetpJES(const double *scales) {
  double PJES = 1.;
  double nfactor = stdfgaus->Eval(0);
  scaledjets.clear();
  TLorentzVector scaledMET = LVMET;
  for (unsigned ijet = 0; ijet < (OrderedJets.size() - 1); ++ijet) {
    PJES *= stdfgaus->Eval(scales[ijet]) / nfactor;
    double jeteta = fabs(OrderedJets[ijet].Eta());
    double jessigma(1), jesmean(0);
    for (unsigned ieta = 0; ieta < etabins.size() - 1; ++ieta) {
      if (jeteta < etabins[ieta + 1]) {
        int ibin = pJet[ieta]->FindBin(OrderedJets[ijet].Pt());
        jesmean = pJet[ieta]->GetBinContent(ibin);
        jessigma = pJet[ieta]->GetBinError(ibin);
        break;
      }
      if (ieta == etabins.size() - 1) {cout <<"Out of Eta Range"<<endl;return 0;}
    }
    scaledjets.push_back(OrderedJets[ijet] * (jesmean + scales[ijet] * jessigma));
    scaledMET += OrderedJets[ijet] - scaledjets[ijet];
  }

  double PMass, res;
  double PNeutrino = SolveNeutrino(scaledMET);
  if ( PNeutrino > 0 ) {
    PMass = GetpMass(scaledjets);
    res = -1.0 * PJES * PMass;
    if (resulting) {
      BestPJES = PJES;
      BestPMass = PMass;
    }
    if (debug){
      cout << Form("Scales: %f,%f,%f,%f",scales[0],scales[1],scales[2],scales[3])<<endl;
      cout <<"Positive, PMass= " << PMass << "; PJES = "<<PJES <<endl;
    }

  }
  else {
    if (resulting) {
      BestPJES = PJES;
      BestPMass = -1;
    }
    PMass = 0;
    res = (-1.0 * PNeutrino);
  }
  // PMass = 1.;
  return res;
}

double Optimizer::GetpMass(vector<TLorentzVector> lvjets) {
  double HadWmass = (lvjets[0]+lvjets[1]).M();
  TLorentzVector LVHadT = lvjets[0]+lvjets[1]+lvjets[2];
  double HadTmass = LVHadT.M();
  TLorentzVector LVLepW1 = LVNeu1 + LVLep;
  TLorentzVector LVLepW2 = LVNeu2 + LVLep;
  double LepTmass1 = (LVLepW1 + lvjets[3]).M();
  double LepTmass2 = (LVLepW2 + lvjets[3]).M();

  bool UseHist = false;
  double PHadWmass, PHadTmass, PLepTmass1, PLepTmass2;
  //Use Filled Histograms.
  if (UseHist){
    TH1F* pGenTMass = (TH1F*) PFile->Get("pGenTMass");
    TH1F* pGenWMass = (TH1F*) PFile->Get("pGenWMass");
    int HadWmassBin = pGenWMass->FindBin(HadWmass);
    int HadTmassBin = pGenTMass->FindBin(HadTmass);
    PHadWmass = pGenWMass->GetBinContent(HadWmassBin);
    PHadTmass = pGenTMass->GetBinContent(HadTmassBin);

    // cout << Form("W mass = %5.1f, Bin: %d, Value: %5.1f; T Mass = %5.1f, Bin: %d, Value: %5.1f;  P = %f", HadWmass, HadWmassBin, PHadWmass, HadTmass, HadTmassBin, PHadTmass, PHadMass) << endl;
    int LepTmassBin1 = pGenTMass->FindBin(LepTmass1);
    int LepTmassBin2 = pGenTMass->FindBin(LepTmass2);
    PLepTmass1 = pGenTMass->GetBinContent(LepTmassBin1);
    PLepTmass2 = pGenTMass->GetBinContent(LepTmassBin2);
  }

  //Use Standard Mass distribution functions.
  else {
    PHadWmass = WMassDis->Eval(HadWmass)/WMassDis->GetMaximum();
    PHadTmass = TopMassDis->Eval(HadTmass)/TopMassDis->GetMaximum();
    PLepTmass1 = TopMassDis->Eval(LepTmass1)/TopMassDis->GetMaximum();
    PLepTmass2 = TopMassDis->Eval(LepTmass2)/TopMassDis->GetMaximum();
  }
  TH1F* pGenWPdPhi = (TH1F*) PFile->Get("pGenWPdPhi");
  double dPhi = LVHadT.DeltaPhi(OrderedJets[4]);
  int dPhiBin = pGenWPdPhi->FindBin(dPhi);
  double PPhi = pGenWPdPhi->GetBinContent(dPhiBin);

  double PHadMass = PHadWmass * PHadTmass * PPhi;
  double PLepMass;
  if (PLepTmass1 > PLepTmass2) {
    PLepMass = PLepTmass1;
    LVNeu = LVNeu1;
  }
  else {
    PLepMass = PLepTmass2;
    LVNeu = LVNeu2;
  }
  double PMass = PHadMass * PLepMass;

  // cout << Form("T1 mass = %5.1f, Bin: %d, Value: %5.1f; T2 Mass = %5.1f, Bin: %d, Value: %5.1f;  P = %f; TMET Mass: %5.1f", LepTmass1, LepTmassBin1, PLepTmass1, LepTmass2, LepTmassBin2, PLepTmass2, PLepMass, LVLepTMET.M()) << endl;
  return PMass;
}

double Optimizer::OptimizeThisPerm(vector<double> &ThisScale_) {
  for (unsigned iscale = 0; iscale < (OrderedJets.size() - 1); ++iscale) {
    mini->SetLimitedVariable(iscale,Form("Scale_%i",iscale),0.0,0.05,-2.0,2.0);
  }
  mini->Minimize();
  ThisScale_.clear();
  if (!(mini->Status())) {
    for (unsigned iscale = 0; iscale < (OrderedJets.size() - 1); ++iscale) {
      ThisScale_.push_back(mini->X()[iscale]);
    }
    if (debug) {
      cout << "Minimizer Worked" <<endl;
    }

    return (-1.0 * (mini->MinValue()) );
  }
  if (debug) {
    cout <<"No Minimizer" <<endl;
  }
  return 0;
}

void Optimizer::GradientOptimize() {

}

void Optimizer::FindBestPerm() {
  ROOT::Math::Functor f(&GetpJES,4);
  mini->SetPrintLevel(0);
  mini->SetStrategy(2);
  mini->SetMaxFunctionCalls(400);
  mini->SetMaxIterations(400);
  mini->SetTolerance(0.01);
  mini->SetErrorDef(0.5);
  mini->SetFunction(f);
  if (debug) {
    cout << "BestP:" << BestP <<endl;
  }

  for (unsigned iperm = 0; iperm < Permutations.size(); ++iperm) {
    vector<int> ThisPerm = Permutations.at(iperm);
    if (debug) {
      cout <<Form("ThisPerm: %i,%i,%i,%i, %i \n" , ThisPerm[0], ThisPerm[1], ThisPerm[2], ThisPerm[3],ThisPerm[4]);
    }

    double pType = GetpType(ThisPerm);
    OrderedJets.clear();
    for (unsigned ijet = 0; ijet < ThisPerm.size(); ++ijet) {
      OrderedJets.push_back(LVJets[ThisPerm[ijet]]);
    }
    vector<double> ThisScale;
    double ThisP = OptimizeThisPerm(ThisScale) * pType;
    if (ThisP > BestP) {
      BestPType = pType;
      BestP = ThisP;
      BestPerm = ThisPerm;
      BestScales = ThisScale;
    }
  }
  // cout << Form("P: %f, BestPerm: %i,%i,%i,%i; scales: %f, %f, %f ,%f",BestP, BestPerm[0], BestPerm[1], BestPerm[2], BestPerm[3], BestScales[0],BestScales[1],BestScales[2],BestScales[3])<<endl;
}

void Optimizer::Optimize() {
  pJet.clear();
  for (unsigned ieta = 0; ieta < etabins.size() - 1; ++ieta) {
    TString pjetname = Form("pJet_%d",ieta);
    pJet.push_back((TH1F*)PFile->Get(pjetname));
  }
  FindBestPerm();
}

double Optimizer::GetBestP() {
  // cout << "BESTP = " << BestP <<endl;
  return BestP;
}

vector<int> Optimizer::GetBestPerm() {
  return BestPerm;
}

vector<double> Optimizer::GetBestScales() {
  return BestScales;
}

vector<TLorentzVector> Optimizer::GetBestLVJets() {
  OrderedJets.clear();
  for (unsigned ijet = 0; ijet < BestPerm.size(); ++ijet) {
    OrderedJets.push_back(LVJets[BestPerm[ijet]]);
  }
  double scales[BestScales.size()];
  for (unsigned i = 0; i < BestScales.size(); ++ i) {
    scales[i] = BestScales[i];
  }
  // cout << Form("scales: %f, %f, %f ,%f ,%f",scales[0],scales[1],scales[2],scales[3],scales[4])<<endl;
  //Getting the intermediate quantities.
  resulting = true;
  GetpJES(scales);
  resulting = false;

  BestNeutrino = LVNeu;
  vector<TLorentzVector> out = scaledjets;
  out.push_back(OrderedJets[4]);
  return out;
}

double Optimizer::GetBestPJES() {
  return BestPJES;
}

double Optimizer::GetBestPMass() {
  return BestPMass;
}

double Optimizer::GetBestPType() {
  return BestPType;
}

TLorentzVector Optimizer::GetBestNeutrino() {
  return BestNeutrino;
}

TLorentzVector Optimizer::GetWPrime() {
  TLorentzVector lvwprime = scaledjets[0] + scaledjets[1] + scaledjets[2] + OrderedJets[4];
  return lvwprime;
}
