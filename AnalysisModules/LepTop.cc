void LepTop_init(Analyzer *a) {
  a->CDOut();
  a->AddPlot(new TH1F("LepTopMatch","Leptonic Top Matching; Matched, Not Matched",2,0,2));
  a->AddPlot(new TH1F("LepTopMatchScaled","Leptonic Top Matching After Scaling;Matched, Not Matched",2,0,2));
  a->AddPlot(new TH1F("LepTopdR","Leptonic Top DeltaR with Gen",50,0.,5.));
  a->AddPlot(new TH1F("LepTopdRScaled","Leptonic Top DeltaR with Gen After Scaling",50,0.,5.));
}

int LepTop_loop(Analyzer *a) {
  if (a->AssignGenParticles() == -1) return -1;
  if (a->RecoPass == -1) return -1;
  map<string, TH1F*> p1d = a->Plots1D;

  a->CountEvent("PassPreSelection");

  vector<TLorentzVector> Jets = a->LVJets;
  vector<bool> BTags = a->BTags;
  if (Jets.size() != BTags.size()) {
    cout << "Jets size is different from BTags size" <<endl;
    return -1;
  }
  TLorentzVector Lepton = a->LVLeptons.at(0);
  TLorentzVector MET = a->LVMET;
  a->RM->SetLep(Lepton,MET);

  TLorentzVector GenLepT = a->LVGenLepT;
  double PLep = 0;
  TLorentzVector LepTopB = TLorentzVector();
  TLorentzVector LepTop = TLorentzVector();
  double PLepScaled = 0;
  TLorentzVector LepTopBScaled = TLorentzVector();
  TLorentzVector LepTopScaled = TLorentzVector();

  for (unsigned ij = 0; ij < Jets.size(); ++ij) {
    double PBtag;
    if (BTags.at(ij)) PBtag = 0.7;
    else PBtag = 0.3;
    TLorentzVector Jet = Jets.at(ij);
    TLorentzVector Neutrino = TLorentzVector();
    double JetPLep = a->JT->CalcPLep(Jet, Lepton, MET, Neutrino);
    if (JetPLep > PLep) {
      PLep = JetPLep;
      LepTopB = Jet;
      LepTop = Jet + Lepton + Neutrino;
    }
    Neutrino = TLorentzVector();
    double LepBScale = 0;
    double JetPLepScaled = a->RM->MinimizePLep(Jet,LepBScale);
    if (JetPLepScaled > PLepScaled) {
      PLepScaled = JetPLepScaled;
      LepTopBScaled = Jet * LepBScale;
      TLorentzVector ScaledMET = MET + Jet - LepTopBScaled;
      a->JT->CalcPLep(LepTopBScaled, Lepton, ScaledMET, Neutrino);
      LepTopScaled = LepTopBScaled + Lepton + Neutrino;
    }
  }

  if (PLep > 0) {
    double dr = LepTop.DeltaR(GenLepT);
    if (dr < 0.2) p1d["LepTopMatch"]->Fill(0);
    else p1d["LepTopMatch"]->Fill(1);
    p1d["LepTopdR"]->Fill(dr);
  }
  else p1d["LepTopMatch"]->Fill(1);

  if (PLepScaled > 0) {
    double dr = LepTopScaled.DeltaR(GenLepT);
    if (dr < 0.2) p1d["LepTopMatchScaled"]->Fill(0);
    else p1d["LepTopMatchScaled"]->Fill(1);
    p1d["LepTopdRScaled"]->Fill(dr);
  }
  else p1d["LepTopMatchScaled"]->Fill(1);

  return 0;
}
