void LepTop_init(Analyzer *a) {
  a->CDOut();
  a->AddPlot(new TH1F("LepTopMatch","Leptonic Top Matching; Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepTopMatchScaled","Leptonic Top Matching After Scaling;Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepTopdR","Leptonic Top DeltaR with Gen",50,0.,5.));
  a->AddPlot(new TH1F("LepTopdRScaled","Leptonic Top DeltaR with Gen After Scaling",50,0.,5.));
  a->AddPlot(new TH1F("LepbMatch","Leptonic b Matching; Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepbMatchScaled","Leptonic b Matching After Scaling;Matched, Not Matched",3,0,3));
  a->AddPlot(new TH1F("LepbdR","Leptonic b DeltaR with Gen",50,0.,5.));
  a->AddPlot(new TH1F("LepbdRScaled","Leptonic b DeltaR with Gen After Scaling",50,0.,5.));

  a->AddPlot(new TH2F("NeutrinodRVsLepbdR","Neutrino dR Vs Lepb dR;LepbdR;Neutrino dR",50,0.,5.,50,0.,5.));
  a->AddPlot(new TH2F("NeutrinodRVsLepbdRScaled","Neutrino dR Vs Lepb dR After Scale;LepbdR;Neutrino dR",50,0.,5.,50,0.,5.));
  a->AddPlot(new TH2F("LepbPtRatioVsdR","LepbPtRatio Vs Lepb dR;LepbdR;Pt Ratio",50,0.,5.,80,0.,4.));
  a->AddPlot(new TH2F("LepbPtRatioVsdR_PtLimited","LepbPtRatio Vs Pt Limited Lepb dR;LepbdR;Pt Ratio",50,0.,5.,40,0.,2.));
  a->AddPlot(new TH2F("WbdPhiVsMatching","dPhi of Lep W and b Vs matching;Matched, Not Matched, Not Found, Matched Scaled, Not Matched Scaled, Not Found Scaled",6,0,6,40,0.,4.));
  a->AddPlot(new TH2F("PtopVsMatching","Ptop Vs Matching;Matched, Not Matched, Not Found, Matched Scaled, Not Matched Scaled, Not Found Scaled",6,0,6,201,-20.,0.1));
}

int LepTop_loop(Analyzer *a) {
  if (a->AssignGenParticles() == -1) return -1;
  if (a->RecoPass == -1) return -1;
  map<string, TH1F*> p1d = a->Plots1D;
  map<string, TH2F*> p2d = a->Plots2D;

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
  TLorentzVector GenLepB = a->LVGenLepB;
  TLorentzVector GenNeu = a->LVGenNeu;
  double PLep = 0;
  TLorentzVector LepTopB = TLorentzVector();
  TLorentzVector LepTop = TLorentzVector();
  double PLepScaled = 0;
  TLorentzVector LepTopBScaled = TLorentzVector();
  TLorentzVector LepTopScaled = TLorentzVector();
  TLorentzVector Neutrino = TLorentzVector();
  TLorentzVector NeutrinoScaled = TLorentzVector();
  TLorentzVector MatchGenb = TLorentzVector();
  TLorentzVector MatchGenbLimited = TLorentzVector();
  double MatchGenbdR = 6.;
  double MatchGenbdRLimited = 6.;

  for (unsigned ij = 0; ij < Jets.size(); ++ij) {
    double PBtag;
    if (BTags.at(ij)) PBtag = 0.7;
    else PBtag = 0.3;
    TLorentzVector Jet = Jets.at(ij);
    TLorentzVector Neutrinotemp = TLorentzVector();
    double JetPLep = a->JT->CalcPLep(Jet, Lepton, MET, Neutrinotemp);
    if (JetPLep > PLep) {
      PLep = JetPLep;
      LepTopB = Jet;
      LepTop = Jet + Lepton + Neutrinotemp;
      Neutrino = Neutrinotemp;
    }
    Neutrinotemp = TLorentzVector();
    double LepBScale = 0;
    double JetPLepScaled = a->RM->MinimizePLep(Jet,LepBScale);
    if (JetPLepScaled > PLepScaled) {
      PLepScaled = JetPLepScaled;
      LepTopBScaled = Jet * LepBScale;
      TLorentzVector ScaledMET = MET + Jet - LepTopBScaled;
      a->JT->CalcPLep(LepTopBScaled, Lepton, ScaledMET, Neutrinotemp);
      LepTopScaled = LepTopBScaled + Lepton + Neutrinotemp;
      NeutrinoScaled = Neutrinotemp;
    }

    double GendR = Jet.DeltaR(GenLepB);
    double GenPtRatio =  Jet.Pt() / GenLepB.Pt();
    bool PtInRange = (GenPtRatio <= 1.5 && GenPtRatio >= 0.5);
    if (GendR < MatchGenbdR) {
      MatchGenbdR = GendR;
      MatchGenb = Jet;
    }
    if (GendR < MatchGenbdRLimited && PtInRange) {
      MatchGenbdRLimited = GendR;
      MatchGenbLimited = Jet;
    }
  }

  p2d["LepbPtRatioVsdR"]->Fill(MatchGenbdR, MatchGenb.DeltaR(GenLepB));
  p2d["LepbPtRatioVsdR_PtLimited"]->Fill(MatchGenbdRLimited, MatchGenbLimited.DeltaR(GenLepB));

  if (PLep > 0) {
    double tdr = LepTop.DeltaR(GenLepT);
    // if (tdr < 0.2) p1d["LepTopMatch"]->Fill(0);
    // else p1d["LepTopMatch"]->Fill(1);
    p1d["LepTopMatch"]->Fill(tdr < 0.2 ? 0 : 1);
    double bdr = LepTopB.DeltaR(GenLepB);
    int bmatch = (bdr < 0.2 ? 0 : 1);
    p1d["LepbMatch"]->Fill(bmatch);
    p1d["LepTopdR"]->Fill(tdr);
    p1d["LepbdR"]->Fill(bdr);
    double neudr = Neutrino.DeltaR(GenNeu);
    p2d["NeutrinodRVsLepbdR"]->Fill(neudr, bdr);
    double wbdphi = fabs((Lepton + Neutrino).DeltaPhi(LepTopB));
    p2d["WbdPhiVsMatching"]->Fill(bmatch,wbdphi);
    p2d["PtopVsMatching"]->Fill(bmatch, log10(PLepScaled));

  }
  else p1d["LepTopMatch"]->Fill(2);

  if (PLepScaled > 0) {
    double tdr = LepTopScaled.DeltaR(GenLepT);
    if (tdr < 0.2) p1d["LepTopMatchScaled"]->Fill(0);
    else p1d["LepTopMatchScaled"]->Fill(1);
    double bdr = LepTopBScaled.DeltaR(GenLepB);
    int bmatch = (bdr < 0.2 ? 0 : 1);
    p1d["LepbMatchScaled"]->Fill(bmatch);
    p1d["LepTopdRScaled"]->Fill(tdr);
    p1d["LepbdRScaled"]->Fill(bdr);
    double neudr = NeutrinoScaled.DeltaR(GenNeu);
    p2d["NeutrinodRVsLepbdRScaled"]->Fill(neudr, bdr);
    double wbdphi = fabs((Lepton + NeutrinoScaled).DeltaPhi(LepTopBScaled));
    p2d["WbdPhiVsMatching"]->Fill(bmatch + 3, wbdphi);
    p2d["PtopVsMatching"]->Fill(bmatch+3,log10(PLepScaled));
  }
  else p1d["LepTopMatchScaled"]->Fill(2);

  return 0;
}
