void TemplateTest2d_init(Analyzer *a) {
  a->CDOut();
  a->AddPlot(new TH2F("WPrimeMasses2D","W\' Masses; Mass if FL; Mass if LL",1000,0,1000,1000,0,1000));
  a->AddPlot(new TH1F("PTTbar","Best Probability of ttbar hypothesis; Log(P)", 101,-10,0.1));
  a->AddPlot(new TH2F("Probabilities","Best Probabilities; Log(Prob) if FL; Log(Prob) if LL",101,-10,0.1,100,-10,0.1));
  a->AddPlot(new TH1F("WPrimeMassFL","W\' mass Reconstructed As FL",1000,0,1000));
  a->AddPlot(new TH1F("WPrimeMassLL","W\' mass Reconstructed As LL",1000,0,1000));
  a->AddPlot(new TH1F("WPrimeMassFlagFL","W\' mass Flagged As FL",1000,0,1000));
  a->AddPlot(new TH1F("WPrimeMassFlagLL","W\' mass Flagged As LL",1000,0,1000));
}

int TemplateTest2d_loop(Analyzer *a, int& DoReco, int LeadingAsWPB = 0) {
  if (a->AssignGenParticles() == -1) return -1;
  if (a->RecoPass == -1) return -1; //lepton != 1, jet < 5 : discard
  map<string, TH1F*> p1d = a->Plots1D;
  map<string, TH2F*> p2d = a->Plots2D;

  vector<TLorentzVector> Jets; //Inputs
  vector<bool> BTags;
  if (DoReco == 0) {
    Jets = a->LVOutPart;
    BTags = a->GenOutBTags;
    a->RM->SetLep(a->LVGenLep, a->LVGenNeu);
  }
  if (DoReco == 1) {
    Jets = a->LVJets;
    BTags = a->BTags;
    a->RM->SetLep(a->LVLeptons.at(0), a->LVMET);
  }

  if (BTags.size() != Jets.size()) {
    cout <<endl;
    cout << "Jets size = " << Jets.size() <<endl;
    cout << "BTags size = " << BTags.size() <<endl;
  }
  if (Jets.size() < 5) return -1;
  a->CountEvent(1); // Event Passed selection

  vector<int> CorrectPerm; // Correct Permutation
  if (DoReco == 0) {
    a->GetGenCorrectPerm();
    CorrectPerm = a->GenCorrectPerm;
  }
  if (DoReco == 1) {
    a->GetRecoCorrectPerm();
    CorrectPerm = a->RecoCorrectPerm;
  }
  TLorentzVector LeadingJet = Jets.at(0);
  int LeadingIndex = 0;
  for (unsigned ij = 1; ij < Jets.size(); ++ij) {
    if (Jets.at(ij).Pt() > LeadingJet.Pt()) {
      LeadingIndex = ij;
      LeadingJet = Jets.at(ij);
    }
  }
  if (LeadingAsWPB) {
    Jets.at(LeadingIndex) = TLorentzVector();
  }

  vector<int> BestPerm;
  pair<double, vector<TLorentzVector> > ttbar = a->SolveTTbar(Jets, BTags, BestPerm);

  double BestP = ttbar.first;
  vector<TLorentzVector> BestParticles = ttbar.second;

  if (BestP == -1.) return -1;

  a->CountEvent("Has TTbar"); //Event has a viable ttbar hypothesis
  p1d["PTTbar"]->Fill(log10(BestP));

  TLorentzVector hadt = BestParticles[0] + BestParticles[1] + BestParticles[2];
  TLorentzVector lept = BestParticles[3] + BestParticles[4] + BestParticles[5];

  // Deterimine the sample by together Top dPhi, Top dPt and W' top Pt.
  double pTopFL(1.), pTopLL(1.);
  double topdphi = fabs(hadt.DeltaPhi(lept));
  double topdpt = hadt.Pt() - lept.Pt();
  pTopFL *= a->GT1->CalcP("TopdPhi",topdphi,0) * a->GT1->CalcP("TopPtDiff",topdpt,0) * a->GT1->CalcP("TopHadPt",hadt.Pt(),0);
  pTopLL *= a->GT1->CalcP("TopdPhi",topdphi,1) * a->GT1->CalcP("TopPtDiff",topdpt,1) * a->GT1->CalcP("TopLepPt",lept.Pt(),1);

  double BestPFL(-1.), BestPLL(-1.);
  TLorentzVector BestFLWPB, BestLLWPB;
  vector<int> WPBCand = a->JT->FindWPB(Jets.size(),BestPerm);
  for (unsigned iwpb = 0; iwpb < WPBCand.size(); ++iwpb) {
    double pWPB = 1.;
    TLorentzVector LVWPb = Jets.at(WPBCand.at(iwpb));
    if (LeadingAsWPB) LVWPb = LeadingJet;
    double pWPBTag = a->JT->CalcBTag(WPBCand.at(iwpb), BTags, true);
    double pWPBFL = a->GT1->CalcP("WPBPt",LVWPb.Pt(),0) * a->GT2->CalcP("WPdPhi", fabs(LVWPb.DeltaPhi(hadt)),0) * pTopFL;
    double pWPBLL = a->GT1->CalcP("WPBPt",LVWPb.Pt(),1) * a->GT2->CalcP("WPdPhi", fabs(LVWPb.DeltaPhi(hadt)),1) * pTopLL;

    if (pWPBFL > BestPFL) {
      BestPFL = pWPBFL;
      BestFLWPB = LVWPb;
    }
    if (pWPBLL > BestPLL) {
      BestPLL = pWPBLL;
      BestLLWPB = LVWPb;
    }
    if (LeadingAsWPB) break;
  }
  TLorentzVector FLWp = BestFLWPB + hadt;
  TLorentzVector LLWp = BestLLWPB + lept;
  double FLWpMass = FLWp.M();
  double LLWpMass = LLWp.M();
  p2d["WPrimeMasses2D"]->Fill(FLWpMass, LLWpMass);
  p2d["Probabilities"]->Fill(log10(BestP*BestPFL),log10(BestP*BestPLL));
  p1d["WPrimeMassFL"]->Fill(FLWpMass);
  p1d["WPrimeMassLL"]->Fill(LLWpMass);
  int SampleFlag = 0;
  if (BestPLL > BestPFL) SampleFlag = 1;
  if (!SampleFlag) {
    p1d["WPrimeMassFlagFL"]->Fill(FLWpMass);
  }
  else {
    p1d["WPrimeMassFlagLL"]->Fill(LLWpMass);
  }
  return 0;
}
