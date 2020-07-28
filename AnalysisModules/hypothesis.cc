TH1F *EventCounter2 = new TH1F("EventCounter2","Number2 of Events;Passcode",10,-0.5,9.5);

void hypothesis_init(Analyzer *a) {
  // Particle Numbers
  a->AddPlot(new TH1F("EventCounter","Number of Events;Passcode",10,-0.5,9.5));
  a->AddPlot(new TH1F("NumberW","Number of W",10,-0.5,9.5));
  a->AddPlot(new TH1F("Numbert","Number of t",10,-0.5,9.5));
  a->AddPlot(new TH1F("Numberb","Number of b",10,-0.5,9.5));
  a->AddPlot(new TH1F("NumberWP","Number of Wp",10,-0.5,9.5));

  // Tops
  a->AddPlot(new TH1F("TopdR","TopdR;dR",80,0,4));
  a->AddPlot(new TH1F("TopPtDiff","W' t minus Other t, hadronic t minus leptonic t (bg)",400,-1000,1000));
  a->AddPlot(new TH1F("WPTPt","W' t Pt",200,0,1000));
  a->AddPlot(new TH1F("OtTPt","Other t Pt",200,0,1000));

  // Ws
  a->AddPlot(new TH1F("WdR","W dR;dR",200,0,10));
  a->AddPlot(new TH1F("WPtDiff","W Pt difference",200,0,1000));

  // Observables
  a->AddPlot(new TH1F("WPBPt","W' b Pt",200,0,1000));
  a->AddPlot(new TH1F("WPTBPt","W' t b Pt",200,0,1000));
  a->AddPlot(new TH1F("OtTBPt","Other t b Pt",200,0,1000));
  a->AddPlot(new TH1F("LFPt","Light flavor Pt",200,0,1000));
  a->AddPlot(new TH1F("WPBLeading","W' b leading?",2,-0.5,1.5));
  a->AddPlot(new TH1F("JetNumber","Number of Jets",20,-0.5,19.5));
  a->AddPlot(new TH1F("BJetNumber","Number of b-Jets",20,-0.5,19.5));
  a->AddPlot(new TH1F("NBJetNumber","Number of Non-b-Jets",20,-0.5,19.5));
  a->AddPlot(new TH2F("NBJetsVsBJets","Non-B-Jets Vs B-Jets; N B-jets; B-jets",20,-0.5,19.5,20,-0.5,19.5));
  a->AddPlot(new TH1F("LeadingIsB","Leading Jet is a B-jet",2,-0.5,1.5));
  a->AddPlot(new TH1F("LeptonPt","Lepton Pt",1000,0,1000));
  a->AddPlot(new TH1F("LeadingJetPt","Leading Jet Pt",1000,0,1000));


  // Additional
  a->AddPlot(new TH1F("WPBDiffLeading","W' b minus other leading parton",400,-1000,1000));
  a->AddPlot(new TH1F("WPBDiffLeadingB","W' b minus other leading b",400,-1000,1000));
  a->AddPlot(new TH1F("WPTBDiffOtB","W'-t-b minus the other t-b",400,-1000,1000));
  a->AddPlot(new TH2F("WPBDiffTB","W' b minus t-b; W'-b - W'-t-b; W'-b - Other t-b",200,-1000,1000,200,-1000,1000));
  a->AddPlot(new TH2F("TPtVsLepPt","Leptonic t Pt Vs lepton Pt;lepton Pt;t Pt",100,0,1000,100,0,1000));
  a->AddPlot(new TH2F("TPtVsMETPt","Leptonic t Pt Vs MET Pt;MET Pt;t Pt",100,0,1000,100,0,1000));
  a->AddPlot(new TH1F("TPtMinusLepPt","Leptonic t Pt minus Lepton Pt",200,0,1000));
  a->AddPlot(new TH1F("LeadingParticleID","Leading Particle Index",6,-0.5,5.5));

}

int hypothesis_loop(Analyzer *a) {
  int SampleType = a->SampleType;
  map<string, TH1F*> p1d = a->Plots1D;
  map<string, TH2F*> p2d = a->Plots2D;

  p1d["EventCounter"]->Fill(0);
  p1d["NumberW"]->Fill(a->GenW.size());
  p1d["Numbert"]->Fill(a->GenT.size());
  p1d["Numberb"]->Fill(a->GenB.size());
  p1d["NumberWP"]->Fill(a->GenWP.size());
  p1d["JetNumber"]->Fill(a->Jets.size());
  p1d["BJetNumber"]->Fill(a->BJets.size());
  p1d["NBJetNumber"]->Fill(a->NBJets.size());
  p2d["NBJetsVsBJets"]->Fill(a->NBJets.size(),a->BJets.size());

  if (a->Jets.size() == 0) return -1;
  if (a->LVLeptons.size() != 1) {
    p1d["EventCounter"]->Fill(1);
    return -1;
  }

  if (a->Jets.size() >= 5) {
    p1d["EventCounter"]->Fill(2);
    if (a->BJets.size() >=2 ) p1d["EventCounter"]->Fill(3);
  }

  if (a->GenPass < 0) return -1;

  a->AssignGenParticles();

  //Tops
  p1d["TopdR"]->Fill(a->LVGenHadT.DeltaR(a->LVGenLepT));

  if (SampleType != 2) {
    p1d["WPTPt"]->Fill(a->LVGenWPT.Pt());
    p1d["OtTPt"]->Fill(a->LVGenOtT.Pt());
    p1d["TopPtDiff"]->Fill(a->LVGenWPT.Pt() - a->LVGenOtT.Pt());
  }
  else {
    p1d["WPTPt"]->Fill(a->LVGenHadT.Pt());
    p1d["OtTPt"]->Fill(a->LVGenLepT.Pt());
    p1d["TopPtDiff"]->Fill(a->LVGenHadT.Pt() - a->LVGenLepT.Pt());
  }

  //Ws
  p1d["WdR"]->Fill(a->LVGenHadW.DeltaR(a->LVGenLepW));
  p1d["WPtDiff"]->Fill(fabs(a->LVGenHadW.Pt() - a->LVGenLepW.Pt()));

  //observables
  p1d["LFPt"]->Fill(a->LVGenLFJet.at(0).Pt());
  p1d["LFPt"]->Fill(a->LVGenLFJet.at(1).Pt());
  if (SampleType != 2) {
    p1d["WPBPt"]->Fill(a->LVGenWPB.Pt());
    p1d["WPTBPt"]->Fill(a->LVGenWPTB.Pt());
    p1d["OtTBPt"]->Fill(a->LVGenOtB.Pt());
    bool leading = (a->LVGenWPB.Pt() > a->LVGenHadB.Pt() && a->LVGenWPB.Pt() > a->LVGenLepB.Pt() && a->LVGenWPB.Pt() > a->LVGenLFJet.at(0).Pt() && a->LVGenWPB.Pt() > a->LVGenLFJet.at(1).Pt() );
    if (leading) p1d["WPBLeading"]->Fill(1);
    else p1d["WPBLeading"]->Fill(0);
  }
  else {
    p1d["WPTBPt"]->Fill(a->LVGenHadB.Pt());
    p1d["OtTBPt"]->Fill(a->LVGenLepB.Pt());
  }

  map<double,int> JetOrder;
  for (unsigned it = 0; it < a->Jets.size(); ++it) {
    JetOrder.insert(pair<double,int>(a->Jets.at(it)->PT,it));
  }
  int nit = 0;
  for (auto it = JetOrder.rbegin(); it != JetOrder.rend(); ++it) {
    if (nit != (*it).second) a->logfile << "Jetorder was wrong at " << nit << " th jet of Event:" << a->iEntry <<endl;
    ++nit;
  }
  p1d["LeadingIsB"]->Fill(a->Jets.at((*(JetOrder.rbegin())).second)->BTag);

  p1d["LeptonPt"]->Fill(a->LVLeptons[0].Pt());

  //Additional
  vector<double> ptout;
  for (unsigned iout = 0; iout < a->LVGenOutSort.size(); ++iout) {
    ptout.push_back(a->LVGenOutSort.at(iout).Pt());
  }
  double wpbpt = a->LVGenWPB.Pt();
  double leadingpt = *max_element(ptout.begin(),ptout.begin()+4);
  double leadingbpt = *max_element(ptout.begin()+2,ptout.begin()+4);
  p1d["LeadingJetPt"]->Fill(*max_element(ptout.begin(),ptout.end()));
  p1d["LeadingParticleID"]->Fill(int(max_element(ptout.begin(),ptout.end()) - ptout.begin()));

  if (SampleType != 2) {
    p1d["WPBDiffLeading"]->Fill(wpbpt - leadingpt);
    p1d["WPBDiffLeadingB"]->Fill(wpbpt - leadingbpt);
    p1d["WPTBDiffOtB"]->Fill(a->LVGenWPTB.Pt() - a->LVGenOtB.Pt());
    p2d["WPBDiffTB"]->Fill(wpbpt - a->LVGenWPTB.Pt(), wpbpt - a->LVGenOtB.Pt());
  }
  else {
    p1d["WPTBDiffOtB"]->Fill(a->LVGenHadB.Pt() - a->LVGenLepB.Pt());
  }
  p2d["TPtVsLepPt"]->Fill(a->LVGenLep.Pt(), a->LVGenLepT.Pt());
  p2d["TPtVsMETPt"]->Fill(a->LVMET.Pt(), a->LVGenLepT.Pt());
  p1d["TPtMinusLepPt"]->Fill(a->LVGenLepT.Pt() - a->LVGenLep.Pt());

  return 0;
}

void hypothesis_end(Analyzer *a) {

}
