#ifndef ANALYZER_CC
#define ANALYZER_CC

#include "Analyzer.hh"

Analyzer::Analyzer(int SampleType_, int irun_, double pt = 30) {
  // gSystem->Load("libDelphes");
  verbose = false;
  SampleType = SampleType_;
  irun = irun_;
  vector<TString> fnames = SearchFiles(SampleType, irun, true);
  chain_ = new TChain("Delphes");
  for(unsigned fIdx=0; fIdx<fnames.size(); ++fIdx) chain_->Add(fnames[fIdx]);
  treeReader = new ExRootTreeReader(chain_);
  branchJet = treeReader->UseBranch("Jet");
  branchMuon = treeReader->UseBranch("Muon");
  branchElectron = treeReader->UseBranch("Electron");
  branchMPT = treeReader->UseBranch("MissingET");
  branchGen = treeReader->UseBranch("Particle");
  branchGenJet = treeReader->UseBranch("GenJet");
  nEntries = treeReader->GetEntries();
  cout << "Total Number of Events in this dataset is: " << nEntries << endl;
  StartEntry = 0;
  EndEntry = nEntries;
  JetPtThreshold = pt;
  LepIso = 0.2;
  debug = -2;
  RecoPass = 0;
  GenPass = 0;

};

void Analyzer::SetStartEntry(int startentry){
  StartEntry = startentry;
  if (StartEntry > nEntries - 1 || StartEntry > EndEntry -1 ) {
    cout << "Invalid StartEntry. StartEntry is set to end of sample and code terminating." <<endl;
    StartEntry = nEntries;
  }
  if (StartEntry < 0) {
    cout << "StartEntry Cannot be set to negative. StartEntry is set to 0." <<endl;
    StartEntry = 0;
  }

};

void Analyzer::ProcessEntries(int nprocess){
  if (nprocess > 0) EndEntry = StartEntry + nprocess;
  if (EndEntry > nEntries) {
    cout << "Too many ProcessEntries. EndEntry is set to end of sample." << endl;
    EndEntry = nEntries;
  }

};


void Analyzer::SetEndEntry(int endentry) {
  EndEntry = endentry;
}

void Analyzer::SetOutput(TString outputfolder_, TString outputname_) {
  if      (SampleType == 0) outputname_ += "_FL";
  else if (SampleType == 1) outputname_ += "_LL";
  else                      outputname_ += "_BG";
  if (irun != 0) outputname_ += Form("_%i",irun);
  if (irun != 0) outputfolder_ += "massresult/";
  outputname = outputname_;
  outputfolder = outputfolder_;
  OutputName = outputfolder+outputname;
  OutputLogName = outputfolder+"logs/"+ outputname;
  TString ofilename = OutputName+".root";
  TString logname = OutputLogName+".txt";
  ofile = new TFile(ofilename,"RECREATE");
  AddPlot(new TH1F("EventCounter","Number of Events;Passcode",100,0,100));
  logfile.open(logname);
  logfile << "Log Started" <<endl;
}

void Analyzer::CDOut() {
  ofile->cd();
}

void Analyzer::SetJetPtThreshold(double pt) {
  JetPtThreshold = pt;
}

Long64_t Analyzer::GetEntries() {
  return nEntries;
}

Long64_t Analyzer::GetStartEntry() {
  return StartEntry;
}

Long64_t Analyzer::GetEndEntry() {
  return EndEntry;
}

void Analyzer::DebugMode(int debug_ = -2) {
  debug = debug_;
  if (debug == -1) {
    SetStartEntry(0);
    ProcessEntries(1000);
  }
  else if (debug > -1){
    SetStartEntry(debug);
    ProcessEntries(1);
  }
}

void Analyzer::SaveOutput() {
  Plots1D["EventCounter"]->GetXaxis()->SetRange(1,CounterLabels.size()+1);
  ofile->Write();
  ofile->Save();
}

// ------------------------------------------------------------------------------------------------

// Initialize Event
void Analyzer::ReadEvent(Int_t ievt) {
  GenPass = 0;
  RecoPass = 0;
  if (ievt == StartEntry) cout << Form("Processing Events from: %lli to %lli :", StartEntry, EndEntry - 1) << endl;
  iEntry = ievt;
  treeReader->ReadEntry(ievt);
  int dividen = 1;
  if (debug == -2) dividen = 100;
  PrintProgress(ievt,EndEntry,dividen);
  Analyzer::GetInfos();

  // Preselection
  if (LVLeptons.size() != 1 || Jets.size() < 5) RecoPass = -1;
  if ((GenWP.size() != 1 && SampleType != 2) || GenW.size() != 2) GenPass = -1;

  CountEvent("All_Evts");
}

void Analyzer::CountEvent(TString lbl) {
  Plots1D["EventCounter"]->Fill(lbl,1);
  if (iEntry == StartEntry) {
    CounterLabels.insert(lbl);
  }
}

void Analyzer::GetInfos() {
  GenParticles.clear();
  for (int iGen = 0; iGen < branchGen->GetEntries(); ++iGen){
    GenParticle* genp = (GenParticle*) branchGen->At(iGen);
    GenParticles.push_back(genp);
  }
  GenParticleTypes();

  LVOutPart.clear();
  for (unsigned ihad = 0; ihad < GenOutQuark.size(); ++ihad) {
    LVOutPart.push_back(GenParticles[ToBeHadron(GenOutQuark[ihad])]->P4());
  }
  for (unsigned igl = 0; igl < GenOutGluon.size(); ++igl) {
    LVOutPart.push_back(GenParticles[ToBeHadron(GenOutGluon[igl])]->P4());
    GenOutBTags.push_back(false);
  }

  AllJets.clear();
  LVAllJets.clear();
  Jets.clear();
  LVJets.clear();
  BTags.clear();
  BJets.clear();
  LVBJets.clear();
  NBJets.clear();
  LVNBJets.clear();
  for (int it = 0; it < branchJet->GetEntries(); ++it){
    Jet* jet = (Jet*) branchJet->At(it);
    AllJets.push_back(jet);
    LVAllJets.push_back(jet->P4());
    if (jet->PT < JetPtThreshold) continue;
    Jets.push_back(jet);
    LVJets.push_back(jet->P4());
    BTags.push_back(jet->BTag);
    if (jet->BTag) {
      BJets.push_back(jet);
      LVBJets.push_back(jet->P4());
    }
    else {
      NBJets.push_back(jet);
      LVNBJets.push_back(jet->P4());
    }
  }
  GenJets.clear();
  LVGenJets.clear();
  for (int it = 0; it < branchGenJet->GetEntries(); ++it) {
    Jet* genjet = (Jet*) branchGenJet->At(it);
    GenJets.push_back(genjet);
    LVGenJets.push_back(genjet->P4());
  }
  Muons.clear();
  LVMuons.clear();
  LVSoftLep.clear();
  LVLeptons.clear();
  for (int it = 0; it < branchMuon->GetEntries(); ++it) {
    Muon* mu = (Muon*) branchMuon->At(it);
    Muons.push_back(mu);
    LVMuons.push_back(mu->P4());
    if (mu->IsolationVar > LepIso || mu->PT < 30) {
      LVSoftLep.push_back(mu->P4());
      continue;
    }
    LVLeptons.push_back(mu->P4());
  }

  Electrons.clear();
  LVElectrons.clear();
  for (int it = 0; it < branchElectron->GetEntries(); ++it){
    Electron* e = (Electron*) branchElectron->At(it);
    Electrons.push_back(e);
    LVElectrons.push_back(e->P4());
    if (e->IsolationVar > LepIso || e->PT < 30) {
      LVSoftLep.push_back(e->P4());
      continue;
    }
    LVLeptons.push_back(e->P4());
  }
  MET = (MissingET*) branchMPT->At(0);
  LVMET = MET->P4();
}

void Analyzer::GenParticleTypes() {
  GenD.clear();
  GenU.clear();
  GenS.clear();
  GenC.clear();
  GenB.clear();
  GenT.clear();
  GenE.clear();
  GenNuE.clear();
  GenMu.clear();
  GenNuMu.clear();
  GenTau.clear();
  GenG.clear();
  GenGamma.clear();
  GenW.clear();
  GenWP.clear();
  GenOutQuark.clear();
  GenOutGluon.clear();
  GenOutBTags.clear();
  for (int igen = 0; igen < branchGen->GetEntries(); ++igen) {
    GenParticle* genp = GenParticles[igen];
    if (genp->Status > 20 && genp->Status < 30) {
      int pid = abs(genp->PID);
      //Leptons are not with status of 2x, status =1 instead.
      if      (pid == 0) continue;
      // else if (pid == 1) GenD.push_back(igen);
      // else if (pid == 2) GenU.push_back(igen);
      // else if (pid == 3) GenS.push_back(igen);
      // else if (pid == 4) GenC.push_back(igen);
      else if (pid == 5) GenB.push_back(igen);
      else if (pid == 6) GenT.push_back(igen);
      else if (pid == 11) GenE.push_back(igen);
      else if (pid == 12) GenNuE.push_back(igen);
      else if (pid == 13) GenMu.push_back(igen);
      else if (pid == 14) GenNuMu.push_back(igen);
      // else if (pid == 15) GenTau.push_back(igen);
      // else if (pid == 21) GenG.push_back(igen);
      // else if (pid == 22) GenGamma.push_back(igen);
      else if (pid == 24) GenW.push_back(igen);
      else if (pid == 34) GenWP.push_back(igen);
      if ( (genp->Status == 23 || genp->Status ==24) ) {
        if (pid <7 ) {
          GenOutQuark.push_back(igen);
          if (pid == 5) GenOutBTags.push_back(true);
          else GenOutBTags.push_back(false);
        }
        if (pid == 21) {
          GenOutGluon.push_back(igen);
        }
      }
    }
  }
}

TString Analyzer::GenParticleInfo(int iGen) {
  GenParticle* genp = GenParticles[iGen];
  TLorentzVector p4 = genp->P4();
  TString out = Form(" %7d %7d %6d M(%5d, %5d) D(%5d, %5d) Pt: %7.2f M: %7.2f    (%8.2f,%8.2f,%8.2f,%8.2f)",iGen, genp->PID, genp->Status, genp->M1, genp->M2, genp->D1, genp->D2, genp->PT, genp->Mass, p4.X(), p4.Y(), p4.Z(), p4.T() );
  return out;
}

void Analyzer::PrintGenParticles(string type = "", bool verb = false) {
  TString suffix;
  if (type == "") suffix = Form("_Evt%lli.txt",iEntry);
  else suffix = Form("_%s_Evt%lli.txt",type.c_str(),iEntry);
  ofstream evtlog(OutputLogName+suffix);
  evtlog << "Entry    GenPar     PID Status      M1     M2       D1     D2          PT       Mass    (P4)" << endl;
  if (verb) {
    for (unsigned it = 0; it < GenParticles.size(); ++ it ) {
      if (GenParticles[it]->Status < 30 && GenParticles[it]->Status > 20) evtlog << "  ***  ";
      else evtlog << "       ";
      evtlog << Analyzer::GenParticleInfo(it) << endl;
      // evtlog << Form("%5d", (int) iEntry) << Analyzer::GenParticleInfo(it) << endl;
    }
  }
  else {
    for (unsigned it = 0; it < GenParticles.size(); ++it ){
      if (GenParticles[it]->Status >= 30 || GenParticles[it]->Status <= 20) continue;
      int AsDaughter = BeRealDaughter(it);
      int AsMother = BeRealMother(it);
      evtlog << Form("D_%04i ",it) << GenParticleInfo(AsDaughter) << endl;
      evtlog << Form("**%2s **",GetParName(GenParticles[it]->PID).c_str());
      evtlog << GenParticleInfo(it) << endl;
      evtlog << Form("M_%04i ",it) << GenParticleInfo(AsMother) << endl<<endl;
    }
  }
  evtlog <<endl<<endl<< "Summary:"<<endl;
  evtlog <<Form("Number of Particles: t:%d, b:%d, w:%d, W':%d",int(GenT.size()),int(GenB.size()),int(GenW.size()),int(GenWP.size()))<<endl<<endl<<endl;
  if (!verb) {
    for (unsigned it = 0; it < GenParticles.size(); ++ it ) {
      if (GenParticles[it]->Status < 30 && GenParticles[it]->Status > 20) evtlog << "  ***  ";
      else evtlog << "       ";
      evtlog << Analyzer::GenParticleInfo(it) << endl;
      // evtlog << Form("%5d", (int) iEntry) << Analyzer::GenParticleInfo(it) << endl;
    }
  }
}

int Analyzer::BeRealMother(int igen) {
  GenParticle* genp = GenParticles[igen];
  int d1 = genp->D1;
  int d2 = genp->D2;
  int pid = genp->PID;
  if (d1 > -1) {
    if (GenParticles[d1]->PID == pid) return BeRealMother(d1);
  }
  if (d2 > -1) {
    if (GenParticles[d2]->PID == pid) return BeRealMother(d2);
  }
  return igen;
}

int Analyzer::BeRealDaughter(int igen) {
  GenParticle* genp = GenParticles[igen];
  int m1 = genp->M1;
  int m2 = genp->M2;
  int pid = genp->PID;

  if (m1 > -1) {
    if (GenParticles[m1]->PID == pid) return BeRealDaughter(m1);
  }
  if (m2 > -1) {
    if (GenParticles[m2]->PID == pid) return BeRealDaughter(m2);
  }
  return igen;
}

vector<int> Analyzer::WeHadAMother(vector<int> &list, int &mother) {
  vector<int> d;
  for (unsigned it = 0; it < list.size(); ++it) {
    GenParticle* p = GenParticles[list[it]];
    int m1(p->M1), m2(p->M2), m(BeRealMother(mother));
    if (m1 == m || m2 == m) d.push_back(list[it]);
  }
  return d;
}

vector<int> Analyzer::WeHadADaughter(vector<int> &list, int &daughter) {
  vector<int> m;
  for (unsigned it = 0; it < list.size(); ++it) {
    GenParticle* p = GenParticles[BeRealMother(list[it])];
    int d1(p->D1), d2(p->D2), d(BeRealDaughter(daughter));
    if (d1 == d || d2 == d) m.push_back(list[it]);
  }
  return m;
}

vector<int> Analyzer::MyDaughters(int &mother) {
  vector<int> d;
  GenParticle* p = GenParticles[BeRealMother(mother)];
  int d1(p->D1), d2(p->D2);
  if (d1 > -1 && d2 > -1) {
    // Case 2
    if (d1 == d2) {
      cout << Form("Entry : %lli ,GenParticle %i, BeRealMother's daughters equals. Log of this event printed as NRM", iEntry, mother);
      PrintGenParticles("NRM");
    }
    // Case 4
    else if (d1 < d2) {
      for (int id = d1; id <= d2; ++id) {
        d.push_back(id);
      }
    }
    // Case 5
    else if (d2 < d1) {
      d.push_back(d1);
      d.push_back(d2);
    }
  }
  // Case 3
  else if (d1 > -1 && d2 == -1) d.push_back(d1);
  return d;
}

vector<int> Analyzer::MyMothers(int &daughter) {
  vector<int> m;
  GenParticle* p = GenParticles[BeRealDaughter(daughter)];
  int m1(p->M1), m2(p->M2);
  if (m1 > -1 && m2 > -1) {
    // Case 2
    if (m1 == m2) {
      cout << Form("Entry : %lli ,GenParticle %i, BeRealDaughter's mothers equals. Log of this event printed as NRD", iEntry, daughter);
      PrintGenParticles("NRD");
    }
    // Case 4 requires Status 81-86, thus not considered here.
    // Case 5 & Case 6
    else if (m1 != m2) {
      m.push_back(m1);
      m.push_back(m2);
    }
  }
  return m;
}

vector<int> Analyzer::WayToBeHadron(int igen) {
  vector<int> out;
  out.push_back(igen);
  GenParticle* p = GenParticles[igen];
  if (p->D1 == -1) return out;
  GenParticle* d1 = GenParticles[p->D1];
  if ((p->D1 == p->D2 || p->D2 == -1) && d1->PID == p->PID) {
    vector<int> dd = WayToBeHadron(p->D1);
    out.insert(out.end(),dd.begin(),dd.end());
  }
  return out;
}

int Analyzer::ToBeHadron(int igen) {
  GenParticle* p = GenParticles[igen];
  if (p->D1 == -1) return igen;
  GenParticle* d1 = GenParticles[p->D1];
  if ((p->D1 == p->D2 || p->D2 == -1) && d1->PID == p->PID) return ToBeHadron(p->D1);
  return igen;
}

void Analyzer::AddPlot(TH1F* hh_) {
  Plots1D[hh_->GetName()] = hh_;
}

void Analyzer::AddPlot(TH2F* hh_) {
  Plots2D[hh_->GetName()] = hh_;
}

int Analyzer::AssignGenParticles() {
  if (GenW.size() !=2 || LVLeptons.size() !=1) return -1;
  if (SampleType == 2) {
    WP = -1;
    GenWPB = -1;
    Gen.WP = TLorentzVector();
    Gen.WPB = TLorentzVector();
  }
  else {
    WP = GenWP[0];
    GenWPB = WeHadAMother(GenB, WP).at(0);
    Gen.WP = GenParticles[WP]->P4();
    Gen.WPB = GenParticles[GenWPB]->P4();
  }
  for (unsigned it = 0; it < GenW.size(); ++it) {
    int genw = BeRealMother(GenW[it]);
    int D1 = GenParticles[genw]->D1;
    int D2 = GenParticles[genw]->D2;
    GenParticle* d1 = GenParticles[D1];
    GenParticle* d2 = GenParticles[D2];
    if ((abs(d1->PID) == 11 || abs(d1->PID) == 13) && (abs(d2->PID) == 12 || abs(d2->PID) == 14)) {GenLep = D1; GenNeu = D2; GenLepW = GenW[it];}
    else if ((abs(d2->PID) == 11 || abs(d2->PID) == 13) && (abs(d1->PID) == 12 || abs(d1->PID) == 14)) {GenLep = D2; GenNeu = D1; GenLepW = GenW[it];}
    else GenHadW = GenW[it];
  }

  GenHadT = WeHadADaughter(GenT,GenHadW).at(0);
  GenLepT = WeHadADaughter(GenT,GenLepW).at(0);
  GenHadB = WeHadAMother(GenB,GenHadT).at(0);
  GenLepB = WeHadAMother(GenB,GenLepT).at(0);
  GenLFJet = WeHadAMother(GenOutQuark, GenHadW);

  Gen.HadW = GenParticles[GenHadW]->P4();
  Gen.HadT = GenParticles[GenHadT]->P4();
  Gen.HadB = GenParticles[GenHadB]->P4();
  Gen.LepW = GenParticles[GenLepW]->P4();
  Gen.LepT = GenParticles[GenLepT]->P4();
  Gen.LepB = GenParticles[GenLepB]->P4();
  Gen.Lep = GenParticles[GenLep]->P4();
  Gen.Neu = GenParticles[GenNeu]->P4();
  vector<TLorentzVector> LVGenLF = vector<TLorentzVector> {TLorentzVector(), TLorentzVector()};
  GenOutSort = vector<int> {-1,-1};

  map<double,int> LFPt;
  LFPt.clear();
  for (unsigned ji = 0; ji < GenLFJet.size(); ++ji) {
    if (abs(GenParticles[GenLFJet[ji]]->PID) < 7) {
      LFPt.insert(pair<double,int>(GenParticles[GenLFJet[ji]]->PT,GenLFJet[ji]));
    }
  }
  if (LFPt.size() > 2) {
    logfile <<iEntry << " LFJet more than 2" <<endl;
  }
  if (LFPt.size() < 2) {
    logfile <<iEntry << " LFJet less than 2" <<endl;
  }
  auto ji = LFPt.end();
  for (int itn = 0; itn < 2; ++itn) {
    if (ji != LFPt.begin()) {
      --ji;
      int genout = (*ji).second;
      GenOutSort[itn] = genout;
      LVGenLF[itn] = GenParticles[genout]->P4();
    }
    else {
      GenOutSort[itn] = -1;
      LVGenLF[itn] = TLorentzVector();
    }
  }
  Gen.LF0 = LVGenLF[0];
  Gen.LF1 = LVGenLF[1];

  //Verification of outgoing particles
  int temp = WeHadAMother(GenOutQuark,GenHadT).at(0);
  if (GenHadB != temp) logfile <<iEntry << " HadB found wrong"<<endl;
  GenOutSort.push_back(temp);

  temp = WeHadAMother(GenOutQuark,GenLepT).at(0);
  if (GenLepB != temp) logfile <<iEntry << " LepB found wrong"<<endl;
  GenOutSort.push_back(temp);

  if (SampleType != 2) {
    temp = WeHadAMother(GenOutQuark,WP).at(0);
    if (GenWPB != temp) logfile <<iEntry << " WPB found wrong"<<endl;
    GenOutSort.push_back(temp);
  }

  LVGenOutSort.clear();
  for (unsigned i = 0; i < GenOutSort.size(); ++i){
    if (GenOutSort[i] == -1) LVGenOutSort.push_back(TLorentzVector());
    else LVGenOutSort.push_back(GenParticles[ToBeHadron(GenOutSort[i])]->P4());
  }
  // GenOutSort Sequence: Light flavor parton (0,1), Hadronic B(2), Leptonic B(3), W' B(4)

  if (SampleType == 0) {
    GenWPT = GenHadT;
    GenWPTB = GenHadB;
    GenWPTW = GenHadW;
    GenOtT = GenLepT;
    GenOtB = GenLepB;
    GenOtW = GenLepW;

    // LVGenWPT = Gen.HadT;
    // LVGenWPTB = Gen.HadB;
    // LVGenWPTW = Gen.HadW;
    // LVGenOtT = Gen.LepT;
    // LVGenOtB = Gen.LepB;
    // LVGenOtW = Gen.LepW;
  }
  else if (SampleType == 1) {
    GenWPT = GenLepT;
    GenWPTB = GenLepB;
    GenWPTW = GenLepW;
    GenOtT = GenHadT;
    GenOtB = GenHadB;
    GenOtW = GenHadW;

    // LVGenWPT = Gen.LepT;
    // LVGenWPTB = Gen.LepB;
    // LVGenWPTW = Gen.LepW;
    // LVGenOtT = Gen.HadT;
    // LVGenOtB = Gen.HadB;
    // LVGenOtW = Gen.HadW;
  }
  return 0;
}

TLorentzVector Analyzer::GetGenParticleP4(int igen) {
  return GenParticles[igen]->P4();
}

void Analyzer::SetupGenTools() {
  GT1 = new GenTools("PFile", "GenFileL1");
  GT2 = new GenTools("PFile", "GenFileL2");
}

void Analyzer::SetupJESTools() {
  JT = new JESTools();
  JESFile = new TFile("PFile/JESFile.root");
  JT->ReadJESPlots(JESFile);
}

void Analyzer::SetupROOTMini() {
  SetupGenTools();
  SetupJESTools();
  RM = new ROOTMini(JT);
  RM->SetDebug(0);
}

void Analyzer::MatchJets(){
  AdvJetMatchMap.clear();
  AdvJetMatchMap = AdvJetMatch(LVOutPart, LVJets, JetMatchMaxDeltaR, 2.0, true, false, true);
  JetMatchMap.clear();
  JetMatchMap = JetMatch(LVOutPart, LVJets, JetMatchMaxDeltaR);
  // AdvJetMatch(genLV, recoLV, max deltaR in the match, cut at larger than AlgodR = 0.2, if skip all gen LV with pt < 30, if skip all reco pt < 30)
}

void Analyzer::GetGenCorrectPerm() {
  GenCorrectPerm.clear();
  for (unsigned i = 0; i < GenOutSort.size(); ++i) {
    int hypo = GenOutSort.at(i);
    if (hypo == -1) {
      GenCorrectPerm.push_back(-1);
      continue;
    }
    for (unsigned j = 0; j < GenOutQuark.size(); ++j) {
      if (GenOutQuark[j] == hypo) {
        GenCorrectPerm.push_back(j);
        break;
      }
    }
  }
}

void Analyzer::GetRecoCorrectPerm() {
  RecoCorrectPerm.clear();
  GetGenCorrectPerm();
  MatchJets();
  for (unsigned it = 0; it < GenCorrectPerm.size(); ++it) {
    RecoCorrectPerm.push_back(JetMatchMap[GenCorrectPerm.at(it)]);
  }
}

void Analyzer::GetRecoHypothesis(){
  RecoHypothesis.clear();
  MatchedHypo = 0;
  // GenOutSort ordered as LFjet1, LFjet2, HadB, LepB, W'B
  MatchJets();
  vector<int> HypoInOutPart; // vector of LVGenOutSort particle index in LVOutPart

  for (unsigned i1 = 0; i1 < LVGenOutSort.size(); ++i1) {
    TLorentzVector v1 = LVGenOutSort.at(i1);
    int found = 0;
    for (unsigned i2 = 0; i2 < LVOutPart.size(); ++i2) {
      TLorentzVector v2 = LVOutPart.at(i2);
      if (v1 == v2) {
        found = 1;
        HypoInOutPart.push_back(i2);
        break;
      }
    }
    if (found == 0) {
      HypoInOutPart.push_back(-1);
      cout << endl << i1 << "th LVGenOutSort is not found in LVOutPart" <<endl;
    }
  }
  for (unsigned it = 0; it < HypoInOutPart.size(); ++it){
    auto mapit = AdvJetMatchMap.find(HypoInOutPart.at(it));
    if (mapit != AdvJetMatchMap.end()) {
      ++MatchedHypo;
      int match1 = (*mapit).second.at(0);
      int match2 = (*mapit).second.at(1);
      if (match1 == match2) {
        RecoHypothesis.push_back(LVJets.at(match1));
      }
      else {
        RecoHypothesis.push_back(LVJets.at(match1) + LVJets.at(match2));
      }
    }
  }
}

pair<double, vector<TLorentzVector> > Analyzer::SolveTTbar(vector<TLorentzVector> jets_, vector<bool> btags_, vector<int>& bestperm_ ) {
  vector< vector<int> > perms = JT->MakePermutations(jets_.size());
  double BestP = -1.0;
  vector<int> BestPerm;
  vector<TLorentzVector> BestParticles;

  for (unsigned iperm = 0; iperm < perms.size(); ++iperm) {
    vector<int> thisperm = perms.at(iperm);
    double PBTags = JT->CalcPFlavor(thisperm, btags_);
    vector<TLorentzVector> PermJets = JT->GetPermutationLV(thisperm, jets_);
    if (find(PermJets.begin(),PermJets.end(),TLorentzVector()) != PermJets.end()) continue;
    double PJes = RM->MinimizeP(PermJets);
    if (PJes < 0) continue;
    double PPerm = PJes * PBTags;
    if (PPerm > BestP) {
      BestP = PPerm;
      BestPerm = thisperm;
      RM->ReCalcP(BestParticles);
    }
  }
  bestperm_ = BestPerm;
  return pair<double, vector<TLorentzVector> > (BestP, BestParticles);
} // This function is not very repetitive, so it needs to be removed

void Analyzer::Tree_Init(int SaveTreeLevel = 3) {
  // TString ofilename = ofile->GetName();
  // TreeFile = new TFile(outputfolder+"Tree_"+outputname+".root","RECREATE");
  TreeFile = new TFile(outputfolder + "TruthTree.root","RECREATE");
  t = new TTree("t0","Event Tree");
  CDOut();
  bool savegen = false;
  bool savereco = false;
  if (SaveTreeLevel > 0) savegen = true;
  if (SaveTreeLevel > 2) savereco = true;

  Gen.BookBranches(t,"Gen",savegen);
  Reco.BookBranches(t,"Reco",savereco);

}

void Analyzer::Tree_Reco() {
  vector<TLorentzVector> genvec = Gen.Observables();
  vector<TLorentzVector> recovec = vector<TLorentzVector>(5);
  double mdr_;
  map<int,int>recomap = JetMatch(genvec, LVJets, mdr_, 0.4, 0.5);
  for (auto it = recomap.begin(); it != recomap.end(); ++it) {
    recovec.at(it->first) = LVJets.at(it->second);
  }
  for (unsigned ir = 0; ir < 5; ++ir) {
    *(Reco.pObservables().at(ir)) = recovec.at(ir);
  }
  Reco.Lep = LVLeptons[0];
  Reco.Neu = LVMET;
  Reco.Neu.SetZ(Gen.Neu.Z());
  Reco.Calculate(SampleType);
}

void Analyzer::Tree_Fill() {
  t->Fill();
}

void Analyzer::Tree_Save() {
  TreeFile->Write();
  TreeFile->Save();
}

#endif
