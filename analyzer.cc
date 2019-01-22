#ifndef ANALYZER_CC
#define ANALYZER_CC

#include "analyzer.hh"

Analyzer::Analyzer(vector<TString> basepaths, vector<TString> inputfolders, double pt) {
  // gSystem->Load("libDelphes");
  verbose = false;
  vector<TString> fnames = SearchFiles(basepaths, inputfolders, false);
  chain_ = new TChain("Delphes");
  for(unsigned fIdx=0; fIdx<fnames.size(); ++fIdx) chain_->Add(fnames[fIdx]);
  treeReader = new ExRootTreeReader(chain_);
  branchJet = treeReader->UseBranch("Jet");
  branchMuon = treeReader->UseBranch("Muon");
  branchElectron = treeReader->UseBranch("Electron");
  branchMPT = treeReader->UseBranch("MissingET");
  branchGen = treeReader->UseBranch("Particle");
  nEntries = treeReader->GetEntries();
  cout << "Total Number of Events in this dataset is: " << nEntries << endl;
  StartEntry = 0;
  EndEntry = nEntries;
  JetPtThreshold = pt;

  double etaarr[] = {0.,1.3,2.5,3.0,5.2};
  etabins.clear();
  for (unsigned i = 0; i < (sizeof(etaarr)/sizeof(etaarr[0])); ++i) {
    etabins.push_back(etaarr[i]);
  }

  double pttemp[4][23] = {
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130.,150.,180.,220.,260.,300.,350.,400.,500.,1000.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,300.,6000.,0,0,0,0,0},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,6000.,0,0,0,0,0,0},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.,0,0,0,0,0,0,0,0,0}
  };
  ptbins.clear();
  npt.clear();
  for (unsigned ieta = 0 ; ieta < etabins.size() - 1; ++ieta){
    npt.push_back(sizeof(pttemp[ieta])/sizeof(pttemp[ieta][0]));
    vector<double> ptbineta;
    for (int ipt = 0; ipt < npt[ieta] ; ++ipt) {
      if (pttemp[ieta][ipt] == 0. ) {npt[ieta] = ipt;break;};
      ptbineta.push_back(pttemp[ieta][ipt]);
    }
    ptbins.push_back(ptbineta);
  }

}

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
}

void Analyzer::ProcessEntries(int nprocess){
  if (nprocess > 0) EndEntry = StartEntry + nprocess;
  if (EndEntry > nEntries) {
    cout << "Too many ProcessEntries. EndEntry is set to end of sample." << endl;
    EndEntry = nEntries;
  }
}
void Analyzer::SetEndEntry(int endentry) {
  EndEntry = endentry;
}
void Analyzer::SetOutput(TString outputfolder, TString outputname) {
  OutputName = outputfolder+outputname;
  OutputLogName = outputfolder+"logs/"+ outputname;
  TString ofilename = OutputName+".root";
  TString logname = OutputLogName+".txt";
  ofile = new TFile(ofilename,"RECREATE");
  logfile.open(logname);
  logfile << "Log Started" <<endl;
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
//Event preselection criteria are in here.
int Analyzer::ReadEvent(Int_t ievt, bool debug) {
  if (ievt == StartEntry) cout << Form("Processing Events from: %lli to %lli :", StartEntry, EndEntry - 1) << endl;
  iEntry = ievt;
  treeReader->ReadEntry(ievt);
  int dividen = 0;
  if (debug) dividen = 1;
  else dividen = 100;
  printprogress(ievt,EndEntry,dividen);
  Analyzer::GetInfos();

  // Preselection
  if (GenWP.size() != 1 || LVLeptons.size() != 1 || Jets.size() < 5) return -1;
  return 0;
}

void Analyzer::GetInfos() {
  SortInitialize();

  for (int iGen = 0; iGen < branchGen->GetEntries(); ++iGen){
    GenParticle* genp = (GenParticle*) branchGen->At(iGen);
    GenParticles.push_back(genp);
  }
  Analyzer::SortGenParticle();

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

  for (unsigned ihad = 0; ihad < GenOut.size(); ++ihad) {
    LVGenJets.push_back(GenParticles[ToBeHadron(GenOut[ihad])]->P4());
  }
  AllJetMatchMap = JetMatch(LVGenJets, LVJets, AllJetMatchMaxDeltaR, false);

  for (int it = 0; it < branchMuon->GetEntries(); ++it){
    Muon* mu = (Muon*) branchMuon->At(it);
    Muons.push_back(mu);
    LVMuons.push_back(mu->P4());
    if (mu->IsolationVar > 0.2 || mu->PT < 30) {
      LVSoftLep.push_back(mu->P4());
      continue;
    }
    LVLeptons.push_back(mu->P4());
  }

  for (int it = 0; it < branchElectron->GetEntries(); ++it){
    Electron* e = (Electron*) branchElectron->At(it);
    Electrons.push_back(e);
    LVElectrons.push_back(e->P4());
    if (e->IsolationVar > 0.2 || e->PT < 30) {
      LVSoftLep.push_back(e->P4());
      continue;
    }
    LVLeptons.push_back(e->P4());
  }
  if (LVLeptons.size() == 1) {
    MET = (MissingET*) branchMPT->At(0);
    LVMET = MET->P4();
  }
}

void Analyzer::SortInitialize() {
  GenParticles.clear();
  Jets.clear();
  Muons.clear();
  Electrons.clear();
  BTags.clear();
  BJets.clear();
  NBJets.clear();

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
  GenOut.clear();
  GenOutSort.clear();
  LVGenWPTWJ.clear();

  LVJets.clear();
  LVBJets.clear();
  LVNBJets.clear();
  LVElectrons.clear();
  LVMuons.clear();
  LVLeptons.clear();
  LVSoftLep.clear();
  LVGenJets.clear();
  LVGenOutSort.clear();
  AllJetMatchMap.clear();
  OutJetMatchMap.clear();
  LVJetSort.clear();
}

void Analyzer::SortGenParticle() {
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
      // else if (pid == 11) GenE.push_back(igen);
      // else if (pid == 12) GenNuE.push_back(igen);
      // else if (pid == 13) GenMu.push_back(igen);
      // else if (pid == 14) GenNuMu.push_back(igen);
      // else if (pid == 15) GenTau.push_back(igen);
      // else if (pid == 21) GenG.push_back(igen);
      // else if (pid == 22) GenGamma.push_back(igen);
      else if (pid == 24) GenW.push_back(igen);
      else if (pid == 34) GenWP.push_back(igen);
      if ( (genp->Status == 23 || genp->Status ==24) &&  pid< 7 /*&& genp->PT > 20. && fabs(genp->Eta) < 5.*/) {
        GenOut.push_back(igen);
      }
    }
  }
}

void Analyzer::AssignGenParticles() {
  WP = GenWP[0];
  GenWPB = WeHadAMother(GenB, WP).at(0);
  GenWPT = WeHadAMother(GenT, WP).at(0);
  GenWPTB = WeHadAMother(GenB, GenWPT).at(0);
  GenWPTW = WeHadAMother(GenW, GenWPT).at(0);
  GenWPTWJ = WeHadAMother(GenOut, GenWPTW);
  for (unsigned it = 0; it < GenW.size(); ++it) {
    if (GenW[it] != GenWPTW) {
      GenOtW = GenW[it];
      int GenOtWM = BeRealMother(GenOtW);
      int D1 = GenParticles[GenOtWM]->D1;
      int D2 = GenParticles[GenOtWM]->D2;
      GenParticle* d1 = GenParticles[D1];
      GenParticle* d2 = GenParticles[D2];
      if ((abs(d1->PID) == 11 || abs(d1->PID) == 13) && (abs(d2->PID) == 12 || abs(d2->PID) == 14)) {GenLep = D1; GenNeu = D2;}
      else if ((abs(d2->PID) == 11 || abs(d2->PID) == 13) && (abs(d1->PID) == 12 || abs(d1->PID) == 14)) {GenLep = D2; GenNeu = D1;}
      else if (verbose){
        cout << "LepW found wrong one!!!  ";
        cout << "LepW :" << GenOtW << "HadW:"<< GenWPTW << endl;
        PrintGenParticles("LepW");
        continue;}
      break;
    }
  }
  GenOtT = WeHadADaughter(GenT, GenOtW).at(0);
  GenOtB = WeHadAMother(GenB,GenOtT).at(0);

  LVWP = GenParticles[WP]->P4();
  LVGenWPB = GenParticles[GenWPB]->P4();
  LVGenWPT = GenParticles[GenWPT]->P4();
  LVGenWPTB = GenParticles[GenWPTB]->P4();
  LVGenWPTW = GenParticles[GenWPTW]->P4();
  LVGenOtW = GenParticles[GenOtW]->P4();
  LVGenOtT = GenParticles[GenOtT]->P4();
  LVGenOtB = GenParticles[GenOtB]->P4();
  LVGenLep = GenParticles[GenLep]->P4();
  LVGenNeu = GenParticles[GenNeu]->P4();
  for (unsigned ji = 0; ji < GenWPTWJ.size(); ++ji) {
    if (abs(GenParticles[GenWPTWJ[ji]]->PID) < 7) {
      LVGenWPTWJ.push_back(GenParticles[GenWPTWJ[ji]]->P4());
      GenOutSort.push_back(GenWPTWJ[ji]);
    }
  }
  if (LVGenWPTWJ.size() > 2) cout <<iEntry << " WPTWJ more than 2" <<endl;
  LVGenWPTWJ1 = LVGenWPTWJ[0];
  LVGenWPTWJ2 = LVGenWPTWJ[1];

  //Verification of outgoing particles
  int temp = WeHadAMother(GenOut,GenWPT).at(0);
  if (GenWPTB != temp) cout <<iEntry << " WPTB found wrong"<<endl;
  GenOutSort.push_back(temp);

  temp = WeHadAMother(GenOut,GenOtT).at(0);
  if (GenOtB != temp) cout <<iEntry << " OtB found wrong"<<endl;
  GenOutSort.push_back(temp);

  temp = WeHadAMother(GenOut,WP).at(0);
  if (GenWPB != temp) cout <<iEntry << " WPB found wrong"<<endl;
  GenOutSort.push_back(temp);

  for (unsigned i = 0; i < GenOutSort.size(); ++i){
    LVGenOutSort.push_back(GenParticles[ToBeHadron(GenOutSort[i])]->P4());
  }
  OutJetMatchMap = JetMatch(LVGenOutSort,LVJets,OutJetMatchDeltaR,false);
  for (map<int,int>::iterator it = OutJetMatchMap.begin(); it != OutJetMatchMap.end(); ++it) {
    LVJetSort.push_back(LVJets[it->second]);
  }
}

void Analyzer::SaveOutput() {
  ofile->Write();
  ofile->Save();
}


TString Analyzer::GenParticleInfo(int iGen) {
  GenParticle* genp = GenParticles[iGen];
  TLorentzVector p4 = genp->P4();
  TString out = Form(" %7d %7d %6d %5d %5d %5d %5d %7.2f %7.2f    (%8.2f,%8.2f,%8.2f,%8.2f)",iGen, genp->PID, genp->Status, genp->M1, genp->M2, genp->D1, genp->D2, genp->PT, genp->Mass, p4.X(), p4.Y(), p4.Z(), p4.T() );
  return out;
}

void Analyzer::PrintGenParticles(string type) {
  TString suffix;
  if (type == "") suffix = Form("_Evt%lli.txt",iEntry);
  else suffix = Form("_%s_Evt%lli.txt",type.c_str(),iEntry);
  ofstream evtlog(OutputLogName+suffix);
  evtlog << "Entry  GenPar     PID Status    M1    M2    D1    D2      PT    Mass    (P4)" << endl;
  for (unsigned it = 0; it < GenParticles.size(); ++ it ) {
    if (GenParticles[it]->Status < 30 && GenParticles[it]->Status > 20) evtlog << " *** ";
    else evtlog << "     ";
    evtlog << Analyzer::GenParticleInfo(it) << endl;
    // evtlog << Form("%5d", (int) iEntry) << Analyzer::GenParticleInfo(it) << endl;
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

GenParticle* Analyzer::GetGenParticle(int igen) {
  return GenParticles[igen];
}
Jet* Analyzer::GetJet(int ijet) {
  return Jets[ijet];
}
Muon* Analyzer::GetMuon(int imuon) {
  return Muons[imuon];
}
Electron* Analyzer::GetElectron(int ie) {
  return Electrons[ie];
}

TLorentzVector Analyzer::RecoParticle(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, double target, int &p1, int &p2) {
  p1 = 0;
  p2 = (v1 == v2 ? 1 : 0);
  TLorentzVector out = v1[p1] + v2[p2];
  for (unsigned iv1 = 0; iv1 < v1.size(); ++iv1) {
    for (unsigned iv2 = (v1 == v2 ? iv1 + 1 : 0); iv2 < v2.size(); ++iv2) {
      TLorentzVector add = v1.at(iv1) + v2.at(iv2);
      if (fabs(add.M() - target) < fabs(out.M() - target) ) {
        out = add;
        p1 = iv1;
        p2 = iv2;
      }
    }
  }
  return out;
}

map<int,int> Analyzer::JetMatch(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, double &maxDeltaR, bool cut) {
  map<int,int> out;
  vector< vector<double> > deltaRs;
  for (unsigned iv1 = 0; iv1 < v1.size(); ++iv1) {
    vector<double> v2p;
    for (unsigned iv2 = 0; iv2 < v2.size(); ++iv2) {
      v2p.push_back(v1[iv1].DeltaR(v2[iv2]));
    }
    deltaRs.push_back(v2p);
  }
  int x = min(v1.size(), v2.size());
  for (int it = 0; it < x; ++it) {
    int minv1(0), minv2(0);
    double minR = deltaRs[0][0];
    for (unsigned iv1 = 0;iv1 < v1.size(); ++iv1) {
      for (unsigned iv2 = 0; iv2 < v2.size(); ++iv2) {
        if (deltaRs[iv1][iv2] < minR) {
          minR = deltaRs[iv1][iv2];
          minv1 = iv1;
          minv2 = iv2;
        }
      }
    }
    if (it > x - 2) maxDeltaR = minR;
    if (cut && minR > 0.2) break;
    for (unsigned iv1 = 0; iv1 < v1.size(); ++iv1) deltaRs[iv1][minv2] = 1000;
    for (unsigned iv2 = 0; iv2 < v2.size(); ++iv2) deltaRs[minv1][iv2] = 1000;
    out.insert(pair<int, int>(minv1, minv2));
  }
  return out;
}

void Analyzer::GetProbability(bool ForceRecreate) {
  TString pfilename = OutputName+"_Probability.root";
  if (!fileexists(pfilename) || ForceRecreate) CreateProbability(pfilename);
  else {
    PFile = new TFile(pfilename,"READ");
    ofile->cd();
  //   pJet.clear();
  //   for (unsigned ieta = 0; ieta < etabins.size() - 1; ++ieta) {
  //     TString pjetname = Form("pJet_%d",ieta);
  //     pJet.push_back((TH1F*)PFile->Get(pjetname));
  //   }
  //   pGenTMass = (TH1F*) PFile->Get("pGenTMass");
  //   pGenWMass = (TH1F*) PFile->Get("pGenWMass");
  }
}

void Analyzer::CreateProbability(TString pfilename) {
  cout << "Creating Probability Histrograms" <<endl;
  PFile = new TFile(pfilename,"RECREATE");
  vector<TH2F*> hJet, pJet2D;
  for (unsigned ieta = 0; ieta < etabins.size() - 1; ++ieta){
    int sizz = ptbins[ieta].size();
    double ptarr[sizz];
    for (unsigned i = 0; i < ptbins[ieta].size(); ++i) {
      ptarr[i] = ptbins[ieta][i];
    }
    // if (ieta == 0) {
    //   double temp[] = {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130.,150.,180.,220.,260.,300.,350.,400.,500.,1000.,6000.};
    //   ptarr = temp;
    // }
    // if (ieta == 1) {
    //   double temp[] = {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,300.,6000.};
    //   ptarr = temp;
    // }
    // if (ieta == 2) {
    //   double temp[] = {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,6000.};
    //   ptarr = temp;
    // }
    // if (ieta == 3) {
    //   double temp[] = {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.};
    //   ptarr = temp;
    // }

    TString pjetname = Form("pJet_%d",ieta);
    TString pjettitle = Form("pJet_%3.1fto%3.1f ; RecoPt ; RecoPt/GenPt",etabins[ieta],etabins[ieta+1]);
    TString pjetname2D = Form("pJet_%d2D",ieta);
    TString pjettitle2D = Form("pJet_%3.1fto%3.1f2D ; RecoPt ; RecoPt/GenPt",etabins[ieta],etabins[ieta+1]);
    TString hjetname = Form("hJet_%d",ieta);
    TString hjettitle = Form("hJet_%3.1fto%3.1f ; RecoPt / GenPt ; DeltaR",etabins[ieta],etabins[ieta+1]);

    pJet.push_back(new TH1F(pjetname,pjettitle, npt[ieta] - 1 , ptarr));
    pJet2D.push_back(new TH2F(pjetname2D,pjettitle2D, npt[ieta] - 1, ptarr,600,0 , 6.));
    hJet.push_back(new TH2F(hjetname,hjettitle,600, 0, 6., 50, 0,5.));
  }
  pGenTMass = new TH1F("pGenTMass","GenTMass",350,0.,350.);
  pGenWMass = new TH1F("pGenWMass","GenWMAss",160,0.,160.);
  pGenWPdPhi = new TH1F("pGenWPdPhi","GenWPdPhi",40,0.,4.);

  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    if (ReadEvent(ievt,false) == -1) continue;
    AssignGenParticles();
    //Jet Resolution;
    map<int,int> match = AllJetMatchMap;
    for (map<int,int>::iterator it = match.begin(); it != match.end(); ++it) {
      TLorentzVector LVGen, LVReco;
      LVGen = LVGenJets[it->first];
      LVReco = LVJets[it->second];
      double deltaR = LVReco.DeltaR(LVGen);
      for (unsigned ieta = 0; ieta < etabins.size() - 1; ++ieta) {
        if (fabs(LVGen.Eta()) < etabins[ieta + 1]) {
          double ratio = LVGen.Pt()/LVReco.Pt();
          hJet[ieta]->Fill(ratio, deltaR);
          if (deltaR > 0.1) continue;
          pJet2D[ieta]->Fill(LVReco.Pt(),ratio);
          break;
        }
      }
    }
    //W and T mass Resolution
    for (unsigned i = 0; i < GenT.size(); ++i) {
      pGenTMass->Fill(GenParticles[GenT[i]]->Mass);
    }
    for (unsigned i = 0; i < GenW.size(); ++i) {
      pGenWMass->Fill(GenParticles[GenW[i]]->Mass);
    }
    pGenWPdPhi->Fill(LVGenWPB.DeltaPhi(LVGenWPT));
  }

  for (unsigned ieta = 0; ieta < etabins.size() - 1; ++ieta) {
    for (int ibin = 1; ibin < npt[ieta]; ++ibin) {
      TF1* f = new TF1("fit","gaus",0.7,1.3);
      TH1D* pj = pJet2D[ieta]->ProjectionY("_py",ibin,ibin);
      if (pj->GetEntries() < 10) continue;
      f->SetParameters(pj->Integral(), 1.0, 0.3);
      pj->Fit(f,"RFQ");
      double p1 = f->GetParameter(1);
      double p2 = f->GetParameter(2);
      pJet[ieta]->SetBinContent(ibin,p1);
      pJet[ieta]->SetBinError(ibin,p2);
    }
  }

  pGenTMass->Scale(1./pGenTMass->GetMaximum());
  pGenWMass->Scale(1./pGenWMass->GetMaximum());
  pGenWPdPhi->Scale(1./pGenWPdPhi->GetMaximum());

  PFile->Write();
  PFile->Save();
  ofile->cd();
  cout << "Probability Histograms done. Starting Analysis." <<endl;
}

double Analyzer::Optimize() {
  return 0;
  // return Optimize(OptiBestPerm, OptiScaledjets, OptiNeutrino, OptiWPrime);
}

double Analyzer::Optimize(vector<int> &BestPerm, vector<TLorentzVector> &scaledjets, TLorentzVector &neutrino, TLorentzVector &WPrime, double &pjes, double &pmass, double &ptype) {
  Optimizer *o = new Optimizer(LVJets, BTags);
  o->SetLepton(LVLeptons[0]);
  o->SetMET(LVMET);
  o->SetPFile(PFile);
  o->SetEtaBins(etabins);
  o->Optimize();

  OptiBestPerm.clear();
  OptiScaledjets.clear();
  OptiP = 0;
  OptiP = o->GetBestP();
  if (OptiP <= 0) {
    logfile <<"This event:"<<iEntry<< " Does not have neutrino solution, Lepton Pt: " << LVLeptons[0].Pt()<<endl;
    return -1;
  }
  BestPerm = o->GetBestPerm();
  scaledjets = o->GetBestLVJets();
  neutrino = o->GetBestNeutrino();
  WPrime = o->GetWPrime();

  OptiBestPerm = BestPerm;
  OptiScaledjets = scaledjets;
  OptiNeutrino = neutrino;
  OptiWPrime = WPrime;
  pjes = o->GetBestPJES();
  pmass = o->GetBestPMass();
  ptype = o->GetBestPType();

  // cout << Form("BestP: %f, scaledjets: %i, neutrino:(%f,%f,%f,%f),WPrimeMass: %f",p, scaledjets.size(), neutrino.X(),neutrino.Y(),neutrino.Z(),neutrino.T(),WPrime.M())<<endl;
  return OptiP;
}

int Analyzer::OptiJetMatch(double &dR) {
  int MatchedJets = 0;
  int mwptwj(0), mwptb(0), mlepb(0), mwpb(0);
  OptiJetMatchMap = JetMatch(LVGenOutSort, OptiScaledjets, dR, false);

  if (OptiJetMatchMap[0] == 0 || OptiJetMatchMap[0] == 1) mwptwj++;
  if (OptiJetMatchMap[1] == 0 || OptiJetMatchMap[1] == 1) mwptwj++;
  if (OptiJetMatchMap[2] == 2) mwptb = 1;
  if (OptiJetMatchMap[3] == 3) mlepb = 1;
  if (OptiJetMatchMap[4] == 4) mwpb = 1;

  MatchedJets = mwptwj + mwptb + mlepb + mwpb;

  return MatchedJets;

}

#endif
