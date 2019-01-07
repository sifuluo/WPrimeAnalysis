#ifndef TOOLS_CC
#define TOOLS_CC

#include <TROOT.h>
#include <TString.h>
#include <TVector2.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <TProfile.h>
#include <TEfficiency.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

using namespace std;

bool fileexists (const TString& name) {
  struct stat buffer;
  return ( stat(name, &buffer) == 0);
}

bool fileexists2 (const TString& name) {
  ifstream f(name);
  return f.good();
}

void printprogress (int entry, int entries, int dividen = 100) {
  entries--;
  if (entry % dividen ==0 || entry == entries ) {
    cout << Form("\rProgress: %i / %i ", entry, entries) <<flush;
  }
  if (entry == entries){
    cout <<endl << "Done."<<endl;
  }
}

//Get all sample filenames within given path.
vector<TString> SearchFiles(vector<TString> samplebasepaths, vector<TString> folders, bool report = true) {
  vector<TString> fnames;
  for (unsigned ipath = 0; ipath < samplebasepaths.size(); ++ipath){
    for (unsigned ifolder = 0; ifolder < folders.size(); ++ifolder){
      for (unsigned irun = 1; irun <= 10; ++irun) {
        for (unsigned itag = 1; itag < 3; ++itag) {
          TString filename = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i/tag_%.1i_delphes_events.root",irun, itag );
          if (fileexists(filename)) {
            fnames.push_back(filename);
          }
        }
      }
    }
  }
  if (report) {
    cout << " Found sample files are :" << endl;
    for (unsigned i = 0; i< fnames.size(); ++i) {
        cout << fnames[i]<<endl;
    }
    cout << "Total # of sample files: " << fnames.size() <<endl;
  }
  return fnames;
}

// TLorentzVector ReconParticle(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, const float ParticleMass = 80.4, bool erase = false) {
//   TLorentzVector reconW;
//   float recomass = 999999.;
//   int parent1 = -1;
//   int parent2 = -1;
//   for (unsigned i1 = 0; i1 < v1.size(); ++i1 ) {
//     int i2s = 0;
//     if (v1 == v2) i2s =  i1+1;
//     for (unsigned i2 = i2s ; i2 < v2.size(); ++i2 ) {
//       TLorentzVector w = v1.at(i1) + v2.at(i2);
//       float mass = (w).M();
//       if (fabs(mass-ParticleMass) < fabs(recomass-ParticleMass)) {
//         recomass = mass;
//         reconW = w;
//         parent1 = (int)i1;
//         parent2 = (int)i2;
//       }
//     }
//   }
//   // if (parent1 == parent2) cout << "same parent" << endl;
//   if (parent1 == parent2 && v1 == v2) cout << "SAME PARENT !!!!!" <<endl;
//   if (parent1 == -1) cout << "parent1 unset" <<endl;
//   if (parent2 == -1) cout << "parent2 unset" <<endl;
//   if (erase) {
//     v1.erase(v1.begin()+parent1);
//     v2.erase(v2.begin()+parent2);
//   }
//   return reconW;
// }
//
// int FindGenMo(TClonesArray* branchGen, vector<int> list, int targetMo) {
//   int out = -1;
//   for (unsigned igent = 0; igent < list.size(); ++igent) {
//     GenParticle* gentpar = (GenParticle*) branchGen->At(list[igent]);
//     int M1 = gentpar->M1;
//     int M2 = gentpar->M2 >= 0 ? gentpar->M2 : gentpar->M1;
//     for (int imo = M1; imo <= M2 ; ++imo){
//       if (imo == targetMo) out = list[igent];
//     }
//   }
//   return out;
// }
//
// vector<int> FindGenParticle(TClonesArray* GenParticles, string genparticle, int status = 2) {
//   map<string, int> pdgid;
//   pdgid["d"]=1;
//   pdgid["u"]=2;
//   pdgid["s"]=3;
//   pdgid["c"]=4;
//   pdgid["b"]=5;
//   pdgid["t"]=6;
//   pdgid["e"]=11;
//   pdgid["ve"]=12;
//   pdgid["mu"]=13;
//   pdgid["vmu"]=14;
//   pdgid["tau"]=15;
//   pdgid["g"]=21;
//   pdgid["gamma"]=22;
//   pdgid["w"]=24;
//   pdgid["wprime"]=34;
//
//   int idn = pdgid.find(genparticle)->second;
//   // cout << "looing for pdgid: " <<idn <<endl;
//   vector<int> out;
//   for (int igen = 0; igen < GenParticles->GetEntries(); ++igen) {
//     GenParticle* genp = (GenParticle*) GenParticles->At(igen);
//     if (fabs(genp->PID) == idn && fabs(genp->Status) > status * 10 && fabs(genp->Status) < (status + 1) *10 ) out.push_back(igen);
//   }
//   return out;
//
// }
#endif
