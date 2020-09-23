#ifndef SEARCHFILES_CC
#define SEARCHFILES_CC

// #include <TROOT.h>
#include <TString.h>
// #include <TVector2.h>
// #include <TH1.h>
// #include <TH2.h>
// #include <TClonesArray.h>
// #include <TLorentzVector.h>
#include <vector>
#include <iostream>
// #include <fstream>
// #include <cmath>
// #include <map>
// #include <TProfile.h>
// #include <TEfficiency.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>

using namespace std;

bool fileexists (const TString& name) {
  struct stat buffer;
  return ( stat(name, &buffer) == 0);
}

// bool fileexists2 (const TString& name) {
//   ifstream f(name);
//   return f.good();
// }

vector<TString> SearchFiles(int SampleType = 0, int irun = 0, bool report = true) {
  // Parsing Inputs
  vector<TString> samplebasepaths;
  // samplebasepaths.push_back("/fdata/hepx/store/user/aoverton0342/madGraph/ak4/");
  samplebasepaths.push_back("/eos/user/s/siluo/samples/");
  vector<TString> folders;
  if (SampleType == 0) {
    // folders.push_back("TDual_FormerLeptonic/");
    folders.push_back("FL/");
  }
  else if (SampleType == 1) {
    // folders.push_back("TDual_LatterLeptonic/");
    folders.push_back("LL/");
  }
  else {
    // folders.push_back("ttbar/");
    folders.push_back("BG/");
  }

  // Finding files
  vector<TString> fnames;
  for (unsigned ipath = 0; ipath < samplebasepaths.size(); ++ipath){
    for (unsigned ifolder = 0; ifolder < folders.size(); ++ifolder){
      if (irun != 0) {
        unsigned itrun = irun;
        for (unsigned itag = 1; itag < 2; ++itag) {
          // TString filename = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i/tag_%.1i_delphes_events.root",itrun, itag);
          TString filename = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i.root",itrun);
          if (fileexists(filename)) {
            fnames.push_back(filename);
          }
        }
      }
      else {
        unsigned itrun = 1;
        // TString runname = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i",itrun);
        TString runname = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i.root",itrun);
        int skipcount = 0;
        while (fileexists(runname) || skipcount < 4){
          if (!fileexists(runname)) {
            skipcount++;
            itrun++;
            continue;
          }
          skipcount = 0;
          for (unsigned itag = 1; itag < 2; ++itag) {
            // TString filename = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i/tag_%.1i_delphes_events.root",itrun, itag );
            TString filename = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i.root",itrun);
            if (fileexists(filename)) {
              fnames.push_back(filename);
            }
          }
          itrun++;
          // runname = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i",itrun);
          runname = samplebasepaths.at(ipath) + folders.at(ifolder) + TString::Format("Events/run_%.2i.root",itrun);
        }
      }
    }
  }
  if (fnames.empty()) {
    cout << "!!!   No File Found in SampleType: "<< SampleType << ", run: " <<irun <<"   !!!" <<endl;
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

#endif
