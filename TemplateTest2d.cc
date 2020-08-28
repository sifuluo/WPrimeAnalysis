#include "Utilities/Analyzer.cc"
#include "Utilities/JESTools.cc"
#include "Utilities/GenTools.cc"
#include "Utilities/ROOTMini.cc"

#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TEfficiency.h>
#include <TString.h>

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>

#include "AnalysisModules/TemplateTest2d.cc"

using namespace std;

void TemplateTest2d(int SampleType = 0, int irun = 1, int OptionCode = 0, int debug = 0) {
  int DoReco = OptionCode%10;
  int LeadingAsWPB = OptionCode/10%10;

  cout << "start" <<endl;
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "results/";
  TString savename;
  if (DoReco == 0) savename = "GenTemplate2d";
  if (DoReco == 1) savename = "RecoTemplate2d";
  if (LeadingAsWPB == 1) savename += "_WPB";
  cout << "Running as " << savename <<endl;

  a->SetupROOTMini();
  a->SetOutput(savepath,savename);
  a->DebugMode(debug);

  TemplateTest2d_init(a);

  for (Int_t entry = a->GetStartEntry(); entry < a->GetEndEntry(); ++entry) {
    a->ReadEvent(entry);
    TemplateTest2d_loop(a, DoReco, LeadingAsWPB);
  }
  a->SaveOutput();

}
