#include "Utilities/Analyzer.cc"
#include "Utilities/JetMatch.cc"

#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TEfficiency.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>

using namespace std;

void MakePFile(int SampleType = 0, int irun = 0, int debug = -2) {
  Analyzer *a = new Analyzer(SampleType, irun, 30);
  TString savepath = "PFile/";
  TString savename = "PFile";
  a->SetOutput(savepath);
  a->DebugMode(debug);

  vector<double> etabins{0., 1.3, 2.5, 3.0, 5.2};
  vector<vector<double> > ptbins{
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130.,150.,180.,220.,260.,300.,350.,400.,500.,1000.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,300.,6000.,0,0,0,0,0},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,6000.,0,0,0,0,0,0},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.,0,0,0,0,0,0,0,0,0}
  };

  vector<vector <TH1F*> > pJet;

}
