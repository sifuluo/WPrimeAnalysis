#include <TRandom.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TEfficiency.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>

void TreeTest() {
  TFile hfile("htree.root","RECREATE","Demo ROOT file with histograms & trees");

  // Create some histograms and a profile histogram
  TH1F hpx("hpx","This is the px distribution",100,-4,4);
  TH2F hpxpy("hpxpy","py ps px",40,-4,4,40,-4,4);
  TProfile hprof("hprof","Profile of pz versus px",100,-4,4,0,20);

  // Define some simple structures
  typedef struct {Float_t x,y,z;} POINT;
  typedef struct {
     Int_t ntrack,nseg,nvertex;
     UInt_t flag;
     Float_t temperature;
  } EVENTN;
  POINT point;
  EVENTN eventn;
  TLorentzVector lv;
  vector<double>* vec;

  // Create a ROOT Tree
  TTree tree("T","An example of ROOT tree with a few branches");
  tree.Branch("point",&point,"x:y:z");
  tree.Branch("eventn",&eventn,"ntrack/I:nseg:nvertex:flag/i:temperature/F");
  tree.Branch("hpx","TH1F",&hpx,128000,0);
  tree.Branch("lv","TLorentzVector",&lv);
  tree.Branch("vec",&vec);

  Float_t px,py,pz;

  // Here we start a loop on 1000 events
  for ( Int_t i=0; i<1000; i++) {
     gRandom->Rannor(px,py);
     pz = px*px + py*py;
     const auto random = gRandom->Rndm(1);
     lv = TLorentzVector(px,py,pz,0);
     vec->clear();
     vec->push_back(px);
     vec->push_back(py);
     vec->push_back(pz);
     // Fill histograms
     hpx.Fill(px);
     hpxpy.Fill(px,py,1);
     hprof.Fill(px,pz,1);

     // Fill structures
     point.x = 10*(random-1);
     point.y = 5*random;
     point.z = 20*random;

     eventn.ntrack  = Int_t(100*random);
     eventn.nseg    = Int_t(2*eventn.ntrack);
     eventn.nvertex = 1;
     eventn.flag    = Int_t(random+0.5);
     eventn.temperature = 20+random;
     tree.Fill();
     // Fill the tree. For each event, save the 2 structures and 3 objects
     // In this simple example, the objects hpx, hprof and hpxpy are slightly
     // different from event to event. We expect a big compression factor!

  }
  // End of the loop

  tree.Print();

  // Save all objects in this file
  hfile.Write();

  // Close the file. Note that this is automatically done when you leave
  // the application upon file destruction.
  hfile.Close();


}
