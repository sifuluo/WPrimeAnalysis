#include <TROOT.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>

void DrawTH1(TFile* fFL, TFile* fLL, TFile* fBG, string plot, string title, bool doscale = true, bool drawbg = true) {
  const double FLxsection = 1.218;
  const double LLxsection = 1.229;
  const double BGxsection = 127100.;
  const double FLEvt = 81215.;
  const double LLEvt = 81262.;
  const double BGEvt = 900788.;
  const double FLScale = FLxsection / FLEvt; //0.00001499
  const double LLScale = LLxsection / LLEvt; //0.00001512
  const double BGScale = BGxsection / BGEvt; //0.1411

  TCanvas *c1 = new TCanvas(("c"+plot).c_str(), title.c_str(),800,500);
  TH1F* fl = (TH1F*)fFL->Get(plot.c_str());
  fl->SetLineColor(2);
  fl->Scale(FLScale);
  TH1F* ll = (TH1F*)fLL->Get(plot.c_str());
  ll->SetLineColor(2);
  ll->Scale(LLScale);
  TH1F* bg = (TH1F*)fBG->Get(plot.c_str());
  bg->SetLineColor(2);
  bg->Scale(BGScale);
  TLegend* leg = new TLegend(0.8,0.8,1.,1.);
  leg->AddEntry(fl,"FL","l");
  leg->AddEntry(ll,"LL","l");
  leg->AddEntry(bg,"BG","l");
  bg->Draw("");
  fl->Draw("same");
  ll->Draw("same");
  leg->Draw();
  c1->Draw();
  TString savename = "plots/hypothesis/"+plot+".pdf";
  c1->SaveAs(savename);
}
