#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH2.h>

#include <iostream>
#include <cmath>
#include <string>

using namespace std;

void Significance(){
  gStyle->SetOptStat(0);
  TFile* fBG = new TFile("results/hypothesis_BG.root");
  TFile* fLL = new TFile("results/hypothesis_LL.root");
  TFile* fFL = new TFile("results/hypothesis_FL.root");
  TFile* fSignificance = new TFile("plots/hypothesis/Significance.root","RECREATE");
  fSignificance->cd();
  const double FLxsection = 1.218;
  const double LLxsection = 1.229;
  const double BGxsection = 127100.;
  const double FLEvt = 81215.;
  const double LLEvt = 81262.;
  const double BGEvt = 900789.;
  const double FLScale = FLxsection / FLEvt; //0.00001499
  const double LLScale = LLxsection / LLEvt; //0.00001512
  const double BGScale = BGxsection / BGEvt; //0.1411

  TH2F* ExSigFL = new TH2F("ExSigFL","Exclusive Significance of FL; ## N b-jets; ## b-jets",20,-0.5,19.5,20,-0.5,19.5);
  TH2F* ExSigLL = new TH2F("ExSigLL","Exclusive Significance of LL; ## N b-jets; ## b-jets",20,-0.5,19.5,20,-0.5,19.5);
  TH2F* InSigFL = new TH2F("InSigFL","Inclusive Significance of FL; ## N b-jets; ## b-jets",20,-0.5,19.5,20,-0.5,19.5);
  TH2F* InSigLL = new TH2F("InSigLL","Inclusive Significance of LL; ## N b-jets; ## b-jets",20,-0.5,19.5,20,-0.5,19.5);

  TH2F* hBG = (TH2F*) fBG->Get("NBJetsVsBJets");
  TH2F* hFL = (TH2F*) fFL->Get("NBJetsVsBJets");
  TH2F* hLL = (TH2F*) fLL->Get("NBJetsVsBJets");
  hBG->Scale(BGScale);
  hFL->Scale(FLScale);
  hLL->Scale(LLScale);

  int NBinX = hBG->GetXaxis()->GetNbins();
  int NBinY = hBG->GetYaxis()->GetNbins();

  for (int ibinx = 1; ibinx <= NBinX; ++ibinx) {
    for (int ibiny = 1; ibiny <= NBinY; ++ibiny) {
      double ExEvtFL = hFL->GetBinContent(ibinx,ibiny);
      double ExEvtLL = hLL->GetBinContent(ibinx,ibiny);
      double ExEvtBG = hBG->GetBinContent(ibinx,ibiny);
      double InEvtFL = hFL->Integral(ibinx,-1,ibiny,-1);
      double InEvtLL = hLL->Integral(ibinx,-1,ibiny,-1);
      double InEvtBG = hBG->Integral(ibinx,-1,ibiny,-1);

      if (hFL->GetBinError(ibinx,ibiny) / ExEvtFL < 0.3) {
        double exsigFL = ExEvtFL / ( sqrt(ExEvtFL + ExEvtBG) );
        ExSigFL->SetBinContent(ibinx,ibiny,exsigFL);
        double insigFL = InEvtFL / ( sqrt(InEvtFL + InEvtBG) );
        InSigFL->SetBinContent(ibinx,ibiny,insigFL);
      }

      if (hFL->GetBinError(ibinx,ibiny) / ExEvtFL < 0.3) {
        double exsigLL = ExEvtLL / ( sqrt(ExEvtLL + ExEvtBG) );
        ExSigLL->SetBinContent(ibinx,ibiny,exsigLL);
        double insigLL = InEvtLL / ( sqrt(InEvtLL + InEvtBG) );
        InSigLL->SetBinContent(ibinx,ibiny,insigLL);
      }
    }
  }
  fSignificance->Write();
  fSignificance->Save();

}
