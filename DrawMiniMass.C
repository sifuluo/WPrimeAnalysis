{
  TFile* f = new TFile("results/MinimizerTest.root");
  TH1F* GenHadWMass = (TH1F*) f->Get("GenHadWMass");
  TH1F* GenLepWMass = (TH1F*) f->Get("GenLepWMass");
  TH1F* GenHadTMass = (TH1F*) f->Get("GenHadTMass");
  TH1F* GenLepTMass = (TH1F*) f->Get("GenLepTMass");
  TH1F* GenPDis = (TH1F*) f->Get("GenPDis");

  TH1F* SolLepTMass = (TH1F*) f->Get("SolLepTMass");
  TH1F* SolPDis = (TH1F*) f->Get("SolPDis");

  TH1F* MinHadWMass = (TH1F*) f->Get("MinHadWMass");
  TH1F* MinLepWMass = (TH1F*) f->Get("MinLepWMass");
  TH1F* MinHadTMass = (TH1F*) f->Get("MinHadTMass");
  TH1F* MinLepTMass = (TH1F*) f->Get("MinLepTMass");
  TH1F* MinPDis = (TH1F*) f->Get("MinPDis");

  TF1* w = new TF1("WBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,200.0);
  TF1* t = new TF1("TBW","[0]*TMath::BreitWigner(x,[1],[2])",0.0,300.0);
  TCanvas *c1 = new TCanvas();
  t->SetParameters(100.,172.7,6.7);//1.32
  w->SetParameters(100,80.385,2.085);

  TH1F* hw = GenHadWMass;
  hw->Draw();
  w->SetParameter(0,100);
  w->SetParameter(0,100*hw->GetMaximum()/w->Eval(80.385));
  w->Draw("same");
  c1->SaveAs("plots/MiniMass/GenHadWMass.png");

  hw = GenLepWMass;
  hw->Draw();
  w->SetParameter(0,100);
  w->SetParameter(0,100*hw->GetMaximum()/w->Eval(80.385));
  w->Draw("same");
  c1->SaveAs("plots/MiniMass/GenLepWMass.png");

  hw = MinHadWMass;
  hw->Draw();
  w->SetParameter(0,100);
  w->SetParameter(0,100*hw->GetMaximum()/w->Eval(80.385));
  w->Draw("same");
  c1->SaveAs("plots/MiniMass/MinHadWMass.png");

  hw = MinLepWMass;
  hw->Draw();
  w->SetParameter(0,100);
  w->SetParameter(0,100*hw->GetMaximum()/w->Eval(80.385));
  w->Draw("same");
  c1->SaveAs("plots/MiniMass/MinLepWMass.png");

  TH1F* ht = GenHadTMass;
  ht->Draw();
  t->SetParameter(0,100);
  t->SetParameter(0,100*ht->GetMaximum()/t->Eval(172.7));
  t->Draw("same");
  c1->SaveAs("plots/MiniMass/GenHadTMass.png");

  ht = GenLepTMass;
  ht->Draw();
  t->SetParameter(0,100);
  t->SetParameter(0,100*ht->GetMaximum()/t->Eval(172.7));
  t->Draw("same");
  c1->SaveAs("plots/MiniMass/GenLepTMass.png");

  ht = MinHadTMass;
  ht->Draw();
  t->SetParameter(0,100);
  t->SetParameter(0,100*ht->GetMaximum()/t->Eval(172.7));
  t->Draw("same");
  c1->SaveAs("plots/MiniMass/MinHadTMass.png");

  ht = MinLepTMass;
  ht->Draw();
  t->SetParameter(0,100);
  t->SetParameter(0,100*ht->GetMaximum()/t->Eval(172.7));
  t->Draw("same");
  c1->SaveAs("plots/MiniMass/MinLepTMass.png");

  ht = SolLepTMass;
  ht->Draw();
  t->SetParameter(0,100);
  t->SetParameter(0,100*ht->GetMaximum()/t->Eval(172.7));
  t->Draw("same");
  c1->SaveAs("plots/MiniMass/SolLepTMass.png");

  GenPDis->Draw();
  c1->SaveAs("plots/MiniMass/GenPDis.png");

  SolPDis->Draw();
  c1->SaveAs("plots/MiniMass/SolPDis.png");

  MinPDis->Draw();
  c1->SaveAs("plots/MiniMass/MinPDis.png");


}
