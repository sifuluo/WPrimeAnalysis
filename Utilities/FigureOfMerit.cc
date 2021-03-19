#ifndef FIGUREOFMERIT_CC
#define FIGUREOFMERIT_CC

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TMinuit.h"

#include <utility>
#include <iostream>
#include <vector>
#include <cmath>

// #include "JetMETAnalysis/JetUtilities/interface/TProfileMDF.h"

using std::cout;
using std::endl;
using std::vector;
using std::pair;

class FigureOfMerit{

public:

  FigureOfMerit();
  ~FigureOfMerit();

  //Use this when working with arbitrary dimensions and using cuts
  // static double usingShapeFromTemplatesMD(TH1* signal, TH1* background, double relativeBinErrorLimit = 3);

  // This returns false if the Histogram has an empty bin
  // static bool checkBackgroundHistogramMD(TProfileMDF* back, double relativeBinErrorLimit = 3);

  // Use this when finding limits on cross sections, as with the WH analysis
  static double usingShapeFromTemplates(TH1* signal, TH1* background, double relativeBinErrorLimit = 3);
  static double usingShapeFromTemplates2D(TH2* signal, TH2* background, double relativeBinErrorLimit = 3);

  // This computes the sqrt[sum_bin S^2/(S+B+deltaS^2+deltaB^2)]
  static double usingChi2(TH1* signal, TH1* background, double relativeBinErrorLimit = 3);

  // This computes the sqrt[sum_bin S^2/B]
  static double usingS2OverB(TH1* signal, TH1* background);

  // This returns false if the Histogram has an empty bin
  static bool checkBackgroundHistogram(TH1* back, double relativeBinErrorLimit = 3);
  static bool checkBackgroundHistogram2D(TH2* back, double relativeBinErrorLimit = 3);

private:

  //These are used in usingShapeFromTemplates for the fitting procedure
  static TH1* signal;
  static TH1* background;
  static TH2* signal2D;
  static TH2* background2D;
  static void minimizationFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void minimizationFunction2D(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  // static void minimizationFunctionMD(Int_t &npar, Double_t *gin, Double_t &f,
  // Double_t *par, Int_t iflag);

};

//Declare the static members
TH1* FigureOfMerit::signal;
TH1* FigureOfMerit::background;
TH2* FigureOfMerit::signal2D;
TH2* FigureOfMerit::background2D;

// ----------------------------------------------------------------------------
FigureOfMerit::FigureOfMerit(){

}// C'tor


// ----------------------------------------------------------------------------
FigureOfMerit::~FigureOfMerit(){

}// D'tor

// ----------------------------------------------------------------------------
// This checks a histogram for empty bins
// Returns true if no empty bins
bool FigureOfMerit::checkBackgroundHistogram(TH1* back, double relativeBinErrorLimit){

  // Report an error if there is an something in the
  // underflow or overflow bins
  if (back->GetBinContent(0)       > 0 ||
  back->GetBinContent(back->GetNbinsX()+1) > 0 ){

    cout<<" WARNING FigureOfMerit::checkBackgroundHistogram background"
    <<" histogram have underflow/overflow bins."<<endl;
  }

  // Loop over all bins
  vector<pair<int,int> > unfilled_bin_pairs;
  int unfilled_start = 0;
  int unfilled_end = 0;
  int last_unfilled_bin = 0;
  for (int ibin = 1; ibin <= back->GetNbinsX(); ++ibin){

    // If there is no background in this bin...
    if (back->GetBinContent(ibin) == 0){

      if(unfilled_start == 0) {
        last_unfilled_bin = ibin;
        unfilled_start = ibin;
      }
      if (ibin-last_unfilled_bin > 1 || ibin == back->GetNbinsX()) {
        if(ibin == back->GetNbinsX())
        unfilled_end = ibin;
        else
        unfilled_end = last_unfilled_bin;
        last_unfilled_bin = ibin;
        unfilled_bin_pairs.push_back(std::make_pair(unfilled_start,unfilled_end));
        unfilled_start = ibin;
      }
      else if(ibin-last_unfilled_bin == 1) {
        last_unfilled_bin = ibin;
      }
      else
      continue;
      //std::cerr << "ERROR FigureOfMerit::checkBackgroundHistogram bin "<<ibin<<" is not filled."<<endl;
      //return false;

    }
    // If the background contribution is smaller than three time its error return a warning.
    else if (back->GetBinContent(ibin)  < relativeBinErrorLimit * back->GetBinError(ibin) ){

      std::cerr << "WARNING FigureOfMerit::checkBackgroundHistogram "
      <<"  background error in bin "<<ibin<<" is too large!\n";
      return false;
    }

  }//for
  for (unsigned int i=0; i<unfilled_bin_pairs.size(); i++) {
    std::cerr << "ERROR FigureOfMerit::checkBackgroundHistogram bins [" << unfilled_bin_pairs[i].first
    << "-" << unfilled_bin_pairs[i].second << "] are not filled." << endl;
  }

  return true;

}// checkBackgroundHistogram;

bool FigureOfMerit::checkBackgroundHistogram2D(TH2* back, double relativeBinErrorLimit){

  // Report an error if there is an something in the
  // underflow or overflow bins
  vector<vector<pair<int,int> > > unfilled_bin_pairs_2D;
  for (int iy = 1; iy < back->GetNbinsY()+1; ++iy) {
    if (back->GetBinContent(0,iy) > 0 || back->GetBinContent(back->GetNbinsX()+1,iy) > 0 ){

      cout<<" WARNING FigureOfMerit::checkBackgroundHistogram background"
      <<" histogram have underflow/overflow bins."<<endl;
    }

    // Loop over all bins
    vector<pair<int,int> > unfilled_bin_pairs;
    int unfilled_start = 0;
    int unfilled_end = 0;
    int last_unfilled_bin = 0;
    for (int ibin = 1; ibin <= back->GetNbinsX(); ++ibin){
      // If there is no background in this bin...
      if (back->GetBinContent(ibin,iy) == 0){

        if(unfilled_start == 0) {
          last_unfilled_bin = ibin;
          unfilled_start = ibin;
        }
        if (ibin-last_unfilled_bin > 1 || ibin == back->GetNbinsX()) {
          if(ibin == back->GetNbinsX())
          unfilled_end = ibin;
          else
          unfilled_end = last_unfilled_bin;
          last_unfilled_bin = ibin;
          unfilled_bin_pairs.push_back(std::make_pair(unfilled_start,unfilled_end));
          unfilled_start = ibin;
        }
        else if(ibin-last_unfilled_bin == 1) {
          last_unfilled_bin = ibin;
        }
        else
        continue;
        //std::cerr << "ERROR FigureOfMerit::checkBackgroundHistogram bin "<<ibin<<" is not filled."<<endl;
        //return false;
      }
      // If the background contribution is smaller than three time its error return a warning.
      else if (back->GetBinContent(ibin, iy)  < relativeBinErrorLimit * back->GetBinError(ibin, iy) ){
        std::cerr << "WARNING FigureOfMerit::checkBackgroundHistogram ";
        std::cerr <<"  background error in bin "<<ibin<<" , " << iy << " is too large!\n";
        std::cerr <<"  BinContent = " << back->GetBinContent(ibin, iy)<<"; Bin Error = " << back->GetBinError(ibin, iy) <<endl;
        return false;
      }
    }
    unfilled_bin_pairs_2D.push_back(unfilled_bin_pairs);
  }
  for (unsigned iy = 0; iy < unfilled_bin_pairs_2D.size(); ++iy) {
    for (unsigned int i=0; i<unfilled_bin_pairs_2D[iy].size(); i++) {
      std::cerr << "ERROR FigureOfMerit::checkBackgroundHistogram BinY = "<< iy+1 << " , bins [" << unfilled_bin_pairs_2D[iy][i].first
      << "-" << unfilled_bin_pairs_2D[iy][i].second << "] are not filled." << endl;
    }
  }

  return true;

}// checkBackgroundHistogram;

// ----------------------------------------------------------------------------
// This is the one to use when estimating cross section as one would do calculate
// limits on the cross section.
// Note that if trying to find the cross section of the signal (sigma_sig) the
// figure of merit (fom) returned here is sigma_sig=1/fom. This is, the fom is
// the inverse of the error on the cross section.
double FigureOfMerit::usingShapeFromTemplates(TH1* _signal, TH1* _background, double relativeBinErrorLimit){

  // The returning figure of merit
  double fom = 0;

  //Try to run only of the background histo does not have empty bins
  if (checkBackgroundHistogram(_background, relativeBinErrorLimit)){

    // Assign both histograms to static functions
    signal     = _signal;
    background = _background;

    // Create the minimization fitter with one parameter.
    // Do not create minuit as an object, always create in the heap,
    // otherwise it seems that running twice gives different results!!
    TMinuit  * minuit = new TMinuit(1);

    // Set the minimization function
    minuit->SetFCN(minimizationFunction);

    Double_t arglist[10];
    Int_t ierflg = 0;

    // Make Minuit run in a quiet (i.e. non-verbose) mode
    minuit->SetPrintLevel(-1);
    arglist[0] = 0;
    gMinuit->mnexcm("SET NOWarnings",&arglist[0],1,ierflg);
    arglist[0] = -1;
    gMinuit->mnexcm("SET PRint",&arglist[0],1,ierflg);


    // Set the error as half a point around the NLL
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist ,1,ierflg);

    // Set starting values and step sizes for parameters
    //            (par# , name  start val , step  , min val,  max val, ierflg);
    //minuit->mnparm(    0, "a1",   0.50    , 0.1   ,  0     ,   200  ,
    minuit->mnparm(    0, "a1",   0.50    , 0.1   ,  0     ,   2000  ,
    ierflg);

    // Do the minimization
    arglist[0] = 0.5; // 0.5 for 1 sigma. 2 for 2 sigmas
    minuit->mnexcm("SET NOG",arglist,0,ierflg); //
    minuit->mnexcm("MIGRAD", arglist ,2,ierflg); // -
    arglist[0] = 50;
    arglist[1] = 0.01;
    minuit->mnexcm("MINOS",arglist,2,ierflg);

    // Retrieve the central value (central) and errors (up,down)
    double central, err_down, err_up, junk, junk1, junk2;
    minuit->GetParameter(0, central, junk);
    minuit->mnerrs(0, err_up, err_down, junk1, junk2);

    // Report ?
    // cout<<" central Sig ="<<central<<" +"<<err_up<<"-"<<err_down<<" junk ="<<junk<<endl;
    // cout<<" Number of calls = "<<minuit->fNfcn<<endl;
    // Clean up
    delete minuit;

    // Do a basic check
    if (err_up ==0){
      cout<<" ERROR  FigureOfMerit::usingShapeFromTemplates err_up= "<<err_up<<endl
      <<"   This can sometimes happen when the fom is < 1/max_val"<<endl
      <<"   To solve it just increase the value of max_val! " <<endl;
      return 0;
    }

    // The figure of merit is one over the error
    fom = 1/err_up;

  }// background does not have empty bins

  return fom;

}// usingShapeFromTemplates

double FigureOfMerit::usingShapeFromTemplates2D(TH2* _signal, TH2* _background, double relativeBinErrorLimit){

  // The returning figure of merit
  double fom = 0;
  //Try to run only of the background histo does not have empty bins
  if (checkBackgroundHistogram2D(_background, relativeBinErrorLimit)){

    // Assign both histograms to static functions
    signal2D     = _signal;
    background2D = _background;

    // Create the minimization fitter with one parameter.
    // Do not create minuit as an object, always create in the heap,
    // otherwise it seems that running twice gives different results!!
    TMinuit  * minuit = new TMinuit(1);

    // Set the minimization function
    minuit->SetFCN(minimizationFunction2D);

    Double_t arglist[10];
    Int_t ierflg = 0;

    // Make Minuit run in a quiet (i.e. non-verbose) mode
    minuit->SetPrintLevel(-1);
    arglist[0] = 0;
    gMinuit->mnexcm("SET NOWarnings",&arglist[0],1,ierflg);
    arglist[0] = -1;
    gMinuit->mnexcm("SET PRint",&arglist[0],1,ierflg);


    // Set the error as half a point around the NLL
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist ,1,ierflg);

    // Set starting values and step sizes for parameters
    //            (par# , name  start val , step  , min val,  max val, ierflg);
    //minuit->mnparm(    0, "a1",   0.50    , 0.1   ,  0     ,   200  ,
    minuit->mnparm(    0, "a1",   0.50    , 0.1   ,  0     ,   2000  ,
    ierflg);

    // Do the minimization
    arglist[0] = 0.5; // 0.5 for 1 sigma. 2 for 2 sigmas
    minuit->mnexcm("SET NOG",arglist,0,ierflg); //
    minuit->mnexcm("MIGRAD", arglist ,2,ierflg); // -
    arglist[0] = 50;
    arglist[1] = 0.01;
    minuit->mnexcm("MINOS",arglist,2,ierflg);

    // Retrieve the central value (central) and errors (up,down)
    double central, err_down, err_up, junk, junk1, junk2;
    minuit->GetParameter(0, central, junk);
    minuit->mnerrs(0, err_up, err_down, junk1, junk2);

    // Report ?
    // cout<<" central Sig ="<<central<<" +"<<err_up<<"-"<<err_down<<" junk ="<<junk<<endl;
    // cout<<" Number of calls = "<<minuit->fNfcn<<endl;
    // Clean up
    delete minuit;

    // Do a basic check
    if (err_up ==0){
      cout<<" ERROR  FigureOfMerit::usingShapeFromTemplates err_up= "<<err_up<<endl
      <<"   This can sometimes happen when the fom is < 1/max_val"<<endl
      <<"   To solve it just increase the value of max_val! " <<endl;
      return 0;
    }

    // The figure of merit is one over the error
    fom = 1/err_up;

  }// background does not have empty bins

  return fom;

}// usingShapeFromTemplates2D
// ----------------------------------------------------------------------------
// The minimization function used in ::usingShapeFromTemplates
void FigureOfMerit::minimizationFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  //calculate chisquare
  Double_t sum_nll = 0;
  for (int ibin = 1 ; ibin <= signal->GetNbinsX() ; ibin++) {

    // skip bins in which both signal and background are zero
    if (background->GetBinContent(ibin) == 0 ){
      if ( signal->GetBinContent(ibin) ==0)
      continue;
      else{
        cout<<" ERROR FigureOfMerit::minimizationFunction "
        <<" at bin="<<ibin<<" background is zero but signal is NOT!"
        <<" This can give infinite sensitivity!"<<endl;
      }
    }//if

    double bin_num  = par[0]*signal->GetBinContent(ibin) * par[0]*signal->GetBinContent(ibin);
    double bin_den  = par[0]*signal->GetBinContent(ibin) + background->GetBinContent(ibin)
    + par[0]*signal->GetBinError(ibin) * par[0]*signal->GetBinError(ibin)
    + background->GetBinError(ibin) * background->GetBinError(ibin);

    sum_nll +=  bin_num/bin_den;

  }//for

  //cout<<" par[0]="<<par[0]<<" sum_nll="<<sum_nll<<endl;

  f = sqrt(sum_nll);

}// minimizationFunction

void FigureOfMerit::minimizationFunction2D(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  //calculate chisquare
  Double_t sum_nll = 0;
  for (int iy = 1; iy <= signal2D->GetNbinsY(); ++iy) {
    for (int ibin = 1 ; ibin <= signal2D->GetNbinsX() ; ibin++) {

      // skip bins in which both signal and background are zero
      if (background2D->GetBinContent(ibin,iy) == 0 ){
        continue;
        if ( signal2D->GetBinContent(ibin,iy) ==0)
        continue;
        else{
          cout<<" ERROR FigureOfMerit::minimizationFunction "
          <<" at bin="<<ibin<<" background is zero but signal is NOT!"
          <<" This can give infinite sensitivity!"<<endl;
        }
      }//if

      double bin_num  = par[0]*signal2D->GetBinContent(ibin, iy) * par[0]*signal2D->GetBinContent(ibin, iy);
      double bin_den  = par[0]*signal2D->GetBinContent(ibin, iy) + background2D->GetBinContent(ibin, iy)
      + par[0]*signal2D->GetBinError(ibin, iy) * par[0]*signal2D->GetBinError(ibin, iy)
      + background2D->GetBinError(ibin, iy) * background2D->GetBinError(ibin, iy);

      sum_nll +=  bin_num/bin_den;

    }//for
  }


  //cout<<" par[0]="<<par[0]<<" sum_nll="<<sum_nll<<endl;

  f = sqrt(sum_nll);

}// minimizationFunction


// ----------------------------------------------------------------------------
// This method computes the sensitivity as
//     S = sqrt(sum_bin S_bin^2/ (S_bin+B_bin+DeltaS_bin^2+DeltaB_bin^2))
// The sensitivity is set to zero if any given bin does not have a background
// contribution, or  if that contribution is smaller than 3 times the error
// on that contribution
double FigureOfMerit::usingChi2(TH1* signal, TH1* background, double relativeBinErrorLimit){

  // Check that the number of bins is the same
  if (signal->GetNbinsX() != background->GetNbinsX()){
    cout<<"   FigureOfMerit::usingChi2 was given two histograms of different sizes!"<<endl;
    return 0;
  }

  //Check that this histo passes the basic cuts
  if (checkBackgroundHistogram(background, relativeBinErrorLimit)){

    // Get the number of bins
    unsigned int nBins = signal->GetNbinsX();

    // Loop over bins
    double total = 0;
    for (unsigned ibin = 1; ibin <= nBins; ++ibin){

      double dataStat2 = signal->GetBinContent(ibin) + background->GetBinContent(ibin) ;
      double mcStat2   = signal->GetBinError(ibin)     * signal->GetBinError(ibin)
      + background->GetBinError(ibin) * background->GetBinError(ibin);

      // Compute this bin's sensitivity
      double signal2         = signal->GetBinContent(ibin)*signal->GetBinContent(ibin);
      double bin_sensitivity = 0;
      if((dataStat2+mcStat2)!=0)
      bin_sensitivity = signal2 / ( dataStat2 + mcStat2 );
      else
      bin_sensitivity = 0;

      // Add it to the total
      total += bin_sensitivity;

    }//for bins

    //return the square of the total
    return sqrt(total);

  }

  // if we are still here return zero
  return 0;

}//using Chi2


// ----------------------------------------------------------------------------
// This method computes the sensitivity as S = sqrt(sum_bin S_bin^2/B_bin)
// The sensitivity is set to zero if any given bin does not have a background
// contribution, or  if that contribution is smaller than 3 times the error
// on that contribution
double FigureOfMerit::usingS2OverB(TH1* signal, TH1* background){

  // Check that the number of bins is the same
  if (signal->GetNbinsX() != background->GetNbinsX()){
    cout<<"   FigureOfMerit::usingS2OverB was given two histograms of different sizes!"<<endl;
    return 0;
  }

  //Check that this histo passes the basic cuts
  if (checkBackgroundHistogram(background)){

    // Get the number of bins
    unsigned int nBins = signal->GetNbinsX();

    // Loop over bins
    double total = 0;
    for (unsigned ibin = 1; ibin <= nBins; ++ibin){

      // Compute this bin's sensitivity
      double signal2         = signal->GetBinContent(ibin)*signal->GetBinContent(ibin);
      double bin_sensitivity = signal2 / background->GetBinContent(ibin) ;

      // Add it to the total
      total += bin_sensitivity;

    }//for bins

    //return the square of the total
    return sqrt(total);

  }

  // if we are still here return zero
  return 0;

}//using S2/B

/*
// ----------------------------------------------------------------------------
// This checks a histogram for empty bins
// Returns true if no empty bins
bool FigureOfMerit::checkBackgroundHistogramMD(TProfileMDF* back, double relativeBinErrorLimit) {

// Loop over all bins
for (int ibin = 0; ibin<back->GetNbins(); ++ibin){
// Report an error if there is an something in the
// underflow or overflow bins
if(back->IsBinUnderflow(ibin)||back->IsBinOverflow(ibin)) {
if(back->GetBinContent(ibin) > 0) {
cout<<" WARNING FigureOfMerit::checkBackgroundHistogram background"
<<" histogram have underflow/overflow bins."<<endl;
}
}//is overflow/underflow
if(!back->IsBinUnderflow(ibin) && !back->IsBinOverflow(ibin)) {
// If there is no background in this bin...
if (back->GetBinContent(ibin) == 0){

std::cerr << "ERROR FigureOfMerit::checkBackgroundHistogram bin "<<ibin<<" is not filled."<<endl;
//return false;

}
// If the background contribution is smaller than three time its error return a warning.
else if (back->GetBinContent(ibin)  < relativeBinErrorLimit * back->GetBinError(ibin) ){

std::cerr << "WARNING FigureOfMerit::checkBackgroundHistogram "
<<"  background error in bin "<<ibin<<" is too large!\n";
return false;
}
}//not overflow/underflow
}//for

return true;

}//checkBackgroundHistogramMD

// ----------------------------------------------------------------------------
// The minimization function used in ::usingShapeFromTemplates
void FigureOfMerit::minimizationFunctionMD(Int_t &npar, Double_t *gin, Double_t &f,
Double_t *par, Int_t iflag){

//calculate chisquare
Double_t sum_nll = 0;
for (int ibin = 0 ; ibin <= ((TProfileMDF*)signal)->GetNbins() ; ibin++) {
if(!((TProfileMDF*)background)->IsBinUnderflow(ibin) && !((TProfileMDF*)background)->IsBinOverflow(ibin)) {
// skip bins in which both signal and background are zero
if (((TProfileMDF*)background)->GetBinContent(ibin) == 0 ){
if (((TProfileMDF*)signal)->GetBinContent(ibin) == 0 )
continue;
else{
cout<<" ERROR FigureOfMerit::minimizationFunction "
<<" at bin="<<ibin<<" background is zero but signal is NOT!"
<<" This can give infinite sensitivity!"<<endl;
}
}//if

double bin_num  = par[0]*((TProfileMDF*)signal)->GetBinContent(ibin) * par[0]*((TProfileMDF*)signal)->GetBinContent(ibin);
double bin_den  = par[0]*((TProfileMDF*)signal)->GetBinContent(ibin) + ((TProfileMDF*)background)->GetBinContent(ibin)
+ par[0]*((TProfileMDF*)signal)->GetBinError(ibin) * par[0]*((TProfileMDF*)signal)->GetBinError(ibin)
+ ((TProfileMDF*)background)->GetBinError(ibin) * ((TProfileMDF*)background)->GetBinError(ibin);

sum_nll +=  bin_num/bin_den;
}//not overflow/underflow
}//for

//cout<<" par[0]="<<par[0]<<" sum_nll="<<sum_nll<<endl;

f = sqrt(sum_nll);

}// minimizationFunction

// ----------------------------------------------------------------------------
// This is the one to use when estimating cross section as one would do calculate
// limits on the cross section.
// Note that if trying to find the cross section of the signal (sigma_sig) the
// figure of merit (fom) returned here is sigma_sig=1/fom. This is, the fom is
// the inverse of the error on the cross section.
double FigureOfMerit::usingShapeFromTemplatesMD(TH1* _signal, TH1* _background, double relativeBinErrorLimit) {
// The returning figure of merit
double fom = 0;

//Try to run only of the background histo does not have empty bins
if (checkBackgroundHistogram((TProfileMDF*)_background, relativeBinErrorLimit)){

// Assign both histograms to static functions
signal     = _signal;
background = _background;

// Create the minimization fitter with one parameter.
// Do not create minuit as an object, always create in the heap,
// otherwise it seems that running twice gives different results!!
TMinuit  * minuit = new TMinuit(1);

// Set the minimization function
minuit->SetFCN(minimizationFunction);

Double_t arglist[10];
Int_t ierflg = 0;

// Make Minuit run in a quiet (i.e. non-verbose) mode
minuit->SetPrintLevel(-1);
arglist[0] = 0;
gMinuit->mnexcm("SET NOWarnings",&arglist[0],1,ierflg);
arglist[0] = -1;
gMinuit->mnexcm("SET PRint",&arglist[0],1,ierflg);


// Set the error as half a point around the NLL
arglist[0] = 1;
minuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
//            (par# , name  start val , step  , min val,  max val, ierflg);
minuit->mnparm(    0, "a1",   0.50    , 0.1   ,  0     ,   200  ,
ierflg);

// Do the minimization
arglist[0] = 0.5; // 0.5 for 1 sigma. 2 for 2 sigmas
minuit->mnexcm("SET NOG",arglist,0,ierflg); //
minuit->mnexcm("MIGRAD", arglist ,2,ierflg); // -
arglist[0] = 50;
arglist[1] = 0.01;
minuit->mnexcm("MINOS",arglist,2,ierflg);

// Retrieve the central value (central) and errors (up,down)
double central, err_down, err_up, junk, junk1, junk2;
minuit->GetParameter(0, central, junk);
minuit->mnerrs(0, err_up, err_down, junk1, junk2);

// Report ?
// cout<<" central Sig ="<<central<<" +"<<err_up<<"-"<<err_down<<" junk ="<<junk<<endl;
// cout<<" Number of calls = "<<minuit->fNfcn<<endl;
// Clean up
delete minuit;

// Do a basic check
if (err_up ==0){
cout<<" ERROR  FigureOfMerit::usingShapeFromTemplates err_up= "<<err_up<<endl
<<"   This can sometimes happen when the fom is < 1/max_val"<<endl
<<"   To solve it just increase the value of max_val! " <<endl;
return 0;
}

// The figure of merit is one over the error
fom = 1/err_up;

}// background does not have empty bins

return fom;

}//usingShapeFromTemplatesMD
*/

#endif
