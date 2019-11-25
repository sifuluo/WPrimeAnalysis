// #include "Analyzer.cc"
#include "Utilities/JESTools.cc"
#include "Utilities/ROOTMini.cc"


void includetest() {

  JESTools *b  = new JESTools();
  cout <<"Successfully included JESTools"<<endl;

  ROOTMini *m = new ROOTMini(b);
  cout <<"Successfully loaded ROOTMini" <<endl;
}
