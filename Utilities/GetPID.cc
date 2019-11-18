#ifndef GETPID_CC
#define GETPID_CC

#include<cmath>
#include<string>

#include <TROOT.h>
#include <TString.h>

string GetParName(int pid_) {
  int pid = fabs(pid_);
  string par;
  switch (pid) {
    case 1: par = "d";
    break;
    case 2: par = "u";
    break;
    case 3: par = "s";
    break;
    case 4: par = "c";
    break;
    case 5: par = "b";
    break;
    case 6: par = "t";
    break;
    case 11: par = "e";
    break;
    case 12: par = "nue";
    break;
    case 13: par = "mu";
    break;
    case 14: par = "numu";
    break;
    case 15: par = "tau";
    break;
    case 16: par = "nutau";
    break;
    case 21: par = "g";
    break;
    case 22: par = "gamma";
    break;
    case 24: par = "W";
    break;
    case 34: par = "W'";
    break;
    default: par = "??";
    break;
  }
  return par;
}


#endif
