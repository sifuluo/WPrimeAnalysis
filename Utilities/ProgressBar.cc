#ifndef PROGRESSBAR_CC
#define PROGRESSBAR_CC

#include <cmath>
#include <iostream>
#include <string>
#include <TString.h>

void PrintProgress (int entry, int entries, int dividen = 100) {
  entries--;
  if (entry % dividen ==0 || entry == entries ) {
    cout << Form("\rProgress: %i / %i ", entry, entries) <<flush;
  }
  if (entry == entries){
    cout <<endl << "Done."<<endl;
  }
}

#endif
