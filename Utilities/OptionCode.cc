#ifndef OPTIONCODE_CC
#define OPTIONCODE_CC

#include <cmath>
#include <iostream>

vector<int> ToBinary(int n) {
  vector<int> out;
  out.clear();
  if (n / 2 != 0) {
    out = ToBinary(n / 2);
  }
  out.push_back(n%2);
  // cout << "pushing back " << n%2 <<endl;
  return out; // out[0] is the largest rank, so printing out will be the converted binary
}

void OptionCode(int n) {
  vector<int> out = ToBinary(n);
  cout << "Converting " << n << " to binary is:";
  for (unsigned i = 0; i < out.size(); ++i) {
    cout <<" "<< out[i];
  }
  cout <<endl;
}

#endif
