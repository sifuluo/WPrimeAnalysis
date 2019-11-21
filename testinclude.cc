#include "Utilities/EtaPtBins.cc"

#include <vector>
#include <iostream>
#include <cmath>

void testinclude(int a, int c){
  EtaPtBins *b = new EtaPtBins();
  double aa = double(a);
  double cc = double(c);
  int ieta = b->iBin(aa,cc).first;
  int ipt = b->iBin(aa,cc).second;
  cout << "Begins" <<endl;
  cout << ieta<<"  " <<ipt <<endl;
  cout << "ends" <<endl;
}
