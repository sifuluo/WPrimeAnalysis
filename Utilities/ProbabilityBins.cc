#ifndef PROBABILITYBINS_CC
#define PROBABILITYBINS_CC

#include <vector>

void ProbabilityBins(vector< vector<double> > &ptbins, vector<int> &npt, vector<double> &etabins) {
  double etaarr[] = {0.,1.3,2.5,3.0,5.2};
  etabins.clear();
  for (unsigned i = 0; i < (sizeof(etaarr)/sizeof(etaarr[0])); ++i) {
    etabins.push_back(etaarr[i]);
  }

  double pttemp[4][23] = {
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,130.,150.,180.,220.,260.,300.,350.,400.,500.,1000.,6000.},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,300.,6000.,0,0,0,0,0},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,180.,220.,260.,6000.,0,0,0,0,0,0},
    {30.,32.,34.,37.,40.,45.,50.,57.,65.,75.,90.,110.,150.,6000.,0,0,0,0,0,0,0,0,0}
  };
  ptbins.clear();
  npt.clear();
  for (unsigned ieta = 0 ; ieta < etabins.size() - 1; ++ieta){
    npt.push_back(sizeof(pttemp[ieta])/sizeof(pttemp[ieta][0]));
    vector<double> ptbineta;
    for (int ipt = 0; ipt < npt[ieta] ; ++ipt) {
      if (pttemp[ieta][ipt] == 0. ) {npt[ieta] = ipt;break;};
      ptbineta.push_back(pttemp[ieta][ipt]);
    }
    ptbins.push_back(ptbineta);
  }
}
#endif
