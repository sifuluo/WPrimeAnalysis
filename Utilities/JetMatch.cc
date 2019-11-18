#ifndef JETMATCH_CC
#define JETMATCH_CC

#include "TLorentzVector.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <map>

double JetPtThreshold = 30;
double AlgodR = 0.2;

pair<int,int> FindMatrixMin(vector< vector<double> > &matrix, vector<int> &list1, vector<int> &list2, double &MatrixMin) {
  int ml1(0), ml2(0);
  MatrixMin = -1;
  for (unsigned il1 = 0; il1 < list1.size(); ++il1) {
    if (list1[il1] == 0) continue;
    for (unsigned il2 = 0; il2 < list2.size(); ++il2) {
      if (list2[il2] == 0) continue;
      if (matrix[il1][il2] < MatrixMin || MatrixMin == -1) {
        MatrixMin = matrix[il1][il2];
        ml1 = il1;
        ml2 = il2;
      }
    }
  }
  return pair<int,int>(ml1,ml2);
}

map<int,int> JetMatch(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, double &maxDeltaR, bool cut, bool FilterSoftGen, bool FilterSoft) {
  bool testmodule = false;
  vector<int> l1, l2;
  for (unsigned il = 0; il < v1.size(); ++il) {
    if ((v1[il].Pt() < JetPtThreshold && FilterSoftGen)||v1[il].Pt() == 0) l1.push_back(0);
    else l1.push_back(1);
  }
  for (unsigned il = 0; il < v2.size(); ++il) {
    if ((v2[il].Pt() < JetPtThreshold && FilterSoft)||v2[il].Pt() == 0) l2.push_back(0);
    else l2.push_back(1);
  }
  map<int,int> out;
  vector< vector<double> > deltaRs;

  // if (testmodule) cout << endl<< " Entry: " << iEntry <<endl;
  for (unsigned iv1 = 0; iv1 < v1.size(); ++iv1) {
    vector<double> v2p;
    if (testmodule) cout << "|";
    for (unsigned iv2 = 0; iv2 < v2.size(); ++iv2) {
      double dr12 = v1[iv1].DeltaR(v2[iv2]);
      v2p.push_back(dr12);
      if (testmodule) cout << Form(" %5.3f |",dr12);
    }
    if (testmodule) cout <<endl;
    deltaRs.push_back(v2p);
  }
  if (testmodule) cout <<endl;

  int x = min(v1.size(), v2.size());
  for (int it = 0; it < x; ++it) {
    double minR;
    pair<int,int> matchpair = FindMatrixMin(deltaRs,l1,l2,minR);
    if (testmodule) {
      int minv1 = matchpair.first;
      int minv2 = matchpair.second;
      cout << Form("V1:%d,pT: %5.2f,eta: %5.2f,phi: %5.2f | V2:%d,pT: %5.2f,eta: %5.2f,phi: %5.2f | dR: %5.3f",minv1,v1[minv1].Pt(),v1[minv1].Eta(),v1[minv1].Phi(),minv2,v2[minv2].Pt(),v2[minv2].Eta(),v2[minv2].Phi(),minR) << endl;
    }

    if (cut && minR > AlgodR) {
      if (testmodule) cout <<Form("******Cut at minR = %5.3f",minR) << endl;
      break;
    }
    maxDeltaR = minR;
    out.insert(matchpair);
    l1[matchpair.first] = 0;
    l2[matchpair.second] = 0;
  }
  return out;
}

map<int,vector<int> > AdvJetMatch(vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, double &maxDeltaR, bool cut, bool FilterSoftGen, bool FilterSoft) {
  bool testmodule = false;
  int nl1(0), nl2(0);
  vector<int> l1, l2;
  for (unsigned il = 0; il < v1.size(); ++il) {
    if ((v1[il].Pt() < JetPtThreshold && FilterSoftGen)||v1[il].Pt() == 0) l1.push_back(0);
    else {
      l1.push_back(1);
      ++nl1;
    }
  }
  for (unsigned il = 0; il < v2.size(); ++il) {
    if ((v2[il].Pt() < JetPtThreshold && FilterSoft)||v2[il].Pt() == 0) l2.push_back(0);
    else {
      l2.push_back(1);
      ++nl2;
    }
  }
  map<int,vector<int> > out;
  vector< vector< vector<double> > > deltaRs;

  // if (testmodule) cout << endl<< " Entry: " << iEntry <<endl;
  for (unsigned iv1 = 0; iv1 < v1.size(); ++iv1) {
    vector< vector<double> > v2p;
    if (testmodule) cout << "|";
    for (unsigned iv2a = 0; iv2a < v2.size(); ++iv2a) {
      vector< double > v2merge;
      if (testmodule) cout << "|";
      for (unsigned iv2b = 0; iv2b <= iv2a; ++iv2b) {
        TLorentzVector Mergedv2;
        if (iv2a == iv2b) Mergedv2 = v2[iv2a];
        else Mergedv2 = v2[iv2a] + v2[iv2b];
        double dr12 = v1[iv1].DeltaR(Mergedv2);
        v2merge.push_back(dr12);
        if (testmodule) cout << Form(" %5.3f |",dr12);
      }
      v2p.push_back(v2merge);
    }
    if (testmodule) cout <<endl;
    deltaRs.push_back(v2p);
  }
  if (testmodule) cout <<endl;

  // int x = min(v1.size(), v2.size());
  if (testmodule) cout <<Form("\nnl1 = %d, nl2 = %d \n",nl1,nl2);
  while((nl1)&&(nl2)) {
    int ml1(0), ml2a(0), ml2b(0);
    double minR = 0;
    // pair<int,vector<int> > pair;

    for (unsigned il1 = 0; il1 < l1.size(); ++il1) {
      if (l1[il1] == 0) continue;
      for (unsigned il2a = 0; il2a < l2.size(); ++il2a) {
        if (l2[il2a] == 0) continue;
        for (unsigned il2b = 0; il2b <= il2a; ++il2b) {
          if (l2[il2b] ==0) continue;
          if (deltaRs[il1][il2a][il2b] < minR || minR == 0) {
            minR = deltaRs[il1][il2a][il2b];
            ml1 = il1;
            ml2a = il2a;
            ml2b = il2b;
          }
        }
      }
    }


    if (testmodule) {
      cout << Form("V1:%d,pT: %5.2f,eta: %5.2f,phi: %5.2f | V2a:%d,pT: %5.2f,eta: %5.2f,phi: %5.2f | V2b:%d,pT: %5.2f,eta: %5.2f,phi: %5.2f | dR: %5.3f",ml1,v1[ml1].Pt(),v1[ml1].Eta(),v1[ml1].Phi(),ml2a,v2[ml2a].Pt(),v2[ml2a].Eta(),v2[ml2a].Phi(),ml2b,v2[ml2b].Pt(),v2[ml2b].Eta(),v2[ml2b].Phi(),minR) << endl;
    }

    if (cut && minR > AlgodR) {
      if (testmodule) cout <<Form("******Cut at minR = %5.3f",minR) << endl;
      break;
    }
    maxDeltaR = minR;
    vector<int> mergejets;
    mergejets.push_back(ml2a);
    mergejets.push_back(ml2b);
    out.insert( pair<int,vector<int> > (ml1, mergejets) );
    l1[ml1] = 0;
    l2[ml2a] = 0;
    l2[ml2b] = 0;
    nl1--;
    if (ml2a == ml2b) nl2-=1;
    else nl2-=2;
    if ((nl1 < 0 || nl2 < 0) && testmodule) cout << "!!!Lists over-substraction!!!!"<<endl;
  }
  return out;
}

#endif
