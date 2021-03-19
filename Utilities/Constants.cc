#include <vector>
#include <cmath>

using namespace std;

class ScaleToLumi {
public:
  ScaleToLumi(double lumi_ = 1.) {
    SetLumi(lumi_);
  }

  vector<double> ScaleVectors(double lumi_ = -1.) {
    SetLumi(lumi_);
    vector<double>out{ScaleFL, ScaleLL, ScaleBG};
    return out;
  }

  void SetLumi(double lumi_) {
    if (lumi_ != -1) {
      Lumi = lumi_;
      ScaleFL = xsecFL * Lumi / nEvtFL;
      ScaleLL = xsecLL * Lumi / nEvtLL;
      ScaleBG = xsecBG * Lumi / nEvtBG;
    }
  }

  double Lumi;
  double ScaleFL, ScaleLL, ScaleBG;

  const double nEvtFL = 81215.;
  const double nEvtLL = 81262;
  const double nEvtBG = 900789;
  const double xsecFL = 1.218;
  const double xsecLL = 1.229;
  const double xsecBG = 127100.;
  const double cmsLumi = 146.45;
  const double cmsLumi18 = 63.67;
  const double cmsLumi17 = 44.98;
  const double cmsLumi16 = 37.80;
  const double cmsLumi15 = 3.80;
  const double lhcLumi18 = 67.86;
};

namespace Constants {
  const double nEvtFL = 81215.;
  const double nEvtLL = 81262;
  const double nEvtBG = 900789;
  const double xsecFL = 1.218;
  const double xsecLL = 1.229;
  const double xsecBG = 127100.;
  const double Lumi = 1.;
  const double cmsLumi = 146.45;
  const double cmsLumi18 = 63.67;
  const double cmsLumi17 = 44.98;
  const double cmsLumi16 = 37.80;
  const double cmsLumi15 = 3.80;
  const double lhcLumi18 = 67.86;
  const double ScaleFL = xsecFL * Lumi / nEvtFL;
  const double ScaleLL = xsecLL * Lumi / nEvtLL;
  const double ScaleBG = xsecBG * Lumi / nEvtBG;
}
