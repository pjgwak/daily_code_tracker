#include <cmath>

double PtThreshold(double eta)
{
  const double x = eta;

  const int nPts = 10;
  const double mEta[nPts] = {-2.4, -1.7, -1.3, -1.3, -1.0,
                             1.0, 1.3, 1.3, 1.7, 2.4};

  const double mPt[nPts] = {1.0, 1.0, 1.53, 2.1, 3.3,
                            3.3, 2.1, 1.53, 1.0, 1.0};

  int iseg = -1;

  for (int i = 0; i < nPts - 1; i++)
  {
    if (x >= mEta[i] && x < mEta[i + 1])
    {
      iseg = i;
      break;
    }
  }

  if (x == mEta[nPts - 1])
    iseg = nPts - 2;

  if (iseg < 0)
    return 999999.;

  double mSlope = (mPt[iseg + 1] - mPt[iseg]) /
                  (mEta[iseg + 1] - mEta[iseg]);

  double mPtTh = mSlope * (x - mEta[iseg]) + mPt[iseg];

  return mPtTh;
}
bool IsAcceptanceQQ(double pt, double eta)
{
  double th = PtThreshold(eta);
  return pt >= th;
}

