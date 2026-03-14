#include <iostream>

void err2(float ptLow = 14, float ptHigh = 20, float yLow = 1.6, float yHigh = 2.4)
{
  std::cout
      << "[INFO] err2.C is disabled in 03/13_OO_no_pee because the event-by-event error-dependent "
      << "cuts and PDFs were removed. "
      << "No err model is built for pt=[" << ptLow << ", " << ptHigh
      << "], |y|=[" << yLow << ", " << yHigh << "]."
      << std::endl;
}
