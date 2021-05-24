#include <algorithm>
#include <cmath>
#include <vector>

#include "quantile.h"

double Lerp(double v0, double v1, double t)
{
    return (1 - t)*v0 + t*v1;
}


void quantile(double* x, double *probs, double* out, int n, int k) {
  
  std::sort(x, x+n);

  // Loop through the quantiles of interest.
  for (size_t i = 0; i < k; ++i)
  {
    double poi = Lerp(-0.5, n - 0.5, probs[i]);
    
    size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
    size_t right = std::min(int64_t(std::ceil(poi)), int64_t(n - 1));
    
    out[i] = Lerp(x[left], x[right], poi - left);
  }
}
