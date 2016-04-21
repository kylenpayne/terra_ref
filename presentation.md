Terra Ref Meeting
========================================================
author: Kyle N. Payne
date: 

Comparing Two Intervals
========================================================
 
- Same Model Predictions, implies prediction error part should 
be the same
- Different Intervals


Interval 1
========================================================

$$
\hat{Y} \pm z_{\alpha/2}\sqrt{\sqrt{\sigma} + s}
$$

Interval 2
========================================================

$$
\hat{Y} \pm z_{\alpha/2, boot}\sqrt{\sqrt{\sigma} + s}
$$

$$
z_{\alpha/2, boot}
$$ 
is a bootstrap quantile estimate that should, in principal, be inaccurate for a small number of 
bootstrap replicates and get better than the normal approximation.

Consequence
========================================================

Overall Score should get better with higher number of bootstrap samples.





```
Error in results[k, ] <- c(p_norm, p_boot) : 
  number of items to replace is not a multiple of replacement length
```
