Progress Report
========================================================
author: Kyle , Tangwei, Manze
date: 3/17/2016

Prediction Interval
========================================================
Let $C_{n}(x) \subset \mathbb{R}^{1}$ be an interval prediction s.t. for all $x$

$$
C_{n}(x) = C((X_1, Y_1), ..., (X_{n}, Y_n), x)
$$
and 
$$
P(Y_{n+1} \in C_{n}(X_{n+1})) \geq 1-\alpha
$$

We want a consistent method of evaluating both the interval and the point prediction.




Old Idea
========================================================
$$
PRMSPE = \frac{\| Y-\hat{Y} \|}{\| Y \|} + \frac{\sum_{j=1}^{k}(C_{n}^{u}(X_{n+j}) - C_{n}^{l}(X_{n+j}))}{\| Y \|}
$$

Penalize against the average size of the interval 

Problem
========================================================

$$
\frac{\sum_{j=1}^{k}(C_{n}^{u}(X_{n+j}) - C_{n}^{l}(X_{n+j}))}{\| Y \|}
$$

This penalty is way too large relative to the size of $RMSPE$. It also doesn't
take into account the number of times $Y_{n+1} \in C_{n}(X_{n+1})$ in the test set.

Solution
========================================================
- Gneiting, Raferty (2007)
- Defined the idea of Proper Scoring Rules
- Many nice theoretical properties in abstract probability spaces
- Analogous to loss functions from estimation literature.
- Adapting this idea to evaluating intervals, we get:


Solution
========================================================
- lower limit $p_l$ 
- upper limit $p_u$, 
- Goal: prediction interval at level $\alpha$ 

$$
S(C_n(x), y)_{\alpha} = (p_u - p_l) + \frac{2}{\alpha}(p_l - y)\mathbb{1}(y < p_l) + \frac{2}{\alpha}(y - p_u)\mathbb{1}(y > p_u)
$$

- Intuition: We want to have very short prediction intervals 
- Cover the realization $y$, $p_l \leq y \leq p_u$.


Solution
========================================================

$$
S \approx \text{ width of prediction interval} + \text{ adjust for coverage}
$$

$$
PRMSPE = RMSPE + \frac{\bar{S}}{\|Y\|}
$$

Conjecture:
$S$ allows for meaningful comparisons of models with prediction intervals of 
different probability coverages.

Performance
========================================================



```
Error in parse(text = as.character(lab)) : 
  <text>:1:7: unexpected numeric constant
1: Model 2
          ^
```
