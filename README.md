# chick2016

## Goals

The goals of my analysis are to understand the distinctions between 1) the "LOD difference" statistic (my terminology) that Chick et al 2016 use and 2) the Baron-Kenny approach, as described in vanderWeele's text (Chapter 2).

It's not obvious to me, at first glance, that there is an equivalence between the two approaches. vanderWeele describes the Baron & Kenny approach as requiring two regressions, Y ~ M + A and M ~ A, for exposure A, mediator M, and outcome (target) Y.

Chick et al. (2016) use a LOD difference statistic. This requires four regressions. First, they calculate the LOD by comparing likelihoods for these two models:

Y ~ A v Y ~ 1

Then, they fit two additional models:

Y ~ A + M and Y ~ 1 + M

Each pair of models above gives a LOD score. They then calculate the difference of these two LOD scores as their summary statistic.








## Chick, Munger, et al. (2016) data and scripts

Data and R scripts downloaded, on May 4, 2019, from: https://churchill-lab.jax.org/website/data/Chick_Nature_2016/

