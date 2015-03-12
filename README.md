## Nonparametric Estimation using Influence Functions
This is a Matlab Implementation of our Nonparametric Estimators using Influence
Functions. For more details, read our paper on Arxiv:
http://arxiv.org/abs/1411.4342

### Functionals \& Features
- Currently, our implementation covers the following functionals,
    Entropies: Shannon Entropy, Conditional Shannon Entropy
    Mutual Informations: Shannon MI, Conditional Shannon MI
    Divergences: Chi-squared, Hellinger, KL Divergence, Tsallis, Renyi, 
      Conditional KL, Conditional Tsallis
- In addition to the estimators, we also produce asymptotic confidence intervals
  for the estimators.
- Since all estimators are nonparametric, they work well only in under 5-7
  dimensions depending on the problem.

### Installation and Getting Started
- Just add this directory to your Matlab workspace and you are good to go.
- We have created Unit Tests in the files utX.m, utXY.m, utXgivenY.m,
  utXYgivenZ.m. You can use them as demos to learn how to use this library.
- For questions/ bug reports please email kandasamy@cs.cmu.edu 

### Citation
If you use this library in your academic work please cite our paper: "Influence
Functions for Machine Learning: Nonparametric Estimators for Entropies,
Divergences and Mutual Informations", Kirthevasan Kandasamy, Akshay
Krishnamurthy, Barnabas Poczos, Larry Wasserman, James Robins.

### Matlab ITE Toolbox
Our estimators are to be made available as part of the Matlab's ITE toolbox
(https://bitbucket.org/szzoli/ite/).

### License
This software is released under the GNU GPL License (version 3). See LICENSE.txt
for more details.

