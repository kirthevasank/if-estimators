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
- Since all estimators are nonparametric, they work well only in under 4-6 
  dimensions depending on the problem.

### Installation and Getting Started
- Run ifSetup.m with the path to the installation as the argument. This adds all
  subdirectories to the matlab workspace and you are good to go.
- To get started, see the demos directory. demos.m illustrates a very simple use
  case.
- demo1.m ... demo4.m illustrate all functionals in different settings. We re
- See parseCommonParams for setting of hyperparameters. We recommend using the
  default settings unless you are very familiar with the paper.

### Citation
If you use this library in your academic work please cite our paper: "Influence
Functions for Machine Learning: Nonparametric Estimators for Entropies,
Divergences and Mutual Informations", Kirthevasan Kandasamy, Akshay
Krishnamurthy, Barnabas Poczos, Larry Wasserman, James Robins.

### License
This software is released under GNU GPL v3(>=) License. Please read LICENSE.txt for
more information.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

"Copyright 2015 Kirthevasan Kandasamy"


#### Acknowledgements
We thank Zoltán Szabó from UCL for pointing out some bugs in a previous version of
this software.
* Some of the above estimators also appear in the ITE toolbox: 
https://bitbucket.org/szzoli/ite/
* For questions/ bug reports please email kandasamy@cs.cmu.edu


