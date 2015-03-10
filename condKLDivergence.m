function [estim, asympAnalysis, bwXZ, bwYZ] = condKLDivergence(X, Y, Z, ...
  functionalParams, params)
% Estimates the Conditional KL Divergence between X and Y conditioned on Z. I.e.
% KL(X||Y | Z). Note that this is just the KL between XZ and YZ
% This estimator is for the setting when X,Y,Z ~ P_XYZ. If instead you have
% data from the joints XZ, YZ use klDivergence(XZ, YZ, ...) instead.
  XZ = [X Z];
  YZ = [Y Z];
  [estim, asympAnalysis, bwXZ, bwYZ] = klDivergence(XZ, YZ, ...
    functionalParams, params);
end

