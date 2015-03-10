function [estim, asympAnalysis, bwXZ, bwYZ] = condTsallisDivergence(X, Y, Z, ...
  functionalParams, params)
% Estimates the Conditional Tsallis Divergence between X and Y given Z. Note
% that this is the same as I_T(XZ ; YZ). functionalParams should contain a field
% alpha.
% This estimator is for the setting when X,Y,Z ~ P_XYZ. If instead you have
% data from the joints XZ, YZ use klDivergence(XZ, YZ, ...) instead.
  XZ = [X Z];
  YZ = [Y Z];
  [estim, asympAnalysis, bwXZ, bwYZ] = tsallisDivergence(XZ, YZ, ...
    functionalParams, params);
end

