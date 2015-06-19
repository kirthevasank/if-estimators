function [estim, asympAnalysis, bwXZ, bwYZ] = condTsallisDivergence(X, Y, ZX, ZY, ...
  functionalParams, params)
% Estimates the Conditional Tsallis Divergence between X and Y given Z.
% If you have data (X,Y,Z) from the joint P_XYZ, send X = X, Y = Y, ZX = Z and pass
% an empty matrix for ZY. Otherwise if you have data from the marginals (X,Z) form
% P_XZ and (Y,Z) from P_YZ then pass them in X, ZX, Y, ZY respectively.

  if isempty(ZY), 
  % In this case, we have data X,Y,Z from the joint P_XYZ.
    ZY = ZX;
  end

  XZ = [X ZX];
  YZ = [Y ZY];
  [estim, asympAnalysis, bwXZ, bwYZ] = tsallisDivergence(XZ, YZ, ...
    functionalParams, params);
end

