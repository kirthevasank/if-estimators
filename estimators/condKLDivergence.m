function [estim, asympAnalysis, bwXZ, bwYZ] = condKLDivergence(X, Y, ZX, ZY, ...
  functionalParams, params)
% Estimates the Conditional KL Divergence between X and Y conditioned on Z. 
% If you have data (X,Y,Z) from the joint P_XYZ, send X = X, Y = Y, ZX = Z and pass
% an empty matrix for ZY. Otherwise if you have data from the marginals (X,Z) form
% P_XZ and (Y,Z) from P_YZ then pass them in X, ZX, Y, ZY respectively. 

  if isempty(ZY),
    ZY = ZX;
  end

  XZ = [X ZX];
  YZ = [Y ZY];
  [estim, asympAnalysis, bwXZ, bwYZ] = klDivergence(XZ, YZ, ...
    functionalParams, params);
end

