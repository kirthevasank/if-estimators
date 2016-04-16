function [estim, asympAnalysis, bwXZ, bwYZ, bwXYZ, bwZ] = ...
  condShannonMI(X, Y, Z, functionalParams, params)
% Estimates the Shannon Mutual Information between X and Y conditioned on Z.
% I.e. I(X; Y | Z). X, Y, Z could have different number of dimensions (columns)
% but should have the same number of rows.

  [hXYZ, asympXYZ, bwXYZ] = shannonEntropy([X Y Z], functionalParams, params);
  [hXZ, asympXZ, bwXZ] = shannonEntropy([X Z], functionalParams, params);
  [hYZ, asympYZ, bwYZ] = shannonEntropy([Y Z], functionalParams, params);
  [hZ, asympZ, bwZ] = shannonEntropy(Z, functionalParams, params);

  % The Estimator
  estim = hXZ + hYZ - hXYZ - hZ;

  % Asymptotic Variance
  if ~isempty(asympXYZ)
    n = size(X, 1);
    asympVar = asympXZ.asympVar + asympYZ.asympVar + ...
      asympXYZ.asympVar + asympZ.asympVar;
    asympAnalysis = getAsympAnalysis(estim, asympVar, params.alpha, n);
  else
    asympAnalysis = [];
  end

end

