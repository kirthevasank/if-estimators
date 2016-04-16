function [estim, asympAnalysis, bwX, bwY, bwXY] = shannonMI(X, Y, ...
  functionalParams, params)
% Estimates the Shannon Mutual Information I(X;Y) = \int fXY log(fXY/(fX*fY)).
% X and Y can have different dimensions (columns) but should have the same
% number of rows.

  [hXY, asympXY, bwXY] = shannonEntropy([X Y], functionalParams, params);
  [hX, asympX, bwX] = shannonEntropy(X, functionalParams, params);
  [hY, asympY, bwY] = shannonEntropy(Y, functionalParams, params);

  % The Estimator
  estim = max(hX + hY - hXY, 0);

  % Asymptotic Variance
  if ~isempty(asympX) 
    n = size(X, 1);
    asympVar = asympXY.asympVar + asympX.asympVar + asympY.asympVar;
    asympAnalysis = getAsympAnalysis(estim, asympVar, params.alpha, n);
  else
    asympAnalysis = [];
  end

end

