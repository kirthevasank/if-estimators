function [estim, asympAnalysis, bwXY, bwY] = condShannonEntropy(X, Y, ...
  functionalParams, params)
% This estimates the entropy of X conditioned on Y.

  [hXY, asympXY, bwXY] = shannonEntropy([X Y], functionalParams, params);
  [hY, asympY, bwY] = shannonEntropy(Y, functionalParams, params);

  % Estimator
  estim = hY - hXY;

  % Asymptotic Variance
  if ~isempty(asympXY) 
    n = size(X, 1);
    asympVar = asympXY.asympVar + asympY.asympVar;
    asympAnalysis = getAsympAnalysis(estim, asympVar, params.alpha, n);
  else
    asympAnalysis = [];
  end

end

