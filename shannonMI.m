function [estim, asympAnalysis, bwX, bwY, bwXY] = shannonMI(X, Y, ...
  functionalParams, params)
% Estimates the Shannon Mutual Information I(X;Y) = \int fXY log(fXY/(fX*fY)).
% X and Y can have different dimensions (columns) but should have the same
% number of rows.

  [hXY, asympXY, bwXY] = shannonEntropy([X Y], functionalParams, params);
  [hX, asympX, bwX] = shannonEntropy(X, functionalParams, params);
  [hY, asympY, bwY] = shannonEntropy(Y, functionalParams, params);

  % The Estimator
  estim = hX + hY - hXY;

  % Asymptotic Variance
  if params.doAsympAnalysis
    n = size(X, 1);
    asympAnalysis.asympVar = asympXY.asympVar + asympX.asympVar + ...
                               asympY.asympVar;
    asympAnalysis.asympStd = sqrt(asympAnalysis.asympVar);
    width = norminv(1-params.alpha/2) * asympAnalysis.asympStd / sqrt(n);
    asympAnalysis.confInterval(1) = estim - width;
    asympAnalysis.confInterval(2) = estim + width;
  end

end

