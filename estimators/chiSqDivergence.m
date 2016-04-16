function [estim, asympAnalysis, bwX, bwY] = chiSqDivergence(X, Y, ...
  functionalParams, params)
% Estimates the Chi-Squared Divergence between f and g where X comes from f and
% Y comes from g.

  if isempty(functionalParams), functionalParams = struct;
  end
  functionalParams.alpha = 2;

  [estim1, asymp1, bwX, bwY] = fAlphaGBeta(X, Y, functionalParams, params);
  estim = estim1 - 1;
  if ~isempty(asymp1)
    asympAnalysis.confInterval(1) = asympAnalysis.confInterval(1) - 1;
    asympAnalysis.confInterval(2) = asympAnalysis.confInterval(2) - 1;
  else
    asympAnalysis = [];
  end

end

