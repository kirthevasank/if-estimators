function [estim, asympAnalysis, bwX, bwY] = hellingerDivergence(X, Y, ...
  functionalParams, params)
% Estimates the Hellinger Divergence between f and g where X comes from f and Y
% comes from g.
  if isempty(functionalParams), functionalParams = struct;
  end
  functionalParams.alpha = 0.5;

  [estim1, asymp1, bwX, bwY] = fAlphaGBeta(X, Y, functionalParams, params);
  estim = 1 - estim1;
  if ~isempty(asymp1)
    asympAnalysis = asymp1;
    asympAnalysis.confInterval(1) = 1 - asymp1.confInterval(2);
    asympAnalysis.confInterval(2) = 1 - asymp1.confInterval(1);
  else
    asympAnalysis = [];
  end

end
