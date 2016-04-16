function [estim, asympAnalysis, bwX, bwY] = renyiDivergence(X, Y, ...
  functionalParams, params)
% Estimates the Reniy-alpha Divergence between f and g where X comes from f and
% Y from g. functionalParmas should contain a field alpha.

  % To apply delta method
  g = @(t) 1/(functionalParams.alpha-1) * log(t);
  gPrime = @(t) 1/(functionalParams.alpha -1) * 1/t;
    
  [estim1, asymp1, bwX, bwY] = fAlphaGBeta(X, Y, functionalParams, params);
  estim = g(estim1);

  if ~isempty(asymp1) 
    % Asymptotic Analysis
    N = size(X, 1) + size(Y, 1);
    asympAnalysis.asympVar = asymp1.asympVar * (gPrime(estim1))^2;
    asympAnalysis.asympStd = sqrt(asympAnalysis.asympVar);
    width = norminv(1-params.alpha/2) * asympAnalysis.asympStd/sqrt(N);
    asympAnalysis.confInterval(1) = estim - width;
    asympAnalysis.confInterval(2) = estim + width;
  else
    asympAnalysis = [];
  end

end

