function [estim, asympAnalysis, bwX, bwY] = tsallisDivergence(X, Y, ...
  functionalParams, params)
% Estimates the Tsallis Divergences between f and g where X comes from f and Y
% comes from g. functionalParams should contain a field alpha.

  [estim1, asymp1, bwX, bwY] = fAlphaGBeta(X, Y, functionalParams, params);
  estim = 1/(functionalParams.alpha-1) * (estim1 - 1); 

  if ~isempty(asymp1)
    % Now modify the terms in the asymptotic analysis
    N = size(X, 1) + size(Y, 1);
    asympAnalysis.asympVar = asymp1.asympVar / (functionalParams.alpha-1)^2;
    asympAnalysis.asympStd = sqrt(asympAnalysis.asympVar);
    width = norminv(1-params.alpha/2) * asympAnalysis.asympStd/sqrt(N);
    asympAnalysis.confInterval(1) = estim - width;
    asympAnalysis.confInterval(2) = estim + width;
  else
    asympAnalysis = [];
  end

end

