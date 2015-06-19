function [estim, asympAnalysis, bw] = shannonEntropy(X, functionalParams, params)
% This estimates the shannon entropy -\int plog(p)
  params = parseOneDistroParams(params, X);
  [estim, asympAnalysis, bw] = getInfFunAvgs(X, @entropyInfFun, ...
    @entropyAsympVar, params);
end

function infFunVals = entropyInfFun(densX)
% densX is a vector of the values of the density (or its estimate) at the data
% points X.
  infFunVals = -log(densX);
end

function asympVar = entropyAsympVar(densX)
  logDensX = log(densX);
  asympVar = mean(logDensX.^2) - mean(logDensX)^2;
end

