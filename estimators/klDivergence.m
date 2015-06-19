function [estim, asympAnalysis, bwX, bwY] = klDivergence(X, Y, ...
  functionalParams, params)
% Estimates the KL Divergence between f and g where X comes from f and Y comes
% from g. 
  params = parseTwoDistroParams(params, X, Y);
  [estim, asympAnalysis, bwX, bwY] = ...
    getTwoDistroInfFunAvgs(X, Y, @klInfFunX, @klInfFunY, @klAsympVar, params);
end


function infFunVals = klInfFunX(densXatX, densYatX)
% densXatX and densYatX are the densities of X and Y respectively (or their
% estimates) at the X points.
  infFunVals = log( densXatX ./ densYatX );
end


function infFunVals = klInfFunY(densXatY, densYatY)
% densXatY and densYatY are the densities of X and Y respectively (or their
% estimates) at the Y points.
  infFunVals = 1 - densXatY ./ densYatY;
end


function [asympVarX, asympVarY] = klAsympVar( ...
  densXatX, densXatY, densYatX, densYatY)
% asympVarX and asympVarY are the asymptotic variances of the X and Y influence
% functions.
  asympVarY = estimFAGBTemp(densXatX, densXatY, densYatX, densYatY, 2);
  % Compute the following to estimate asympVarX
  fgX = densXatX ./ densYatX;
  fgY = densXatY ./ densYatY;
  logfgX = log(fgX);
  t1 = mean(logfgX.^2 + 2*logfgX) + mean(-fgY .* log(fgY) );
  t2 = 1 + mean(logfgX) - mean(fgY);
  asympVarX = t1 - t2^2;
end

