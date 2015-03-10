function [estim, asympAnalysis, bwX, bwY] = ...
  fAlphaGBeta(X, Y, functionalParams, params)
  % Estimates the integral \int f^alpha g^beta where beta = 1-alpha. X comes
  % from f and Y comes from g. functionalParams should contain the field alpha.
  params = parseTwoDistroParams(params, X, Y);

  infFunX = @(u,v) fAlphaGBetaInfFunX(u, v, functionalParams.alpha);
  infFunY = @(u,v) fAlphaGBetaInfFunY(u, v, functionalParams.alpha);
  asympVarFun = @(u,v,w,z) ...
    fAlphaGBetaAsympVar(u, v, w, z, functionalParams.alpha);
  [estim, asympAnalysis, bwX, bwY] = ...
    getTwoDistroInfFunAvgs(X, Y, infFunX, infFunY, asympVarFun, params);

end


function infFunVals = fAlphaGBetaInfFunX(densXatX, densYatX, alpha)
% densXatX and densYatX are the densities of X and Y respectively (or their
% estimates) at the X points.
  infFunVals = alpha * (densYatX ./ densXatX).^(1-alpha);
end


function infFunVals = fAlphaGBetaInfFunY(densXatY, densYatY, alpha)
% densXatY and densYatY are the densities of X and Y respectively (or their
% estimates) at the Y points.
  infFunVals = (1-alpha) * (densXatY ./ densYatY).^alpha;
end


function [asympVarX, asympVarY] = fAlphaGBetaAsympVar( ...
  densXatX, densXatY, densYatX, densYatY, alpha)
% asympVarX and asympVarY are the asymptotic variances of the X and Y influence
% functions.
  phiAB = estimFAGBTemp(densXatX, densXatY, densYatX, densYatY, alpha);
  phi2Am12B = estimFAGBTemp(densXatX, densXatY, densYatX, densYatY, 2*alpha-1);
  phi2A2Bm1 = estimFAGBTemp(densXatX, densXatY, densYatX, densYatY, 2*alpha);

  asympVarX = alpha^2 * (phi2Am12B - phiAB^2); 
  asympVarY = (1-alpha)^2 * (phi2A2Bm1 - phiAB^2);
end

