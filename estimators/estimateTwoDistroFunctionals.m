function [estim, asympAnalysis, bwX, bwY] = estimateTwoDistroFunctionals...
  (X, Y, ZX, ZY, functional, functionalParams, params)
% X, Y are data from densities f and g respectively.
% functional is a string specifying which functional to estimate. See below.
% functionalParams specifies any parameters specific to the functional (viz.
%   alpha for the alpha-divergences)
% params is a struct with other parameters (viz. whether to data split or not)

  % prelims
  numX = size(X, 1);
  numY = size(Y, 1);
  numDims = size(X, 2);

  if ~exist('params', 'var')
    params = struct;
  end
  if ~exist('functionalParams', 'var')
    functionalParams = struct;
  end

  switch functional

    case 'klDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        klDivergence(X, Y, functionalParams, params);

    case 'renyiDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        renyiDivergence(X, Y, functionalParams, params);

    case 'tsallisDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        tsallisDivergence(X, Y, functionalParams, params);

    case 'hellingerDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        hellingerDivergence(X, Y, functionalParams, params);

    case 'chiSqDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        chiSqDivergence(X, Y, functionalParams, params);

    case 'condKLDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        condKLDivergence(X, Y, ZX, ZY, functionalParams, params);

    case 'condTsallisDiv'
      [estim, asympAnalysis, bwX, bwY] = ...
        condTsallisDivergence(X, Y, ZX, ZY, functionalParams, params);

    case 'condTsallisMI'

    case 'condRenyiDiv'

    case 'condRenyiMI'

  end

end
