function [estim, asympAnylsis] = estimateOneDistroFunctionals...
  (X, functional, functionalParams, params)
% X is data from densities f .
% functional is a string specifying which functional to estimate. See below.
% functionalParams specifies any parameters specific to the functional (viz.
%   alpha for the alpha-divergences) % params is a struct with other parameters (viz. whether to data split or not) 
  % prelims

  numX = size(X, 1);
  numDims = size(X, 2);

  if ~exist('params', 'var') | isempty(params)
    params = struct;
  end
  if ~exist('functionalParams', 'var') | isempty(functionalParams)
    functionalParams = struct;
  end

  switch functional

    case 'shannonEntropy'
      [estim, asympAnylsis] = ...
        shannonEntropy(X, functionalParams, params);
% 
%     case 'pAlpha'
%       [estim, asympAnylsis] = ...
%         pAlpha(X, functionalParams, params);

      
  end

end

