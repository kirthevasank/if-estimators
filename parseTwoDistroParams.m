function params = parseTwoDistroParams(params, X, Y)
% X, Y are data from densities f and g respectively.
% params is a struct with other parameters (viz. whether to data split or not)

  numDims = size(X, 2);
  params = parseCommonParams(params, numDims);

end

