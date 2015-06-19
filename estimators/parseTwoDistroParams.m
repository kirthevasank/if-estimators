function params = parseTwoDistroParams(params, X, Y)
% X, Y are data from densities f and g respectively.
% params is a struct with other parameters (viz. whether to data split or not)

  [n, numDims] = size(X);
  m = size(Y, 2);
  params = parseCommonParams(params, numDims, min(m, n));

end

