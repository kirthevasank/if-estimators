function params = parseOneDistroParams(params, X)
% X is data from the density f.
% params is a struct with some parameters for the estimation. 
% (viz. whether to data split or not)

  numDims = size(X, 2);
  params = parseCommonParams(params, numDims);

end

