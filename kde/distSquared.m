function D2 = distSquared(X, Y, h)
% distanceSquared: Calculates squared distance between two sets of points.
% If h is provided, this scales each dimension by the bandwidth

  if ~exist('h', 'var')
    h = 1;
  end

  [nX, dX] = size(X);
  [nY, dY] = size(Y);
  if dX ~= dY
    error('Dimensions of X, Y do not match')
  end

  X = bsxfun(@rdivide, X, h');
  Y = bsxfun(@rdivide, Y, h');

  D2 = (ones(nY, 1) * sum((X.^2)', 1))' + ...
    ones(nX, 1) * sum((Y.^2)',1) - ...
    2.*(X*(Y'));

  % Rounding errors occasionally cause negative entries in n2
  if any(any(D2<0))
    D2(D2<0) = 0;
  end

end

