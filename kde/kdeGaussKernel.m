function K = kdeGaussKernel(X, Y, h)
% Returns the Kernel Matrix for a Gaussian Kernel of bandwidth h.
% X is an nxd matrix. K is an nxn matrix.
% If Y is nonempty then returns the gaussian kernel for XxY
% h is a column vector of size the dimension of the space - i.e size(X, 1).

  % Prelims
  d = size(X, 2); % dimensions
  if size(h, 1) == 1, h = h'; % if you get a row vector
  end
  if isscalar(h)
    h = h * ones(d, 1);
  end;

  if ~exist('Y', 'var') | isempty(Y)
    Y = X;
  end

  D2 = distSquared(X, Y, h);
  K = 1/(sqrt(2*pi)^d * prod(h)) * exp(-D2/2);

end

