function K = kdeLegendreKernel(X, C, h, order)
% Returns the value of the kernel evaluated at the points X centred at C and
% with bandwidth h.
% Inputs
% X : nxd data matrix
% C : mxd centre matrix. If empty is initialized to zero(1, d)
% h : the bandwidth of the kernel
% order : order of the kernel
% Ouputs
% K : The nxm kernel matrix where K(i,j) = k(X(i,:), C(j,:))
%%% Warning: make sure mxn < 1e6 to avoid crashing
  
  % Prelims
  numDims = size(X, 2);

  if isempty(C)
    C = zeros(1, numDims);
  end

  numData = size(X, 1);
  numCentres = size(C, 1);

  K = ones(numData, numCentres);
  for d = 1:numDims
    K = K .* kernel1D( X(:, d), C(:, d), h, order);
  end
end

% 1 Dimensional Kernel. The d-dimensional kernel is the product kernel.
function ret = kernel1D(x, c, h, order)
% Same as above but now x and c are 1 dimensional (d=1)

  numCentres = size(c, 1);
  % u is a numData x numCentres matrix, u_ij = (x_i - c_j)/h
  u = bsxfun(@minus, repmat(x, 1, numCentres), c')/h;

  ret = zeros(size(u));
  for m = 0:2:order
    % only need to iterate through even m since legPoly(0,m) = 0 for m odd
    ret = ret + legPoly(0, m) * legPoly(u, m);
  end
  % Finally check if u is within the domain of the kernel and divide by h.
  ret = ret .* double(abs(u) < 1)/h ;
end

