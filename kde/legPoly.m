function val = legPoly(V, order)
% Returns the legendre polynomial of order 'order' at each element in V.
% V is an numPts x 1 vector
% This is sqrt( (2*order + 1)/ 2) *legendrePolynomial. The initial constant is
% needed for the kde

  % First identify those that are within [-1, 1];
  v = V(:);
  legalIdxs = abs(v) < 1;
  u = v(legalIdxs);
 
  switch order
    case 0
      valU = ones(size(u));

    case 1
      valU = u;

    case 2
      valU = 1/2 * (3*u.^2 - 1);

    case 3
      valU = 1/2 * (5*u.^3 - 3*u);

    case 4
      valU = 1/8 * (35*u.^4 - 30*u.^2 + 3);

    otherwise
      temp = legendre(order, u);
      valU = temp(1, :)';
%       valU = reshape(valU, size(u));
  end

  % Now fill in the rest with zeros
  val = zeros(size(v));
  val(legalIdxs) = valU;
  val = reshape(val, size(V));
  
  % finally multiply by sqrt( (2*order + 1)/ 2)
  val = sqrt( (2*order + 1)/2 ) * val;
end

