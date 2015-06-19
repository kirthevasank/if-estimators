% This is a unit test for legPoly.
% Compare the outputs of legPoly vs matlab's legendre implementation. For orders
% less than 4 use our version since it doesn't unnecessarily compute all
% legendre functions.

x = 2*rand(4, 3) - 1; % 5 random points in [-1, 1]
X = x;
x(10) = -1.1;
x(7) = 1.2;

for order = 1:6
  val1 = legPoly(x, order);
  val2 = legendre(order, X);
  val2 = sqrt( (2*order+1)/2)* val2(1, :)'; val2 = reshape(val2, size(X));
  val2(10) = 0; val2(7) = 0;
  val1, val2,
  fprintf('Order: %d, error: %f, x: %s\n',order, norm(val1 - val2), mat2str(x));
end
