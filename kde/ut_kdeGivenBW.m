% Unit test for kde01.m

clear all;
close all;
N = 1000;
resolution = 100;
doLogit = false;

% prelims
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/

% Test 1: Beta distribution
% =========================

fprintf('Test 1: 1D Beta Distribution\n');
th1 = [9;40];
th2 = [30;5];
p1 = 0.7;
p2 = 1 - p1;
D1 = dirichlet_sample(th1', N); D1 = D1(:,1);
D2 = dirichlet_sample(th2', N); D2 = D2(:,1);
Z = double(rand(N,1) < p1);
X = Z .* D1 + (1-Z) .* D2;
% estimate the density
h = 0.05;
kde = kdeGivenBW(X, h, 4);

[~, f] = kde01(X, [], doLogit);
% Plot the density.
t = linspace(0,1,resolution)'; 
p = kde(t);
% obtain true density 
true_density = ...
  p1 * t.^(th1(1)-1) .* (1-t).^(th1(2)-1) / beta(th1(1), th1(2))  + ...
  p2 * t.^(th2(1)-1) .* (1-t).^(th2(2)-1) / beta(th2(1), th2(2));
plot(t, p, 'b', t, true_density, 'r'); hold on,
plot(X, 0.2*rand(size(X)), 'kx');
titlestr = sprintf('Estimated(b) vs True(r)\nh');
title(titlestr);
area = numerical_1D_integration(t, p);
fprintf('Area under estimated curve: %f\n', area);

%Test 2: Uniform 2D
fprintf('\n\nTest 2: 2D Distribution\n');
X = rand(2000,2); X = 0.25 + 0.5*X;
% [~, f] = kde01(X, [], doLogit);
f = kdeGivenBW(X, h, 4);
t = linspace(0,1,resolution); [T1, T2] = meshgrid(t);
p = f([T1(:), T2(:)]);
P = reshape(p, resolution, resolution);
figure;
surfc(T1, T2, P);
area = numerical_2D_integration(P, T1, T2);
fprintf('Area under estimated curve: %f\n', area);
area,


