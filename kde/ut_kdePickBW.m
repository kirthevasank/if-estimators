% Unit test for kdePickBW.m

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
params.doBoundaryCorrection = true;
smoothness = 4;
bwLogBounds = []; %log( [0.5 0.51]);
[optBw, f] = kdePickBW(X, smoothness, params, bwLogBounds);
fprintf('Picked h = %f\n', optBw);

% Plot the density.
t = linspace(0,1,resolution)'; 
p = f(t);
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
pause;

% Test 2: Uniform 2D
X = rand(2000,2); X = 0.25 + 0.5*X;
[optBw, f] = kdePickBW(X, smoothness, params, bwLogBounds);
fprintf('Picked h = %f\n', optBw);
t = linspace(0,1,resolution); [T1, T2] = meshgrid(t);
p = f([T1(:), T2(:)]);
P = reshape(p, resolution, resolution);
figure;
surfc(T1, T2, P);
area = numerical_2D_integration(P, T1, T2);
fprintf('Area under estimated curve: %f\n', area);
pause,

% Test 3: A 1D distribution
gamma = 10;
Z = rand(N, 1+gamma); B = double(rand(N, 1) < 0.5);
  % X = 0.5*(Z(:,1) + max(Z(:,2:end), [], 2));
  % trueDensity = @(t) 2*( (t<0.5).*(t/0.5).^gamma + (t>=0.5).*(1-(t/0.5-1).^gamma));
X = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
trueDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
% First estimate the density
[optBw, f] = kdePickBW(X, smoothness, params, bwLogBounds);
fprintf('Picked h = %f\n', optBw);
% Now plot them
th = linspace(0,1,resolution)';
truth = trueDensity(th);
figure;
plot(th, f(th), 'b', th, trueDensity(th), 'r'); hold on;
plot(X, 0.2*rand(size(X)), 'kx');
titlestr = sprintf('Estimated(b) vs True(r)\nh');
title(titlestr);
area = numerical_1D_integration(th, f(th));
fprintf('Area under estimated curve: %f\n', area);
pause;

