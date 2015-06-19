function demo4
% Unit tests for functionals on X, Y where they come from a joint distribution
% X, Y

  close all;
  clear all;

  fprintf('\nSome demos on conditional functionals of two distribution.\n');

  functionals = {'condShannonMI', 'condKLDiv', 'condTsallisDiv'};
  tests = {'1D-UnifUnifUnif', 'Indep-Gaussians', 'Gaussian'};

  % This is for storing parameters specific to the functional (E.g. alpha for the
  % alpha-divergences.) 
  functionalParams = struct;
  functionalParams.alpha = 0.8; % for the Alpha divergences

  % params is for storing the various parameters for estimation. For the most part,
  % we recommend using the default parameters (see estimators/parseCommonParams.m).
  params = struct;
  % If you also need to obtain asymptotic confidence sets set doAsympAnalysis to true
  % (set to false by default) and set the alpha level (i.e. alpha = 0.05 for a 
  % 95% confidence set).
  params.alpha = 0.05;
  params.doAsympAnalysis = true;

  % Test 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for testIdx = 1:numel(tests)

    % First generate the data
    if testIdx == 1
      N = 2000; 
      X = rand(N, 1);
      Y = rand(N, 1);
      Z = rand(N, 1);
      trueVals = zeros(3, 1);

    elseif testIdx == 2
      d1 = 1; d2 = 1; d3 = 4;
      N = 5000; 
      [C1, L1] = getRandomCovar(d1);
      [C2, L2] = getRandomCovar(d2);
      [C3, L3] = getRandomCovar(d3);
      X = randn(N, d1) * L1;
      Y = randn(N, d2) * L2;
      Z = randn(N, d3) * L3;
      trueVals(1) = 0; % Conditional Shannon MI
      trueVals(2) = 0.5*( log(det(C2)) - log(det(C1)) - d1 + trace(C2\C1)); ...
          % Conditional KL 
      trueVals(3) = nan; % Conditional Tsallis

    elseif testIdx == 3
      d1 = 1; d2 = 1; d3 = 4; d = d1+d2+d3;
      N = 5000; 
      [C, L] = getRandomCovar(d);
      CX = C(1:d1, 1:d1);
      CY = C( (d1+1):(d1+d2), (d1+1):(d1+d2));
      CZ = C( (d1+d2+1):end, (d1+d2+1):end );
      XZidxs = [1:d1, (d1+d2+1):d]; CXZ = C(XZidxs, XZidxs);
      YZidxs = [(d1+1):d]; CYZ = C(YZidxs, YZidxs);
      XYZ = randn(N, d) * L;
      X = XYZ(:, 1:d1);
      Y = XYZ(:, (d1+1):(d1+d2));
      Z = XYZ(:, (d1+d2+1):end);
      trueVals(1) = 0.5*( log(det(CXZ)) + log(det(CYZ)) - ...
        log(det(C)) - log(det(CZ)) ); % Conditional Shannon MI
      trueVals(2) = 0.5*( log(det(CYZ)) - log(det(CXZ)) -(d1+d3) + ...
        trace(CYZ\CXZ)); % Conditional KL
      trueVals(3) = nan; % Conditional Tsallis

    end

    % Now do the tests
    doTests(X, Y, Z, functionals, trueVals, tests{testIdx}, params, functionalParams);

  end

end


function doTests(X, Y, Z, functionals, trueVals, test, params, functionalParams)

  fprintf('\n%s\n====================================================\n', test);

  for k = 1:numel(functionals)

    functional = functionals{k};
    trueVal = trueVals(k);
    if isnan(trueVal), trueValStr = 'unknown';
    else trueValStr = sprintf('%.5f', trueVal);
    end
    fprintf('  %s: Truth: %s\n', functional, trueValStr);

    switch functional
      case 'condShannonMI'
        [estDS, asympAnalysis] = condShannonMI(X, Y, Z, functionalParams, ...
          params);

      case 'condKLDiv'
        [estDS, asympAnalysis] = condKLDivergence(X, Y, Z, [], functionalParams, ...
          params);

      case 'condTsallisDiv'
        [estDS, asympAnalysis] = condTsallisDivergence(X, Y, Z, [], ...
          functionalParams, params);
      
    end

    % Now compute the errors
    errDS = abs(trueVal - estDS);
    fprintf('    EstimDS : %.4f,  ErrDS : %.4f, CI: %s\n\n', ...
      estDS, errDS, mat2str(asympAnalysis.confInterval));

  end

end


% Function to obtain a random covariance matrix and its Cholesky
% Decomposition
function [K, L] = getRandomCovar(m, n)
  if ~exist('n', 'var')
    n = m;
  end
  A = rand(m, n);
  K = A * A';
  L = chol(K);
end

