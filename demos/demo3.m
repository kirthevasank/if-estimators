function demo3
% Unit tests for functionals on X, Y where they come from a joint distribution
% X, Y

  close all;
  clear all;

  fprintf('\nSome demos on conditional functionals of one distribution.\n'); 
  functionals = {'shannonMI', 'condShannonEntropy'};
  tests = {'1D-UnifUnif', 'Indep-Gaussians', 'Gaussian'};

  % This is for storing parameters specific to the functional (E.g. alpha for the
  % alpha-divergences.) We will not use it here.
  functionalParams = struct;

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
      trueVals = zeros(2, 1);

    elseif testIdx == 2
      d1 = 1; d2 = 3;
      N = 4000; 
      [C1, L1] = getRandomCovar(d1);
      [C2, L2] = getRandomCovar(d2);
      X = randn(N, d1) * L1;
      Y = randn(N, d2) * L2;
      trueVals(1) = 0; % Shannon MI
      trueVals(2) = -d1/2 * (1 + log(2*pi)) - 0.5*log(det(C1));

    elseif testIdx == 3
      d1 = 1; d2 = 3; d = d1+d2;
      N = 4000; 
      [C, L] = getRandomCovar(d);
      CX = C(1:d1, 1:d1);
      CY = C( (d1+1):end, (d1+1):end );
      XY = randn(N, d) * L;
      X = XY(:, 1:d1);
      Y = XY(:, (d1+1):end);
      trueVals(1) = 0.5*( -log(det(C)) + log(det(CX)) + log(det(CY)) );
      trueVals(2) = -d1/2 * (1 + log(2*pi)) + 0.5*( log(det(CY)) - log(det(C)));

    end

    % Now do the tests
    doTests(X, Y, functionals, trueVals, tests{testIdx}, params, functionalParams);

  end

end


function doTests(X, Y, functionals, trueVals, test, params, functionalParams)

  fprintf('\n%s\n====================================================\n', test);

  for k = 1:numel(functionals)

    functional = functionals{k};
    trueVal = trueVals(k);
    fprintf('  %s: Truth: %0.5f\n', functional, trueVal);

    switch functional
      case 'shannonMI'
        [estDS, asympAnalysis] = shannonMI(X, Y, functionalParams, ...
          params);

      case 'condShannonEntropy'
        [estDS, asympAnalysis] = condShannonEntropy(X, Y, functionalParams, ...
          params);
      
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

