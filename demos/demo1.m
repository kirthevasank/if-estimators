function demo1
% Unit tests for all one distro functionals

  close all;
  clear all;
  clc;

  fprintf('\nSome demos on estimating functionals of a single distribution.\n');

  functionals = {'shannonEntropy'};
  tests = {'1D-Uniform', '1D-Conv', '2D-Gaussian'};

  % This is for storing parameters specific to the functional (E.g. alpha for the
  % alpha-divergences.) We will not use this here.
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
      N = 1000;
      X = rand(N, 1);
      trueVals = [0];

    elseif testIdx == 2
      N = 1000;
      gamma = 10;
      Z = rand(N, 1+gamma); B = double(rand(N, 1) < 0.5);
      X = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
      trueDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
      params.estLowerBound = 0.4;

      % Compute the true Entropy
      entropyFunc = @(t) trueDensity(t) .* log( trueDensity(t) );
      t = linspace(0,1,1000); trueEntropy = -mean(entropyFunc(t));
      
      % True vals
      trueVals(1) = trueEntropy;

    elseif testIdx == 3
      d = 5;
      N = 10000;
      X = randn(N, d);
      trueVals(1) = d/2 *(1 + log(2*pi));
      params.estLowerBound = 0;

    end

    % Now do the tests
    doTests(X, functionals, trueVals, tests{testIdx}, params, functionalParams);

  end

end


function doTests(X, functionals, trueVals, test, params, functionalParams)

  fprintf('\n%s\n====================================================\n', test);

  for k = 1:numel(functionals)

    functional = functionals{k};
    trueVal = trueVals(k);
    fprintf('  %s: Truth: %0.5f\n', functional, trueVal);

    if strcmp(functional, 'shannonEntropy')
      [estDS, asympAnalysis] = shannonEntropy(X, functionalParams, params);
    end
    % Now compute the errors
    errDS = abs(trueVal - estDS);
    fprintf('    Estimate : %.4f,  Error : %.4f, CI: %s\n', ...
      estDS, errDS, mat2str(asympAnalysis.confInterval));

  end

end

