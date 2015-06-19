function ut_oneDistro
% Unit tests for all one distro functionals

  addpath ../kde
  close all;
  clear all;

  functionals = {'hellingerDiv'};
  tests = {'1D-UnifUnif', '1D-UnifConv'};
  functionalParams = struct;
  params = struct;
  params.alpha = 0.05;
  params.doAsympAnalysis = true;

  % Test 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for testIdx = 1:numel(tests)

    % First generate the data
    if testIdx == 1
      N1 = 1000; X = rand(N1, 1);
      N2 = 1200; Y = rand(N2, 1);
      trueVals = [0];

    elseif testIdx == 2
      N1 = 1000; X = rand(N1, 1);
      N2 = 1000;
      gamma = 10;
      Z = rand(N2, 1+gamma); B = double(rand(N2, 1) < 0.5);
      Y = B.* Z(:,1) + (1-B).*max(Z(:,2:end), [], 2);
      trueYDensity = @(t) 0.5 + 0.5*gamma* t.^(gamma-1);
      params.estLowerBound = 0.4;

      % Compute the true Entropy
      hellingerAffinityFunc = @(t) sqrt(trueYDensity(t));
      t = linspace(0,1,1000); trueHellinger = 1-mean(hellingerAffinityFunc(t));
      % True vals
      trueVals(1) = 1-0.921261;

    end

    % Now do the tests
    doTests(X, Y, functionals, trueVals, tests{testIdx}, params, functionalParams);

  end

end


function doTests(X, Y, functionals, trueVals, test, params, functionalParams)

  fprintf('%s\n====================================================\n\n', test);

  for k = 1:numel(functionals)

    functional = functionals{k};
    trueVal = trueVals(k);
    fprintf('  %s: Truth: %0.5f\n', functional, trueVal);

    % First Do a Data split version
%     params.dataSplit = true;
    [estDS, asympAnylsis] = estimateTwoDistroFunctionals...
      (X, Y, functional, functionalParams, params);
    % Now compute the errors
    errDS = abs(trueVal - estDS);
    fprintf('    EstimDS : %.4f,  ErrDS : %.4f, CI: %s\n', ...
      estDS, errDS, mat2str(asympAnylsis.confInterval));

  end

end

