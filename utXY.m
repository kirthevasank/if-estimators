function utXY
% Unit tests for all one distro functionals

  close all;
  clear all;

  functionals = {'hellingerDiv', 'tsallisDiv', 'chiSqDiv', 'renyiDiv', 'klDiv'};
  tests = {'1D-UnifUnif', '1D-UnifConv', 'Gaussian'};
  functionalParams = struct;
  functionalParams.alpha = 0.8; % for the Alpha divergences
  params = struct;
  params.alpha = 0.05;
  params.doAsympAnalysis = true;
  params.doBoundaryCorrection = true;
  params.estLowerBound = 1e-5;
  params.kdePickMethod = 'silverman';
  params.kdePickMethod = 'cv';

  % Test 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for testIdx = 1:numel(tests)

    % First generate the data
    if testIdx == 1
      params.doAsympAnalysis = true;
      N1 = 1000; X = rand(N1, 1);
      N2 = 1200; Y = rand(N2, 1);
      trueVals = zeros(5, 1);

    elseif testIdx == 2
      params.doAsympAnalysis = true;
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
      trueVals(1) = 1-0.921261; % Hellinger Divergence
      trueVals(2) = (1-0.9522)/0.2; % Tsallis Divergence
      trueVals(3) = 0.55653; % Chi-Sq Divergence
      trueVals(4) = -1/0.2*log(0.9522); % Tsallis Divergence
      trueVals(5) = 0.290857;

    elseif testIdx == 3
      params.doBoundaryCorrection = false;
      d = 5;
      N1 = 4000; N2 = N1;
      X = randn(N1, d);
      Y = 1 + randn(N2, d);
      trueVals(1) = 1 - exp(-d/8); % Hellinger Divergence
      trueVals(2) = nan; % Tsallis Divergence
      trueVals(3) = nan; % Chi-sq Divergence
      trueVals(4) = functionalParams.alpha * d/2;
      trueVals(5) = d/2;

    end

    % Now do the tests
    doTests(X, Y, functionals, trueVals, tests{testIdx}, params, functionalParams);

  end

end


function doTests(X, Y, functionals, trueVals, test, params, functionalParams)

  fprintf('%s\n====================================================\n', test);

  for k = 1:numel(functionals)

    functional = functionals{k};
    trueVal = trueVals(k);
    fprintf('  %s: Truth: %0.5f\n', functional, trueVal);

    switch functional
      case 'hellingerDiv'
        [estDS, asympAnalysis] = hellingerDivergence(X, Y, functionalParams, ...
          params);

      case 'tsallisDiv'
        [estDS, asympAnalysis] = tsallisDivergence(X, Y, functionalParams, ...
          params);

      case 'chiSqDiv'
        [estDS, asympAnalysis] = chiSqDivergence(X, Y, functionalParams, ...
          params);

      case 'renyiDiv'
        [estDS, asympAnalysis] = renyiDivergence(X, Y, functionalParams, ...
          params);

      case 'klDiv'
        [estDS, asympAnalysis] = klDivergence(X, Y, functionalParams, ...
          params);
      
    end

    % Now compute the errors
    errDS = abs(trueVal - estDS);
    fprintf('    EstimDS : %.4f,  ErrDS : %.4f, CI: %s\n\n', ...
      estDS, errDS, mat2str(asympAnalysis.confInterval));
%     asympAnalysis,

  end

end


% function printResults(trueVal, estDS, asympAnalysis)
% 
%   errDS = abs(trueVal - estDS);
%   fprintf('    EstimDS : %.4f,  ErrDS : %.4f, CI: %s\n\n', ...
%     estDS, errDS, mat2str(asympAnalysis.confInterval));
%   asympAnalysis,
% 
% end

